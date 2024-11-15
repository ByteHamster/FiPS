#pragma once
#include <cstdint>
#include <vector>
#include <fstream>
#include <span>

#include <ips2ra.hpp>
#include <bytehamster/util/MurmurHash64.h>

namespace fips {
template <size_t _lineSize = 256, typename _offsetType = uint16_t, bool _useUpperRank = true>
class FiPS {
    public:
        union CacheLine {
            static constexpr size_t LINE_SIZE = _lineSize;
            using offset_t = _offsetType;
            static constexpr size_t OFFSET_SIZE = 8 * sizeof(offset_t);
            static constexpr size_t PAYLOAD_BITS = LINE_SIZE - OFFSET_SIZE;
            static_assert(LINE_SIZE % 64 == 0);
            static_assert(LINE_SIZE > 8 * sizeof(offset_t));

            uint64_t bits[LINE_SIZE / 64] = {0};
            struct {
                [[maybe_unused]] offset_t padding[PAYLOAD_BITS / (8 * sizeof(offset_t))];
                offset_t offset;
            } parts;

            [[nodiscard]] inline bool isSet(size_t idx) const {
                return (bits[idx / 64] & (1ull << (idx % 64))) != 0;
            }

            [[nodiscard]] inline size_t rankInWord(uint64_t word, size_t index) const {
                return std::popcount(word & ((1ull << index) - 1));
            }

            /**
             * Rank query looping over array. Needs random conditional jumps in each iteration.
             */
            [[nodiscard]] inline size_t rankLoop(size_t idx) const {
                size_t popcount = 0;
                for (size_t i = 0; i < idx / 64; i++) {
                    popcount += std::popcount(bits[i]);
                }
                popcount += rankInWord(bits[idx / 64], idx % 64);
                return popcount;
            }

            /**
             * Calculates whole prefix sum even if not necessary and then uses an array index to avoid branches.
             * From Pibiri, Kanda: "Rank/select queries over mutable bitmaps"
             */
            [[nodiscard]] inline size_t rank(size_t idx) const {
                size_t prefix[LINE_SIZE / 64];
                prefix[0] = 0;
                if constexpr (LINE_SIZE > 64) {
                    prefix[1] = std::popcount(bits[0]);
                }
                if constexpr (LINE_SIZE > 128) {
                    prefix[2] = prefix[1] + std::popcount(bits[1]);
                    prefix[3] = prefix[2] + std::popcount(bits[2]);
                }
                if constexpr (LINE_SIZE > 256) {
                    prefix[4] = prefix[3] + std::popcount(bits[3]);
                    prefix[5] = prefix[4] + std::popcount(bits[4]);
                    prefix[6] = prefix[5] + std::popcount(bits[5]);
                    prefix[7] = prefix[6] + std::popcount(bits[6]);
                }
                return prefix[idx / 64] + rankInWord(bits[idx / 64], idx % 64);
            }
        };
    private:
        static constexpr size_t UPPER_RANK_SAMPLING = (size_t(std::numeric_limits<typename CacheLine::offset_t>::max()) + 1) / CacheLine::PAYLOAD_BITS;
        std::vector<CacheLine> bitVector;
        std::vector<size_t> levelBases;
        std::vector<size_t> upperRank;
        size_t levels = 0;
    public:
        explicit FiPS(std::span<const std::string> keys, float gamma = 2.0f) {
            std::vector<uint64_t> hashes;
            hashes.reserve(keys.size());
            for (const std::string &key : keys) {
                hashes.push_back(bytehamster::util::MurmurHash64(key));
            }
            construct(hashes, gamma);
        }

        explicit FiPS(std::span<const uint64_t> keys, float gamma = 2.0f) {
            std::vector<uint64_t> modifiableKeys(keys.begin(), keys.end());
            construct(modifiableKeys, gamma);
        }

        explicit FiPS(std::istream &is) {
            uint64_t TAG;
            is.read(reinterpret_cast<char *>(&TAG), sizeof(TAG));
            assert(TAG == 0xf1b5);
            size_t size;
            is.read(reinterpret_cast<char *>(&size), sizeof(size));
            levelBases.resize(size);
            is.read(reinterpret_cast<char *>(levelBases.data()), size * sizeof(size_t));
            is.read(reinterpret_cast<char *>(&size), sizeof(size));
            bitVector.resize(size);
            is.read(reinterpret_cast<char *>(bitVector.data()), size * sizeof(CacheLine));
            levels = size - 1;
        }

        void writeTo(std::ostream &os) {
            uint64_t TAG = 0xf1b5;
            os.write(reinterpret_cast<const char *>(&TAG), sizeof(TAG));
            size_t size = levelBases.size();
            os.write(reinterpret_cast<const char *>(&size), sizeof(size));
            os.write(reinterpret_cast<const char *>(levelBases.data()), size * sizeof(size_t));
            size = bitVector.size();
            os.write(reinterpret_cast<const char *>(&size), sizeof(size));
            os.write(reinterpret_cast<const char *>(bitVector.data()), size * sizeof(CacheLine));
        }

        void construct(std::vector<uint64_t> &remainingKeys, float gamma = 2.0f) {
            size_t levelBase = 0;
            size_t currentCacheLineIdx = 0;
            CacheLine currentCacheLine = {};
            size_t prefixSum = 0;
            levelBases.push_back(0);
            if constexpr (_useUpperRank) {
                upperRank.push_back(0);
            }

            size_t level = 0;
            while (!remainingKeys.empty()) {
                size_t n = remainingKeys.size();
                size_t domain = ((size_t(n * gamma) + CacheLine::PAYLOAD_BITS - 1)
                        / CacheLine::PAYLOAD_BITS) * CacheLine::PAYLOAD_BITS;
                bitVector.reserve((levelBase + domain) / CacheLine::PAYLOAD_BITS);

                std::vector<uint64_t> collision;
                collision.reserve(size_t(float(n) * gamma * exp(-gamma)));
                if (level > 0) {
                    for (size_t i = 0; i < n; i++) {
                        remainingKeys[i] = bytehamster::util::remix(remainingKeys[i]);
                    }
                }
                ips2ra::sort(remainingKeys.begin(), remainingKeys.end());

                for (size_t i = 0; i < n; i++) {
                    size_t fingerprint = bytehamster::util::fastrange64(remainingKeys[i], domain) + levelBase;
                    size_t idx = fingerprint / CacheLine::PAYLOAD_BITS;
                    flushCacheLineIfNeeded(currentCacheLine, currentCacheLineIdx, prefixSum, idx);

                    if (i + 1 < n && fingerprint == bytehamster::util::fastrange64(remainingKeys[i + 1], domain) + levelBase) {
                        do {
                            collision.push_back(remainingKeys[i]);
                            i++;
                        } while (i < n && fingerprint == bytehamster::util::fastrange64(remainingKeys[i], domain) + levelBase);
                        i--;
                    } else {
                        size_t idxInCacheLine = fingerprint % CacheLine::PAYLOAD_BITS;
                        currentCacheLine.bits[idxInCacheLine / 64] |= 1ul << (idxInCacheLine % 64);
                        prefixSum++;
                    }
                }
                levelBase += domain;
                flushCacheLineIfNeeded(currentCacheLine, currentCacheLineIdx, prefixSum, levelBase / CacheLine::PAYLOAD_BITS);
                levelBases.push_back(levelBase);
                remainingKeys = std::move(collision);
                level++;
            }
            levels = levelBases.size() - 1;
        }

        void flushCacheLineIfNeeded(CacheLine &currentCacheLine, size_t &currentCacheLineIdx, size_t &prefixSum, size_t targetIdx) {
            while (currentCacheLineIdx < targetIdx) {
                bitVector.push_back(currentCacheLine);
                currentCacheLineIdx++;
                if (currentCacheLineIdx % UPPER_RANK_SAMPLING == 0) {
                    if constexpr (_useUpperRank) {
                        assert(upperRank.size() == currentCacheLineIdx / UPPER_RANK_SAMPLING);
                        upperRank.push_back(upperRank.back() + prefixSum);
                        prefixSum = 0;
                    } else {
                        throw std::runtime_error("Too many keys, enable useUpperRank or increase offset type");
                    }
                }
                currentCacheLine = {};
                assert(prefixSum < std::numeric_limits<typename CacheLine::offset_t>::max());
                currentCacheLine.parts.offset = prefixSum;
            }
        }

        [[nodiscard]] size_t getBits() const {
            return 8 * (levelBases.size() * sizeof(size_t)
                          + upperRank.size() * sizeof(size_t)
                          + bitVector.size() * sizeof(CacheLine)
                          + sizeof(*this));
        }

        [[nodiscard]] size_t operator()(const std::string &key) const {
            return this->operator()(bytehamster::util::MurmurHash64(key));
        }

        [[nodiscard]] size_t operator()(uint64_t key) const {
            size_t level = 0;
            do {
                const size_t levelBase = levelBases[level];
                const size_t levelSize = levelBases[level + 1] - levelBase;
                const size_t fingerprint = bytehamster::util::fastrange64(key, levelSize) + levelBase;
                const size_t idx = fingerprint / CacheLine::PAYLOAD_BITS;
                const size_t idxInCacheLine = fingerprint % CacheLine::PAYLOAD_BITS;
                const CacheLine &cacheLine = bitVector[idx];
                if (cacheLine.isSet(idxInCacheLine)) {
                    size_t result = cacheLine.parts.offset + cacheLine.rank(idxInCacheLine);
                    if constexpr (_useUpperRank) {
                        result += upperRank[idx / UPPER_RANK_SAMPLING];
                    }
                    return result;
                }
                level++;
                key = bytehamster::util::remix(key);
            } while (level < levels);
            return -1;
        }
};
} // namespace fips
