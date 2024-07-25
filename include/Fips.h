#pragma once
#include <cstdint>
#include <vector>
#include <ips2ra.hpp>
#include <MurmurHash64.h>

namespace fips {
class FiPS {
    private:
        struct CacheLine {
            static constexpr size_t PAYLOAD_BITS = 512 - 32;
            union {
                uint64_t bits[8] = {0, 0, 0, 0, 0, 0, 0, 0};
                struct {
                    uint32_t padding[15];
                    uint32_t offset;
                };
            };
        };
        std::vector<CacheLine> bitVector;
        std::vector<size_t> levelBases;
    public:
        explicit FiPS(const std::vector<std::string> &keys, float gamma = 2.0f) {
            std::vector<uint64_t> hashes;
            hashes.reserve(keys.size());
            for (const std::string &key : keys) {
                hashes.push_back(util::MurmurHash64(key));
            }
            construct(hashes, gamma);
        }

        explicit FiPS(const std::vector<uint64_t> &keys, float gamma = 2.0f) {
            std::vector<uint64_t> modifiableKeys = keys;
            construct(modifiableKeys, gamma);
        }

        void construct(std::vector<uint64_t> &remainingKeys, float gamma = 2.0f) {
            size_t levelBase = 0;
            size_t currentCacheLineIdx = 0;
            CacheLine currentCacheLine = {};
            size_t prefixSum = 0;
            levelBases.push_back(0);

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
                        remainingKeys[i] = util::remix(remainingKeys[i]);
                    }
                }
                ips2ra::sort(remainingKeys.begin(), remainingKeys.end());

                for (size_t i = 0; i < n; i++) {
                    size_t fingerprint = util::fastrange64(remainingKeys[i], domain) + levelBase;
                    size_t idx = fingerprint / CacheLine::PAYLOAD_BITS;
                    flushCacheLineIfNeeded(currentCacheLine, currentCacheLineIdx, prefixSum, idx);

                    if (i + 1 < n && fingerprint == util::fastrange64(remainingKeys[i + 1], domain) + levelBase) {
                        do {
                            collision.push_back(remainingKeys[i]);
                            i++;
                        } while (i < n && fingerprint == util::fastrange64(remainingKeys[i], domain) + levelBase);
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
        }

        void flushCacheLineIfNeeded(CacheLine &currentCacheLine, size_t &currentCacheLineIdx, size_t &prefixSum, size_t targetIdx) {
            while (currentCacheLineIdx < targetIdx) {
                bitVector.push_back(currentCacheLine);
                currentCacheLine = {};
                currentCacheLine.offset = prefixSum;
                currentCacheLineIdx++;
            }
        }

        [[nodiscard]] size_t getBits() const {
            return 8 * (levelBases.size() * sizeof(uint64_t) + sizeof(*this) + bitVector.size() * sizeof(CacheLine));
        }

        [[nodiscard]] size_t operator()(const std::string &key) const {
            return this->operator()(util::MurmurHash64(key));
        }

        [[nodiscard]] size_t operator()(uint64_t key) const {
            size_t level = 0;
            while (level < levelBases.size() - 1) {
                size_t fingerprint = util::fastrange64(key, levelBases[level + 1] - levelBases[level]) + levelBases[level];
                size_t idx = fingerprint / CacheLine::PAYLOAD_BITS;
                size_t idxInCacheLine = fingerprint % CacheLine::PAYLOAD_BITS;
                const CacheLine &cacheLine = bitVector[idx];
                if ((cacheLine.bits[idxInCacheLine / 64] & (1ull << (idxInCacheLine % 64))) != 0) {
                    size_t popcount = 0;
                    for (size_t i = 0; i < idxInCacheLine / 64; i++) {
                        popcount += std::popcount(cacheLine.bits[i]);
                    }
                    popcount += std::popcount(cacheLine.bits[idxInCacheLine / 64] & ((1ull << (idxInCacheLine % 64)) - 1));
                    return cacheLine.offset + popcount;
                }
                level++;
                key = util::remix(key);
            }
            return -1;
        }
};
} // namespace fips
