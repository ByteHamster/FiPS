#pragma once
#include <cstdint>
#include <vector>
#include <ips2ra.hpp>

namespace fips {
class FiPS {
    private:
        struct CacheLine {
            uint64_t bits[7] = {0,0,0,0,0,0,0};
            size_t offset = 0;

            void print() {
                printf("%016lx %016lx %016lx %016lx %016lx %016lx %016lx\n",
                       bits[0], bits[1], bits[2], bits[3], bits[4], bits[5], bits[6]);
            }
        };
        std::vector<CacheLine> bitVector;
        std::vector<size_t> levelBases;
    public:
        explicit FiPS (const std::vector<uint64_t> &keys, float gamma = 2.0f) {
            std::vector<uint64_t> remainingKeys = keys;

            size_t levelBase = 0;
            size_t currentCacheLineIdx = 0;
            CacheLine currentCacheLine = {};
            size_t prefixSum = 0;
            levelBases.push_back(0);

            size_t level = 0;
            while (!remainingKeys.empty()) {
                size_t n = remainingKeys.size();
                size_t domain = ((size_t(n * gamma) + (512 - 64) + 1) / (512 - 64)) * (512 - 64);
                bitVector.reserve((levelBase + domain) / (512 - 64));

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
                    size_t idx = fingerprint / (512 - 64);
                    flushCacheLineIfNeeded(currentCacheLine, currentCacheLineIdx, prefixSum, idx);

                    if (i + 1 < n && fingerprint == util::fastrange64(remainingKeys[i + 1], domain) + levelBase) {
                        do {
                            collision.push_back(remainingKeys[i]);
                            i++;
                        } while (i < n && fingerprint == util::fastrange64(remainingKeys[i], domain) + levelBase);
                        i--;
                    } else {
                        size_t idxInCacheLine = fingerprint % (512 - 64);
                        currentCacheLine.bits[idxInCacheLine / 64] |= 1ul << (idxInCacheLine % 64);
                        prefixSum++;
                    }
                }
                levelBase += domain;
                flushCacheLineIfNeeded(currentCacheLine, currentCacheLineIdx, prefixSum, levelBase / (512 - 64));
                levelBases.push_back(levelBase);
                remainingKeys = std::move(collision);
                level++;
            }
        }

        void flushCacheLineIfNeeded(CacheLine &currentCacheLine, size_t &currentCacheLineIdx,
                                    size_t &prefixSum, size_t targetIdx) {
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

        [[nodiscard]] size_t operator()(uint64_t key) const {
            size_t level = 0;
            while (level < levelBases.size() - 1) {
                size_t fingerprint = util::fastrange64(key, levelBases[level + 1] - levelBases[level]) + levelBases[level];
                size_t idx = fingerprint / (512 - 64);
                size_t idxInCacheLine = fingerprint % (512 - 64);
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
