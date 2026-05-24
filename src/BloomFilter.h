#ifndef BLOOM_FILTER_HEADER
#define BLOOM_FILTER_HEADER

#include "kmer.h"
#include "HashFunction.h"
#include "../include/xxh3.h"
#include "ConcurrentMemoryPool.h"

#include <array>
#include <cstdint>
#include <atomic>

template <uint32_t N>
class ConcurrentDoubleBloomFilter
{
private:
    std::atomic<uint64_t> *filter_bins = nullptr;
    const uint64_t capacity_;
    const uint64_t mod;

public:
    static constexpr uint64_t NUM_HASHES = 3;                          // 使用8个哈希函数
    static constexpr uint64_t BITS_PER_ELEMENT = sizeof(uint64_t) * 8; // 每个元素分配的比特数
    static constexpr uint64_t BITS_MOD = BITS_PER_ELEMENT - 1;         // 用于计算位偏移的掩码

    static constexpr uint64_t SEED_A = 0x9e3779b97f4a7c15ULL; // 黄金比例
    static constexpr uint64_t SEED_B = 0x6a09e667f3bcc909ULL; // √2
    static constexpr uint64_t SEED_C = 0xbb67ae8584caa73bULL; // √3
    static constexpr uint64_t SEED_D = 0x3c6ef372fe94f82bULL; // √5

    explicit ConcurrentDoubleBloomFilter(size_t in_capacity, ConcurrentMemoryPool *memory_pool)
        : capacity_(in_capacity), mod(capacity_ - 1)
    {
        // capacity_ must be power of 2
        filter_bins = reinterpret_cast<std::atomic<uint64_t> *>(memory_pool->allocate_before_init_arenas(sizeof(std::atomic<uint64_t>) * 2 * capacity_));
        // filter_bins = new std::atomic<uint64_t>[2 * capacity_]();
    }

    // ~ConcurrentDoubleBloomFilter()
    // {
    //     delete[] filter_bins;
    // }

    // return whether the element is newly inserted (not exist before)
    bool insert(const kmer<N> &k_mer) noexcept
    {
        // XXH128_hash_t hash_res1 = XXH3_128bits_withSeed(&k_mer, sizeof(k_mer), SEED_A);
        // const uint64_t h1 = hash_res1.low64;
        // const uint64_t h2 = hash_res1.high64 | 1ULL;
        const uint64_t h1 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_A);
        const uint64_t h2 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_B) | 1ULL;
        const uint64_t h3 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_C) | 1ULL;
        const uint64_t block1_idx = 2 * (h1 & mod);
        const uint64_t insert_num1 = calculate_insert_num(h1, h2);

        const uint64_t snapshot1 = filter_bins[block1_idx].load(std::memory_order_relaxed);
        if ((snapshot1 & insert_num1) != insert_num1)
        {
            const uint64_t last_bin1 = filter_bins[block1_idx].fetch_or(insert_num1, std::memory_order_relaxed);
            if ((last_bin1 & insert_num1) != insert_num1)
            {
                return true; // kmer doesn't exist in the first filter
            }
        }

        const uint64_t block2_idx = block1_idx + 1;
        const uint64_t insert_num2 = calculate_insert_num(h1, h3);

        const uint64_t snapshot2 = filter_bins[block2_idx].load(std::memory_order_relaxed);
        if ((snapshot2 & insert_num2) == insert_num2)
        {
            return false; // kmer exists in the second filter, no need to update
        }
        const uint64_t last_bin2 = filter_bins[block2_idx].fetch_or(insert_num2, std::memory_order_relaxed);
        if ((last_bin2 & insert_num2) != insert_num2)
        {
            return true; // kmer doesn't exist in the second filter
        }

        return false; // kmer exists in both filters
    }

    bool exist(const kmer<N> &k_mer) noexcept
    {
        const uint64_t h1 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_A);
        const uint64_t h2 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_B) | 1ULL;
        const uint64_t h3 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_C) | 1ULL;
        const uint64_t block1_idx = 2 * (h1 & mod);
        const uint64_t insert_num1 = calculate_insert_num(h1, h2);

        const uint64_t snapshot1 = filter_bins[block1_idx].load(std::memory_order_relaxed);
        if ((snapshot1 & insert_num1) != insert_num1)
        {
            return false; // kmer doesn't exist in the first filter
        }

        const uint64_t block2_idx = block1_idx + 1;
        const uint64_t insert_num2 = calculate_insert_num(h1, h3);

        const uint64_t snapshot2 = filter_bins[block2_idx].load(std::memory_order_relaxed);
        if ((snapshot2 & insert_num2) != insert_num2)
        {
            return false; // kmer doesn't exist in the second filter
        }

        return true; // kmer exists in both filters
    }

private:
    uint64_t calculate_insert_num(const uint64_t h1, const uint64_t h2)
    {
        uint64_t insert_num = 0;

        for (uint64_t i = 0; i < NUM_HASHES; i++)
        {
            insert_num |= (1ULL << ((h1 + i * h2) & BITS_MOD));
        }
        return insert_num;
    }
};

#endif