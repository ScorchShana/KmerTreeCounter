#ifndef BLOOM_FILTER_HEADER
#define BLOOM_FILTER_HEADER

#include "kmer.h"
#include "HashFunction.h"
#include "../include/xxh3.h"
#include "ConcurrentMemoryPool.h"

#include <array>
#include <cstdint>
#include <atomic>

#ifdef TEST_MODE
#include <x86intrin.h>

inline uint64_t bloom_read_cycles() noexcept
{
    unsigned aux;
    return __rdtscp(&aux);
}

constexpr uint64_t BLOOM_TIMING_SAMPLE_MASK = 63; // 1/64 sampling

struct BloomTimingStats
{
    uint64_t calls = 0;
    uint64_t sampled_calls = 0;
    uint64_t true_count = 0;
    uint64_t second_layer_attempts = 0;
    uint64_t sampled_second_layer_attempts = 0;
    uint64_t first_fetch_or_count = 0;
    uint64_t second_fetch_or_count = 0;
    uint64_t sampled_first_fetch_or_count = 0;
    uint64_t sampled_second_fetch_or_count = 0;

    uint64_t total_cycles = 0;
    uint64_t hash_cycles = 0;
    uint64_t index_mask_cycles = 0;
    uint64_t first_load_cycles = 0;
    uint64_t first_filter_cycles = 0;
    uint64_t first_fetch_or_cycles = 0;
    uint64_t second_load_cycles = 0;
    uint64_t second_filter_cycles = 0;
    uint64_t second_fetch_or_cycles = 0;
};
#endif

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
        // XXH128_hash_t hash_res1 = XXH3_128bits(&k_mer, sizeof(k_mer));
        // const uint64_t h1 = hash_res1.low64;
        // const uint64_t h2 = hash_res1.high64 | 1ULL;
        uint64_t h1, h2;
        double_hash_func(k_mer, h1, h2);
        // const uint64_t h1 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_A);
        // const uint64_t h2 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_B) | 1ULL;
        const uint64_t block_idx = 2 * (h1 & mod);
        const uint64_t insert_num = calculate_insert_num(h1, h2);
        const uint64_t first_snapshot = filter_bins[block_idx].load(std::memory_order_relaxed);

        uint64_t second_tieme_bit = first_snapshot & insert_num;
        if (second_tieme_bit != insert_num)
        {
            second_tieme_bit = filter_bins[block_idx].fetch_or(insert_num, std::memory_order_relaxed) & insert_num;
        }
        if (second_tieme_bit)
        {
            const uint64_t second_snapshot = filter_bins[block_idx + 1].load(std::memory_order_relaxed);
            if ((second_snapshot & second_tieme_bit) != second_tieme_bit)
            {
                filter_bins[block_idx + 1].fetch_or(second_tieme_bit, std::memory_order_relaxed);
            }

        }
        return second_tieme_bit != insert_num; // kmer not exists in first filters
    }

#ifdef TEST_MODE
    bool insert_timed(const kmer<N> &k_mer, BloomTimingStats &stats) noexcept
    {
        stats.calls++;
        const bool sampled = ((stats.calls & BLOOM_TIMING_SAMPLE_MASK) == 0);

        if (!sampled)
        {
            uint64_t h1, h2;
            double_hash_func(k_mer, h1, h2);
            const uint64_t block_idx = 2 * (h1 & mod);
            const uint64_t insert_num = calculate_insert_num(h1, h2);

            const uint64_t first_snapshot = filter_bins[block_idx].load(std::memory_order_relaxed);
            uint64_t second_tieme_bit = first_snapshot & insert_num;

            if (second_tieme_bit != insert_num)
            {
                stats.first_fetch_or_count++;
                second_tieme_bit = filter_bins[block_idx].fetch_or(insert_num, std::memory_order_relaxed) & insert_num;
            }

            if (second_tieme_bit)
            {
                stats.second_layer_attempts++;

                const uint64_t second_snapshot = filter_bins[block_idx + 1].load(std::memory_order_relaxed);
                if ((second_snapshot & second_tieme_bit) != second_tieme_bit)
                {
                    stats.second_fetch_or_count++;
                    filter_bins[block_idx + 1].fetch_or(second_tieme_bit, std::memory_order_relaxed);
                }
            }

            const bool result = second_tieme_bit != insert_num;
            stats.true_count += result;
            return result;
        }

        stats.sampled_calls++;

        const uint64_t t0 = bloom_read_cycles();
        //XXH128_hash_t hash_res1 = XXH3_128bits(&k_mer, sizeof(k_mer));
        const uint64_t t1 = bloom_read_cycles();

        uint64_t h1, h2;
        double_hash_func(k_mer, h1, h2);
        const uint64_t block_idx = 2 * (h1 & mod);
        const uint64_t insert_num = calculate_insert_num(h1, h2);
        const uint64_t t2 = bloom_read_cycles();

        const uint64_t first_load_start = t2;
        const uint64_t first_snapshot = filter_bins[block_idx].load(std::memory_order_relaxed);
        const uint64_t first_load_end = bloom_read_cycles();

        uint64_t second_tieme_bit = first_snapshot & insert_num;
        uint64_t first_done = first_load_end;
        if (second_tieme_bit != insert_num)
        {
            stats.first_fetch_or_count++;
            stats.sampled_first_fetch_or_count++;

            const uint64_t first_fetch_start = bloom_read_cycles();
            second_tieme_bit = filter_bins[block_idx].fetch_or(insert_num, std::memory_order_relaxed) & insert_num;
            const uint64_t first_fetch_end = bloom_read_cycles();

            stats.first_fetch_or_cycles += first_fetch_end - first_fetch_start;
            first_done = first_fetch_end;
        }

        uint64_t done = first_done;
        if (second_tieme_bit)
        {
            stats.second_layer_attempts++;
            stats.sampled_second_layer_attempts++;

            const uint64_t second_load_start = bloom_read_cycles();
            const uint64_t second_snapshot = filter_bins[block_idx + 1].load(std::memory_order_relaxed);
            const uint64_t second_load_end = bloom_read_cycles();

            uint64_t second_done = second_load_end;
            if ((second_snapshot & second_tieme_bit) != second_tieme_bit)
            {
                stats.second_fetch_or_count++;
                stats.sampled_second_fetch_or_count++;

                const uint64_t second_fetch_start = bloom_read_cycles();
                filter_bins[block_idx + 1].fetch_or(second_tieme_bit, std::memory_order_relaxed);
                const uint64_t second_fetch_end = bloom_read_cycles();

                stats.second_fetch_or_cycles += second_fetch_end - second_fetch_start;
                second_done = second_fetch_end;
            }

            stats.second_load_cycles += second_load_end - second_load_start;
            stats.second_filter_cycles += second_done - second_load_start;
            done = second_done;
        }

        stats.hash_cycles += t2 - t0;
        stats.index_mask_cycles += t2 - t1;
        stats.first_load_cycles += first_load_end - first_load_start;
        stats.first_filter_cycles += first_done - first_load_start;
        stats.total_cycles += done - t0;

        const bool result = second_tieme_bit != insert_num;
        stats.true_count += result;
        return result;
    }
#endif

    bool exist(const kmer<N> &k_mer) noexcept
    {
        uint64_t h1, h2;
        double_hash_func(k_mer, h1, h2);
        // const uint64_t h1 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_A);
        // const uint64_t h2 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_B) | 1ULL;
        const uint64_t block_idx = 2 * (h1 & mod);
        const uint64_t insert_num = calculate_insert_num(h1, h2);

        return (filter_bins[block_idx + 1].load(std::memory_order_relaxed) & insert_num) == insert_num; // kmer exists in second filters
    }

    uint64_t get_capacity() const noexcept
    {
        return capacity_;
    }

    std::atomic<uint64_t> *get_filter_bins() const noexcept
    {
        return filter_bins;
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
/*
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
            return false; // kmer doesn't exist in the second filter
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
};*/

template <uint32_t N>
class CountBloomFilter
{

private:
    uint64_t *filter_bins = nullptr;
    const uint64_t capacity_;
    const uint64_t mod;

public:
    static constexpr uint64_t NUM_HASHES = 3;                          // 使用8个哈希函数
    static constexpr uint64_t BITS_PER_ELEMENT = sizeof(uint64_t) * 8; // 每个元素分配的比特数
    static constexpr uint64_t BITS_MOD = BITS_PER_ELEMENT - 1;         // 用于计算位偏移的掩码

    static constexpr uint64_t SEED_A = 0x9e3779b97f4a7c15ULL; // 黄金比例
    static constexpr uint64_t SEED_B = 0x6a09e667f3bcc909ULL; // √2

    explicit CountBloomFilter(size_t in_capacity)
        : capacity_(in_capacity), mod(capacity_ - 1)
    {
        // capacity_ must be power of 2
        filter_bins = new uint64_t[2 * capacity_]();
        // filter_bins = new uint64_t[2 * capacity_]();
    }

    ~CountBloomFilter()
    {
        delete[] filter_bins;
    }

    bool insert(const kmer<N> &k_mer) noexcept
    {
        XXH128_hash_t hash_res1 = XXH3_128bits_withSeed(&k_mer, sizeof(k_mer), SEED_A);
        const uint64_t h1 = hash_res1.low64;
        const uint64_t h2 = hash_res1.high64 | 1ULL;
        // const uint64_t h1 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_A);
        // const uint64_t h2 = XXH3_64bits_withSeed(k_mer.data.data(), sizeof(kmer<N>), SEED_B) | 1ULL;

        const uint64_t block_idx = 2 * (h1 & mod);

        const uint64_t insert_num = calculate_insert_num(h1, h2);
        const uint64_t second_tieme_bit = filter_bins[block_idx] & insert_num;
        filter_bins[block_idx] |= insert_num;
        filter_bins[block_idx + 1] |= second_tieme_bit;
        return second_tieme_bit != insert_num; // kmer not exists in first filters
    }

    bool exist(const kmer<N> &k_mer) noexcept
    {
        XXH128_hash_t hash_res1 = XXH3_128bits(&k_mer, sizeof(k_mer));
        const uint64_t h1 = hash_res1.low64;
        const uint64_t h2 = hash_res1.high64 | 1ULL;

        const uint64_t block_idx = 2 * (h1 & mod);
        const uint64_t insert_num = calculate_insert_num(h1, h2);

        return (filter_bins[block_idx + 1] & insert_num) == insert_num; // kmer exists in second filters
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
