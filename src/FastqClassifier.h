#ifndef FASTQCLASSIFIER_HEADER
#define FASTQCLASSIFIER_HEADER

#include "definition.h"
#include "kmer.h"
#include "RingMemoryPool.h"
#include "NewKmerTree.h"
#include "BloomFilter.h"

#include <array>
#include <cstdint>
#include <bitset>

template <uint32_t N>
class FastqClassifier
{

    // 自旋参数
    static constexpr int YIELD_THRESHOLD = 128;
    static constexpr int MAX_BACKOFF = 128;

    int k_len;
    RingMemoryPool<PARSER_CLASSIFIER_RING_MEMORY_POOL_CAPACITY> *parser_classifier_ring_pool;
    KmerTree<N> *tree;

    std::array<node<N>, 1ULL << (2 * ROOT_BASES)> local_root_nodes{};

    std::array<uint32_t, 1ULL << (2 * ROOT_BASES)> local_block_prefix_counts{};
    std::array<uint32_t, 1ULL << (2 * ROOT_BASES)> local_block_prefix_sums{};
    std::array<kmer<N>, PARSER_CLASSIFIER_RING_MEMORY_POOL_BLOCK_SIZE / sizeof(kmer<N>)> local_block_for_copy{};
    std::bitset<PARSER_CLASSIFIER_RING_MEMORY_POOL_BLOCK_SIZE / sizeof(kmer<N>)> local_block_export_bitmap{};

public:
    explicit FastqClassifier(uint32_t in_k_len, RingMemoryPool<PARSER_CLASSIFIER_RING_MEMORY_POOL_CAPACITY> *in_ring_pool, KmerTree<N> *in_tree)
        : k_len(in_k_len),
          parser_classifier_ring_pool(in_ring_pool), tree(in_tree)
    {
    }

    void classify_and_push()
    {
        content_type content;
        bool not_empty = true;

        int backoff_iterations = 1;
        int spin_count = 0;

        while ((!parser_classifier_ring_pool->producer_finished()) || not_empty)
        {

            not_empty = parser_classifier_ring_pool->consumer_try_dequeue(content);
            if (not_empty)
            {
                kmer<N> *kmer_data = reinterpret_cast<kmer<N> *>(content.data);
                const uint64_t kmer_count = content.length; // length 就是 k-mer数量
                process_block(kmer_data, kmer_count);
                parser_classifier_ring_pool->consumer_enqueue(content.data);
                backoff_iterations = 1;
                spin_count = 0;
            }
            else
            {

                if (spin_count < YIELD_THRESHOLD)
                {
                    // 执行 'backoff_iterations' 次暂停指令
                    for (int i = 0; i < backoff_iterations; ++i)
                    {
                        cpu_relax();
                    }
                    // 指数增加退避时间，直到上限
                    backoff_iterations = std::min(backoff_iterations * 2, MAX_BACKOFF);
                    spin_count++;
                }
                else
                {
                    // 向操作系统让出：如果我们自旋太久，持锁线程
                    // 可能已被抢占。让操作系统调度其他线程。
                    backoff_iterations = 1;
                    spin_count = 0;
                    std::this_thread::yield();
                }
            }
        }
    }

private:
    void process_block(kmer<N> *kmer_data, const uint64_t kmer_count)
    {
        calculate_block_prefix_counts(kmer_data, kmer_count);
        push_kmers_into_local_block_for_copy(kmer_data, kmer_count);
        classify_local_block(kmer_count);
    }

    void calculate_block_prefix_counts(kmer<N> *kmer_data, const uint64_t kmer_count)
    {

        memset(local_block_prefix_counts.data(), 0, local_block_prefix_counts.size() * sizeof(uint32_t));

        // 可以考虑simd加速
        for (uint64_t i = 0; i < kmer_count; ++i)
        {
            const uint64_t prefix = get_root_prefix(kmer_data[i]);
            local_block_prefix_counts[prefix]++;
        }

        local_block_prefix_sums[0] = 0;
        for (uint64_t i = 1; i < local_block_prefix_counts.size(); ++i)
        {
            local_block_prefix_sums[i] = local_block_prefix_sums[i - 1] + local_block_prefix_counts[i - 1];
        }
    }

    void push_kmers_into_local_block_for_copy(kmer<N> *kmer_data, const uint64_t kmer_count)
    {
        for (uint64_t index = 0; index < kmer_count; index++)
        {
            const uint64_t prefix = get_root_prefix(kmer_data[index]);
            const uint64_t pos = local_block_prefix_sums[prefix];
            local_block_for_copy[pos] = kmer_data[index];
            local_block_prefix_sums[prefix]++;
        }
    }

    void classify_local_block(const uint64_t kmer_count) noexcept
    {
        RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY> *export_pool = tree->get_export_ring_pool();

        char *raw_block_ptr = nullptr;
        export_pool->producer_dequeue(raw_block_ptr);
        ExportBlock<N> *export_block_ptr = reinterpret_cast<ExportBlock<N> *>(raw_block_ptr);
        uint64_t export_kmer_count = 0;
        uint64_t local_block_count = 0;

        local_block_export_bitmap.reset();

        uint64_t read_offset = 0;
        for (uint64_t prefix = 0; prefix < local_block_prefix_counts.size(); ++prefix)
        {
            const uint32_t prefix_count = local_block_prefix_counts[prefix];
            if (prefix_count == 0)
            {
                continue;
            }

            ConcurrentDoubleBloomFilter<N> *bloom_filter = tree->get_root_bloom_filter(prefix);

            uint32_t prefix_export_count = 0;
            for (uint32_t i = 0; i < prefix_count; ++i)
            {
                const kmer<N> &val = local_block_for_copy[read_offset];
                /*if (bloom_filter->insert(val))
                {
                   export_block_ptr->k_mers[export_kmer_count++] = val;
                }
                else
                {
                    local_block_for_copy[local_block_count] = val;
                    local_block_count++;
                    local_block_prefix_counts[prefix]++;
                }*/
                if (bloom_filter->insert(val))
                {
                    local_block_export_bitmap.set(read_offset);
                    prefix_export_count++;
                }

                read_offset++;
            }
            local_block_prefix_counts[prefix] -= prefix_export_count;
        }

        for (uint64_t i = 0; i < kmer_count; ++i)
        {
            if (local_block_export_bitmap.test(i))
            {
                export_block_ptr->k_mers[export_kmer_count++] = local_block_for_copy[i];
            }
            else
            {
                local_block_for_copy[local_block_count++] = local_block_for_copy[i];
            }
        }

        export_pool->producer_enqueue({raw_block_ptr, export_kmer_count});

        if (local_block_count > 0) [[likely]]
        {
            tree->main_add_kmer_block_with_local_root_nodes(local_block_for_copy, local_block_prefix_counts, local_root_nodes.data());
        }
    }

    uint64_t get_root_prefix(const kmer<N> &k_mer) const
    {
        constexpr uint32_t shift_bits = 64 - (ROOT_BASES * 2);
        return k_mer.data[0] >> shift_bits;
    }
};

#endif // FASTQCLASSIFIER_HEADER