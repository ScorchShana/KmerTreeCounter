#ifndef FASTQ_PARSER_HEADER
#define FASTQ_PARSER_HEADER

#include "definition.h"
#include "kmer.h"
#include "GetKmer.h"
#include "NewKmerTree.h"
#include "RingMemoryPool.h"
#include "ConcurrentMemoryPool.h"

#include <cstring>
#include <vector>
#include <array>
#include <memory>
#include <thread>
#include <random>

template <uint32_t N>
class FastqParser final
{

    //static constexpr uint64_t KMER_BATCH_CAPACITY = (KMER_BATCH_SIZE - sizeof(uint64_t) * 2) / sizeof(kmer<N>);
    static constexpr uint64_t KMER_BUFFER_CAPACITY = PARSER_KMER_BUFFER_SIZE / sizeof(kmer<N>);

    // 自旋参数
    static constexpr int YIELD_THRESHOLD = 128;
    static constexpr int MAX_BACKOFF = 128;

    const int k;
    uint64_t total_read_kmer = 0;
    RingMemoryPool<RING_MEMORY_POOL_CAPACITY> *ring_memory_pool_ptr;
    KmerTree<N> *tree_ptr;
    ConcurrentMemoryPool *memory_pool_;
    GetKmer<N> get_kmer;
    node<N> *local_root_nodes;
    node<N> *tree_root_nodes;
    uint64_t start_root_node_index = 0;

    std::array<uint32_t, 1ULL << (2 * ROOT_BASES)> local_block_prefix_counts{};
    std::array<uint32_t, 1ULL << (2 * ROOT_BASES)> local_block_prefix_sums{};
    std::array<kmer<N>, KMER_BUFFER_CAPACITY> local_block_for_copy{};
    std::array<kmer<N>, KMER_BUFFER_CAPACITY> kmer_buffer{};
    uint64_t kmer_buffer_count = 0;

public:
    explicit FastqParser(const int in_k,
                         RingMemoryPool<RING_MEMORY_POOL_CAPACITY> *in_ring_memory_pool_ptr,
                         KmerTree<N> *in_tree_ptr,
                         ConcurrentMemoryPool *in_memory_pool_ptr)
        : k(in_k),
          ring_memory_pool_ptr(in_ring_memory_pool_ptr), tree_ptr(in_tree_ptr), memory_pool_(in_memory_pool_ptr), get_kmer(k)
    {

        local_root_nodes = new node<N>[1ULL << (2 * ROOT_BASES)]();
        tree_root_nodes = tree_ptr->root_nodes;

        kmer_buffer_count = 0;

        thread_local std::mt19937 generator(std::random_device{}());
        // 分布对象也被复用
        thread_local std::uniform_int_distribution<int> dist(0, (1ULL << (2 * ROOT_BASES)) - 1);
        start_root_node_index = dist(generator);
    }

    ~FastqParser()
    {
        delete[] local_root_nodes;
    }

    void parse_and_push()
    {
        content_type content;
        bool not_empty = true;

        int backoff_iterations = 1;
        int spin_count = 0;

        while ((!ring_memory_pool_ptr->producer_finished()) || not_empty)
        {

            not_empty = ring_memory_pool_ptr->consumer_try_dequeue(content);
            if (not_empty)
            {
                parse(content.data, content.length);
                ring_memory_pool_ptr->consumer_enqueue(content.data);
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
                    backoff_iterations= std::min(backoff_iterations * 2, MAX_BACKOFF);
                    spin_count++;
                }
                else
                {
                    // 向操作系统让出：如果我们自旋太久，持锁线程
                    // 可能已被抢占。让操作系统调度其他线程。
                    std::this_thread::yield();
                }
            }
        }

        // 清理剩下的kmer
        calculate_block_prefix_counts();
        push_kmers_into_local_block_for_copy();
        flush_block_to_tree();

        flush_local_root_nodes_to_tree();
    }

    // 从 FASTQ 数据块解析 k-mer
    void parse(char const *data_ptr, const uint64_t length)
    {
        get_kmer.clear();
        for (uint64_t now_pos = 0; now_pos < length; now_pos++)
        {
            const char c = data_ptr[now_pos];
            if (get_kmer.get_next_one(c)) [[likely]] // 解析出一条完整的序列
            {
                get_canonical_kmer(); // 转换为规范形式的 k-mer
            }
        }
    }

    inline uint64_t get_total_read_kmer() const noexcept
    {
        return total_read_kmer;
    }

private:
    constexpr uint64_t get_root_prefix(const kmer<N> &k_mer) const
    {
        constexpr uint32_t shift_bits = 64 - (ROOT_BASES * 2);
        return k_mer.data[0] >> shift_bits;
    }

    void calculate_block_prefix_counts()
    {
        /*for(uint64_t i=0;i<(1ULL << (2 * ROOT_BASES));++i)
        {
            root_counts[i].fetch_add(local_block_prefix_counts[i], std::memory_order_relaxed);
        }*/

        local_block_prefix_sums[0] = 0;
        for (uint64_t i = 1; i < local_block_prefix_counts.size(); ++i)
        {
            local_block_prefix_sums[i] = local_block_prefix_sums[i - 1] + local_block_prefix_counts[i - 1];
        }
    }

    void push_kmers_into_local_block_for_copy()
    {
        total_read_kmer += kmer_buffer_count;
        for (uint64_t index = 0; index < kmer_buffer_count; index++)
        {
            const uint64_t prefix = get_root_prefix(kmer_buffer[index]);
            const uint64_t pos = local_block_prefix_sums[prefix];
            local_block_for_copy[pos] = kmer_buffer[index];
            local_block_prefix_sums[prefix]++;
        }
    }

    void get_canonical_kmer() noexcept
    {

        kmer<N> &seq_kmer = get_kmer.seq_kmer;
        kmer<N> &rev_kmer = get_kmer.rev_kmer;

        uint64_t mask = -(seq_kmer < rev_kmer); // 无分支掩码

        

        kmer<N> &canonical_kmer = kmer_buffer[kmer_buffer_count];

        if constexpr (N == 1)
        {
            canonical_kmer.data[0] = (seq_kmer.data[0] & mask) | (rev_kmer.data[0] & (~mask));
        }
        else if constexpr (N == 2)
        {
            canonical_kmer.data[0] = (seq_kmer.data[0] & mask) | (rev_kmer.data[0] & (~mask));
            canonical_kmer.data[1] = (seq_kmer.data[1] & mask) | (rev_kmer.data[1] & (~mask));
        }
        else if constexpr (N == 4)
        {
            canonical_kmer.data[0] = (seq_kmer.data[0] & mask) | (rev_kmer.data[0] & (~mask));
            canonical_kmer.data[1] = (seq_kmer.data[1] & mask) | (rev_kmer.data[1] & (~mask));
            canonical_kmer.data[2] = (seq_kmer.data[2] & mask) | (rev_kmer.data[2] & (~mask));
            canonical_kmer.data[3] = (seq_kmer.data[3] & mask) | (rev_kmer.data[3] & (~mask));
        }
        else
        {
            for (uint32_t i = 0; i < N; ++i)
            {
                canonical_kmer.data[i] = (seq_kmer.data[i] & mask) | (rev_kmer.data[i] & (~mask));
            }
        }

        const uint64_t prefix = canonical_kmer.data[0] >> (64 - (ROOT_BASES * 2));
        local_block_prefix_counts[prefix]++;

        kmer_buffer_count++;

        if (kmer_buffer_count >= KMER_BUFFER_CAPACITY)
        {
            calculate_block_prefix_counts();
            push_kmers_into_local_block_for_copy();
            // 对缓存中的所有的 k-mer 进行多层检查和分流（导出或进树）
            flush_block_to_tree();
        }
    }

    inline void flush_block_to_tree() noexcept
    {
        RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY> *export_pool = tree_ptr->get_export_ring_pool();
        uint64_t high_offset = 0;

        char *raw_block_ptr = nullptr;
        ExportBlock<N> *export_block_ptr = nullptr;
        uint64_t export_kmer_count = 0;

        uint64_t read_offset = 0;
        for (uint64_t prefix = 0; prefix < local_block_prefix_counts.size(); ++prefix)
        {
            const uint32_t prefix_count = local_block_prefix_counts[prefix];
            if (prefix_count == 0)
            {
                continue;
            }

            local_block_prefix_counts[prefix] = 0;
            ConcurrentDoubleBloomFilter<N> *bloom_filter = tree_ptr->get_root_bloom_filter(prefix);

            for (uint64_t i = 0; i < prefix_count; ++i)
            {
                const kmer<N> &val = local_block_for_copy[read_offset + i];
                if (bloom_filter->insert(val))
                {
                    if (export_block_ptr == nullptr)
                    {
                        export_pool->producer_dequeue(raw_block_ptr);
                        export_block_ptr = reinterpret_cast<ExportBlock<N> *>(raw_block_ptr);
                        export_kmer_count = 0;
                    }

                    export_block_ptr->k_mers[export_kmer_count++] = val;
                    if (export_kmer_count >= export_block_ptr->k_mers.size())
                    {
                        export_pool->producer_enqueue({raw_block_ptr, export_kmer_count});
                        export_block_ptr = nullptr;
                        raw_block_ptr = nullptr;
                        export_kmer_count = 0;
                    }
                }
                else
                {
                    kmer_buffer[high_offset++] = val;
                    local_block_prefix_counts[prefix]++;
                }
            }

            read_offset += prefix_count;
        }

        if (export_block_ptr != nullptr)
        {
            export_pool->producer_enqueue({raw_block_ptr, export_kmer_count});
        }

        if (high_offset > 0)
        {
            tree_ptr->main_add_kmer_block_with_local_root_nodes(kmer_buffer, local_block_prefix_counts, local_root_nodes);
        }
        kmer_buffer_count = 0;
        local_block_prefix_counts.fill(0);
    }

    void flush_local_root_nodes_to_tree()
    {
        tree_ptr->flush_local_root_nodes(local_root_nodes, start_root_node_index);
    }
};

#endif