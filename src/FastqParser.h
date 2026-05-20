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
#include <immintrin.h>


/*
template <uint32_t N>
class FastqParser
{

    //static constexpr uint64_t KMER_BATCH_CAPACITY = (KMER_BATCH_SIZE - sizeof(uint64_t) * 2) / sizeof(kmer<N>);
    static constexpr uint64_t KMER_BUFFER_CAPACITY = PARSER_KMER_BUFFER_SIZE / sizeof(kmer<N>);

    // 自旋参数
    static constexpr int YIELD_THRESHOLD = 128;
    static constexpr int MAX_BACKOFF = 128;

    const int k;
    uint64_t total_read_kmer = 0;
    RingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY> *ring_memory_pool_ptr;
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
                         RingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY> *in_ring_memory_pool_ptr,
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

    void flush_block_to_tree() noexcept
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
*/

template <uint32_t N>
class FastqParser
{
    // 自旋参数
    static constexpr int YIELD_THRESHOLD = 256;
    static constexpr int MAX_BACKOFF = 64;

    int k_len;

    uint64_t total_read_kmer = 0;
    RingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY> *reader_parser_ring_pool;
    RingMemoryPool<PARSER_CLASSIFIER_RING_MEMORY_POOL_CAPACITY> *parser_classifier_ring_pool;
    KmerTree<N> *tree;

    GetKmer<N> get_kmer;

public:
    explicit FastqParser(uint32_t in_k_len, RingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY> *in_reader_parser_ring_pool,
                         RingMemoryPool<PARSER_CLASSIFIER_RING_MEMORY_POOL_CAPACITY> *in_parser_classifier_ring_pool, KmerTree<N> *in_tree)
        : k_len(in_k_len), reader_parser_ring_pool(in_reader_parser_ring_pool), parser_classifier_ring_pool(in_parser_classifier_ring_pool), tree(in_tree), get_kmer(in_k_len)
    {
    }

    void parse_and_push()
    {
        content_type reader_parser_content;
        bool not_empty = true;

        int backoff_iterations = 1;
        int spin_count = 0;
        //当队列确实为空且生产者已结束时才退出循环
         while (true)
        {
            if (reader_parser_ring_pool->consumer_try_dequeue(reader_parser_content))
            {
                parse(reader_parser_content.data, reader_parser_content.length);
                reader_parser_ring_pool->consumer_enqueue(reader_parser_content.data);
                backoff_iterations = 1;
                spin_count = 0;
            }else if (reader_parser_ring_pool->producer_finished())
            {
                break;
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

    uint64_t get_total_read_kmer() const noexcept
    {
        return total_read_kmer;
    }

private:
    #if defined(__AVX2__)
        const __m256i MASK_UPPER = _mm256_set1_epi8(0xDF);
        const __m256i A = _mm256_set1_epi8('A');
        const __m256i C = _mm256_set1_epi8('C');
        const __m256i G = _mm256_set1_epi8('G');
        const __m256i T = _mm256_set1_epi8('T');
        const __m256i U = _mm256_set1_epi8('U');
        
        
    #elif defined(__SSE4_2__)
        const __m128i MASK_UPPER = _mm_set1_epi8(0xDF);
        const __m128i A = _mm_set1_epi8('A');
        const __m128i C = _mm_set1_epi8('C');
        const __m128i G = _mm_set1_epi8('G');
        const __m128i T = _mm_set1_epi8('T');
        const __m128i U = _mm_set1_epi8('U');
    #else
    #endif

    void parse(char const *data_ptr, const uint64_t length)
    {
        content_type parser_classifier_content;
        parser_classifier_ring_pool->producer_dequeue(parser_classifier_content.data);
        uint64_t content_kmer_count = 0;
        // 获取 k-mer，把 k-mer写入到parser_classifier_content的data里面，length是k-mer的数量，data最大是64KB
        // 解析出 k-mer 后，total_read_kmer += k-mer数量;
        // 解析出 k-mer 后，写入 parser_classifier_content.data，并设置 parser_classifier_content.length
        kmer<N> *kmer_buffer = reinterpret_cast<kmer<N> *>(parser_classifier_content.data);

    #if defined(__AVX2__)
        // 一次处理32字节
        const char* ptr = data_ptr;
        const char* end = data_ptr + length;
        alignas(32) uint8_t codes[32];

        while (ptr < end) {
            size_t chunk = std::min<size_t>(end - ptr, 32);
            __m256i input = _mm256_loadu_si256((__m256i*)ptr);

            __m256i upper = _mm256_and_si256(input, MASK_UPPER);
            __m256i uA = _mm256_cmpeq_epi8(upper, A);
            __m256i uC = _mm256_cmpeq_epi8(upper, C);
            __m256i uG = _mm256_cmpeq_epi8(upper, G);
            __m256i uTU = _mm256_or_si256(_mm256_cmpeq_epi8(upper, T), _mm256_cmpeq_epi8(upper, U));

            __m256i bit0_mask = _mm256_or_si256(uC, uTU);
            __m256i bit1_mask = _mm256_or_si256(uG, uTU);
            
            __m256i bit0 = _mm256_and_si256(bit0_mask, _mm256_set1_epi8(1));
            __m256i bit1 = _mm256_and_si256(bit1_mask, _mm256_set1_epi8(2));
            __m256i encode = _mm256_or_si256(bit1, bit0);
        
            _mm256_store_si256((__m256i*)codes, encode);

            __m256i base_mask = _mm256_or_si256(
                _mm256_or_si256(uA, uC),
                bit1_mask
            );

            uint32_t base_bits = _mm256_movemask_epi8(base_mask) & ((1ULL << chunk) - 1);

            uint32_t all_valid_bits = base_bits;
            uint32_t invalid_bits = (~all_valid_bits) & ((1ULL << chunk) - 1);

            // 合并所有位（有效位 + 无效位），遍历处理
            uint32_t bits = all_valid_bits | invalid_bits;
            while (bits) {
                int idx = __builtin_ctz(bits);
                uint32_t bit = 1u << idx;

                if (base_bits & bit) {
                    // 找到一段连续碱基运行；计算其长度
                    uint32_t run_bits = base_bits >> idx;
                    int run_len = 1;
                    // 限定最长 16 个碱基
                    while ((idx + run_len < 32) && (run_len < 16) && (run_bits & (1u << run_len))) {
                        ++run_len;
                    }

                    // 将此运行中的所有碱基编码打包为一个 uint32_t（查表 A→00,C→01,G→10,T/U→11）
                    uint32_t packed = 0;
                    for (int i = 0; i < run_len; ++i) {
                        uint8_t byte_val = static_cast<uint8_t>(ptr[idx + i]);
                        packed = (packed << 2) | codes[idx + i];
                    }

                    // 批次插入前确保缓冲区有足够空间（一次最多写入 run_len 个 k-mer）
                    constexpr uint64_t KMER_BUF_CAP = PARSER_CLASSIFIER_RING_MEMORY_POOL_BLOCK_SIZE / sizeof(kmer<N>);
                    if (content_kmer_count + static_cast<uint64_t>(run_len) > KMER_BUF_CAP) [[unlikely]] {
                        parser_classifier_content.length = content_kmer_count;
                        parser_classifier_ring_pool->producer_enqueue(parser_classifier_content);
                        total_read_kmer += content_kmer_count;
                        content_kmer_count = 0;
                        parser_classifier_ring_pool->producer_dequeue(parser_classifier_content.data);
                        kmer_buffer = reinterpret_cast<kmer<N>*>(parser_classifier_content.data);
                    }

                    // 批次插入：一次更新 k-mer 状态，同时拿到新完成的 k-mer
                    uint32_t new_kmers = get_kmer.batch_insert(packed, run_len,
                        &kmer_buffer[content_kmer_count]);
                    content_kmer_count += new_kmers;

                    // 清除已处理的位
                    uint32_t clear_mask = ((1u << run_len) - 1) << idx;
                    bits &= ~clear_mask;
                } else {
                    // 遇到无效字符或换行符，清空 k-mer 状态
                    get_kmer.clear();
                    bits &= bits - 1;
                }
            }
            ptr += chunk;
        }
    #elif defined(__SSE4_2__)
        // 一次处理16字节
        const char* ptr = data_ptr;
        const char* end = data_ptr + length;
        alignas(16) uint8_t codes[16];

        while (ptr < end) {
             size_t chunk = std::min<size_t>(end - ptr, 16);
         __m128i input = _mm_loadu_si128((__m128i*)ptr);

             __m128i upper = _mm_and_si128(input, MASK_UPPER);
            __m128i uA = _mm_cmpeq_epi8(upper, A);
            __m128i uC = _mm_cmpeq_epi8(upper, C);
            __m128i uG = _mm_cmpeq_epi8(upper, G);
            __m128i uTU = _mm_or_si128(_mm_cmpeq_epi8(upper, T), _mm_cmpeq_epi8(upper, U));

            __m128i bit0_mask = _mm_or_si128(uC, uTU);            
            __m128i bit1_mask = _mm_or_si128(uG, uTU);           
            __m128i bit0 = _mm_and_si128(bit0_mask, _mm_set1_epi8(1));
            __m128i bit1 = _mm_and_si128(bit1_mask, _mm_set1_epi8(2));
            __m128i encode = _mm_or_si128(bit1, bit0);            

            __m128i base_mask = _mm_or_si128(_mm_or_si128(uA, uC), bit1_mask);
            uint16_t base_bits = static_cast<uint16_t>(_mm_movemask_epi8(base_mask) & ((1u << chunk) - 1));
            
            _mm_store_si128((__m128i*)codes, encode);

            // 合并有效位与无效位
            uint16_t all_valid = base_bits;
            uint16_t invalid = static_cast<uint16_t>((~all_valid) & ((1u << chunk) - 1));
            uint16_t bits = all_valid | invalid;
            while (bits) {
                int idx = __builtin_ctz(bits);
                uint16_t bit = 1u << idx;

                if (base_bits & bit) {
                    uint16_t run_bits = base_bits >> idx;
                    int run_len = 1;
                    while ((idx + run_len < 16) && (run_len < 16) && (run_bits & (1u << run_len))) {
                        ++run_len;
                    }

                    uint32_t packed = 0;
                    for (int i = 0; i < run_len; ++i) {
                        uint8_t byte_val = static_cast<uint8_t>(ptr[idx + i]);
                        packed = (packed << 2) | codes[idx + i];
                    }

                    constexpr uint64_t KMER_BUF_CAP = PARSER_CLASSIFIER_RING_MEMORY_POOL_BLOCK_SIZE / sizeof(kmer<N>);
                    if (content_kmer_count + static_cast<uint64_t>(run_len) > KMER_BUF_CAP) [[unlikely]] {
                        parser_classifier_content.length = content_kmer_count;
                        parser_classifier_ring_pool->producer_enqueue(parser_classifier_content);
                        total_read_kmer += content_kmer_count;
                        content_kmer_count = 0;
                        parser_classifier_ring_pool->producer_dequeue(parser_classifier_content.data);
                        kmer_buffer = reinterpret_cast<kmer<N>*>(parser_classifier_content.data);
                    }

                    uint32_t new_kmers = get_kmer.batch_insert(packed, run_len,
                        &kmer_buffer[content_kmer_count]);
                    content_kmer_count += new_kmers;

                    uint16_t clear_mask = static_cast<uint16_t>(((1u << run_len) - 1) << idx);
                    bits &= ~clear_mask;
                } else {
                    get_kmer.clear();
                    bits &= bits - 1;
                }
            }
            ptr += chunk;
        }
    #else
        for (int i = 0; i < length; ++i)
        {
            const char c = data_ptr[i];
            if (get_kmer.get_next_one(c))
            {
                kmer_buffer[content_kmer_count++] = get_kmer.canonical_kmer;
                if (content_kmer_count >= (PARSER_CLASSIFIER_RING_MEMORY_POOL_BLOCK_SIZE / sizeof(kmer<N>))) [[unlikely]]
                {
                    parser_classifier_content.length = content_kmer_count;
                    parser_classifier_ring_pool->producer_enqueue(parser_classifier_content);
                    total_read_kmer += content_kmer_count;
                    content_kmer_count = 0;
                    parser_classifier_ring_pool->producer_dequeue(parser_classifier_content.data);
                    kmer_buffer = reinterpret_cast<kmer<N> *>(parser_classifier_content.data);
                }
            }
        }
    #endif
        if (content_kmer_count > 0)
        {
            parser_classifier_content.length = content_kmer_count;
            parser_classifier_ring_pool->producer_enqueue(parser_classifier_content);
            total_read_kmer += content_kmer_count;
        }
        get_kmer.clear();
    }
};

    
#endif