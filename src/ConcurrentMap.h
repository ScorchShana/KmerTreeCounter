#ifndef CONCURRENT_MAP_HEADER
#define CONCURRENT_MAP_HEADER

#include "definition.h"
#include "kmer.h"
#include "ConcurrentMemoryPool.h"

#include <array>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <new>
#include <vector>


template <uint32_t N>
struct alignas(CACHE_LINE_SIZE) concurrent_node
{
    static_assert(std::is_trivially_copyable_v<kmer<N>>, "kmer<N> must be trivially copyable");
    kmer<N> k_mer;
    concurrent_node *next;
    std::atomic<uint32_t> count;
};

inline bool check_is_prime(uint64_t num)
{
    for (uint64_t i = 2; i * i <= num; i++)
    {
        if (num % i == 0)
            return false;
    }
    return true;
}

inline uint64_t get_min_prime(uint64_t num)
{
    if (num % 2 == 0)
        num++;
    while (!check_is_prime(num))
        num += 2;
    return num;
}

template <uint32_t N>
class ConcurrentMap
{

    static_assert(KMER_BLOCK_SIZE >= sizeof(concurrent_node<N>), "KMER_BLOCK_SIZE too small for concurrent_node");

    static constexpr size_t align_up(const size_t value, const size_t alignment)
    {
        return (value + alignment - 1) & ~(alignment - 1);
    }

    static constexpr size_t NODE_SLOT_STRIDE = align_up(sizeof(concurrent_node<N>), alignof(concurrent_node<N>));
    static constexpr size_t SLOTS_PER_BLOCK = KMER_BLOCK_SIZE / NODE_SLOT_STRIDE;

    static_assert(SLOTS_PER_BLOCK > 0, "KMER_BLOCK_SIZE too small for aligned concurrent_node slots");

public:
    uint64_t capacity;

    // 按元素数量分配桶，并记录线程数
    explicit ConcurrentMap(ConcurrentMemoryPool *memory_pool, const uint64_t element_num)
        : map_size(0), memory_pool_(memory_pool)
    {
        capacity = get_min_prime(element_num);
        fastmod_multiplier_ = static_cast<uint64_t>((((1ULL << 32) + capacity - 1) / capacity));
        nodes.reset(new std::atomic<concurrent_node<N> *>[capacity]);
        for (uint64_t i = 0; i < capacity; ++i)
        {
            nodes[i].store(nullptr, std::memory_order_relaxed);
        }
    }

    // 析构：释放桶数组
    ~ConcurrentMap() = default;

    // 哈希：将 k-mer 映射到桶索引
    uint64_t hash_func(const kmer<N> &k_mer) const
    {
        uint64_t h = k_mer.data[0];
        for (uint32_t i = 1; i < N; i++)
        {
            h ^= k_mer.data[i];
            h *= 0xc4ceb9fe1a85ec53ULL;
        }
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdULL;
        h ^= h >> 33;
        return h;
    }


    void increment(const kmer<N> &k_mer, const uint32_t &count = 1)
    {

        const uint64_t index = fastmod_reduce(hash_func(k_mer));
        std::atomic<concurrent_node<N> *> &bucket = bucket_head(index);

        concurrent_node<N> *old_chain_first_node = bucket.load(std::memory_order_acquire);
        concurrent_node<N> *last_find_first_node = nullptr;

        for (concurrent_node<N> *cur_node = old_chain_first_node; cur_node != last_find_first_node; cur_node = cur_node->next)
        {
            if (k_mer == cur_node->k_mer)
            {
                cur_node->count.fetch_add(count, std::memory_order_relaxed);
                return;
            }
        }
        last_find_first_node = old_chain_first_node;

        concurrent_node<N> *new_chain_first_node = nullptr;

        new_chain_first_node = allocate_node();

        new_chain_first_node->k_mer = k_mer;
        new_chain_first_node->count.store(count, std::memory_order_relaxed);

        while (true)
        {
            new_chain_first_node->next = old_chain_first_node;
            if (bucket.compare_exchange_strong(old_chain_first_node, new_chain_first_node,
                                               std::memory_order_release, std::memory_order_acquire))
            {
                map_size.fetch_add(1, std::memory_order_relaxed);
                return;
            }

            for (concurrent_node<N> *cur_node = old_chain_first_node; cur_node != last_find_first_node; cur_node = cur_node->next)
            {
                if (k_mer == cur_node->k_mer)
                {
                    cur_node->count.fetch_add(count, std::memory_order_relaxed);
                    release_node(new_chain_first_node);
                    return;
                }
            }
            last_find_first_node = old_chain_first_node;
        }
    }

    // 当前元素数量
    inline uint64_t size() const
    {
        return map_size.load(std::memory_order_relaxed);
    }

private:

    [[nodiscard]] inline uint64_t fastmod_reduce(const uint64_t hash_value) const noexcept
    {
        const uint64_t divisor = capacity;
        if (divisor == 0)
        {
            return 0;
        }

        const uint64_t product = hash_value * fastmod_multiplier_;
        const uint64_t quotient = product >> 32;
        uint64_t remainder = hash_value - quotient * divisor;

        if (remainder >= divisor)
        {
            remainder -= divisor;
            if (remainder >= divisor) [[unlikely]]
            {
                remainder %= divisor;
            }
        }

        return remainder;
    }


    // 获取桶头指针
    inline std::atomic<concurrent_node<N> *> &bucket_head(const uint64_t index) const noexcept
    {
        return nodes[index];
    };

    // 新的 thread-local 状态（顺序分配 + 单槽回收）
    inline static thread_local std::byte* cur_block_ = nullptr;  // 当前 4KB block
    inline static thread_local uint32_t cur_slot_ = 0;           // 当前 block 已用槽位
    inline static thread_local concurrent_node<N>* waste_slot_ = nullptr;  // CAS 失败回收槽

    // 从线程本地 handle 分配节点
    [[nodiscard]] inline concurrent_node<N> *allocate_node()
    {
        // 1. 优先复用 CAS 失败回收的节点
        if (waste_slot_ != nullptr) {
            concurrent_node<N> *node = waste_slot_;
            waste_slot_ = nullptr;
            return node;
        }
        
        // 2. 从 cur_block 顺序分配
        if (cur_block_ == nullptr || cur_slot_ >= SLOTS_PER_BLOCK) {
            cur_block_ = static_cast<std::byte *>(memory_pool_->allocate());
            cur_slot_ = 0;
        }
        
        std::byte *slot = cur_block_ + (cur_slot_ * NODE_SLOT_STRIDE);
        cur_slot_++;
        return reinterpret_cast<concurrent_node<N> *>(slot);
    };

    // 归还节点到线程本地 handle
    inline void release_node( concurrent_node<N> *node) const noexcept
    {
        // 单槽位回收：CAS 竞争失败时暂存，供下次分配复用
        if (waste_slot_ == nullptr) {
            waste_slot_ = node;
        }
        // 如果槽位已被占用，丢弃节点（理论上不会发生，因为 CAS 失败是串行的）
    }
    ConcurrentMemoryPool *memory_pool_ = nullptr;
    uint64_t fastmod_multiplier_;
    std::unique_ptr<std::atomic<concurrent_node<N> *>[]> nodes;
    alignas(64) std::atomic<uint64_t> map_size{};
};

#endif