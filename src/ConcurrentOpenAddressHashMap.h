#ifndef CONCURRENT_OPEN_ADDRESS_HASH_MAP_HEADER
#define CONCURRENT_OPEN_ADDRESS_HASH_MAP_HEADER

#include "definition.h"
#include "ConcurrentMemoryPool.h"
#include "kmer.h"
#include "../include/rapidhash.h"

#include <atomic>
#include <algorithm>
#include <cstdint>
#include <cstring>

template <uint32_t N>
class ConcurrentOpenAddressHashMap
{

    static constexpr size_t align_up(const size_t value, const size_t alignment)
    {
        return (value + alignment - 1) & ~(alignment - 1);
    }

    constexpr static uint8_t EMPTY = 0x00;
    constexpr static uint8_t INSERTING = 0x01;
    constexpr static double LOAD_FACTOR = 0.8;

    static constexpr std::size_t MAP_SIZE = align_up(sizeof(ConcurrentOpenAddressHashMap<N>), CACHE_LINE_SIZE);
    static constexpr std::size_t MAP_NUM_PER_BLOCK = KMER_BLOCK_SIZE / MAP_SIZE;

    std::atomic<int64_t> size;
    uint64_t capacity;
    std::atomic<uint8_t>* ctrls;
    kmer<N>* keys;
    std::atomic<uint32_t>* counts;
    std::atomic<bool> sealed;
    std::atomic<bool> building_next;
    std::atomic<ConcurrentOpenAddressHashMap<N>*> next_map;

    inline static thread_local std::size_t map_num_in_block = 0;
    inline static thread_local char* map_mem_block = nullptr;
    inline static thread_local ConcurrentMemoryPool* pool = nullptr;

public:

    enum class InsertResult
    {
        INSERTED,
        UPDATED,
        FULL
    };

#ifdef CONCURRENT_OPEN_ADDRESS_HASH_MAP_TESTING
    enum class DebugEvent
    {
        BEFORE_SEALED_CHECK,
        AFTER_SEALED_OBSERVED,
        AFTER_OPEN_CHECK_BEFORE_CTRL_CAS,
        AFTER_SIZE_RESERVED_BEFORE_PUBLISH
    };

    using DebugHook = void (*)(DebugEvent,
        const ConcurrentOpenAddressHashMap<N>*,
        const kmer<N>&);

    inline static std::atomic<DebugHook> debug_hook{ nullptr };

    static void set_debug_hook(DebugHook hook)
    {
        debug_hook.store(hook, std::memory_order_release);
    }

    static void notify_debug_hook(DebugEvent event,
        const ConcurrentOpenAddressHashMap<N>* map,
        const kmer<N>& key)
    {
        DebugHook hook = debug_hook.load(std::memory_order_acquire);
        if (hook != nullptr)
        {
            hook(event, map, key);
        }
    }
#endif

    ConcurrentOpenAddressHashMap(uint64_t capacity)
        : size(0), capacity(capacity), sealed(false), building_next(false), next_map(nullptr)
    {
        char* mem_area = reinterpret_cast<char*>(pool->allocate_large(get_mem_size(capacity)));
        ctrls = reinterpret_cast<std::atomic<uint8_t>*>(mem_area);
        keys = reinterpret_cast<kmer<N>*>(ctrls + capacity);
        counts = reinterpret_cast<std::atomic<uint32_t>*>(keys + capacity);
        for (uint64_t i = 0; i < capacity; ++i)
        {
            new (ctrls + i) std::atomic<uint8_t>(EMPTY);
        }
        for (uint64_t i = 0; i < capacity; ++i)
        {
            new (counts + i) std::atomic<uint32_t>(0);
        }
    }

    static constexpr uint64_t get_mem_size(const uint64_t capacity) {
        return capacity * sizeof(std::atomic<uint8_t>) + capacity * sizeof(kmer<N>) + capacity * sizeof(std::atomic<uint32_t>);
    }

    static uint64_t hash_key(const kmer<N>& key)
    {
        const uint64_t h = rapidhash(&key, sizeof(kmer<N>));
        return h;
    }

    // Extract fingerprint from hash (high 7 bits, most significant bit is 1)
    static uint8_t fingerprint(const uint64_t h)
    {
        uint8_t fp = static_cast<uint8_t>(h >> 57);
        fp |= 0x80; // Set the most significant bit to 1
        return fp;
    }

    uint64_t get_capacity() const
    {
        return capacity;
    }

    ConcurrentOpenAddressHashMap<N>* get_next_map() const
    {
        return next_map.load(std::memory_order_acquire);
    }

    // 遍历当前 map 及所有 next_map segment 中的条目。
    // Visitor 签名为 void(const kmer<N>& key, uint32_t count)。
    // 调用方需保证所有写入已结束（如 final drain 阶段），期间不能并发修改。
    template <typename Visitor>
    void for_each_entry(Visitor&& visitor) const
    {
        const ConcurrentOpenAddressHashMap<N>* map = this;
        while (map != nullptr)
        {
            for (uint64_t i = 0; i < map->capacity; ++i)
            {
                const uint8_t ctrl = map->ctrls[i].load(std::memory_order_acquire);
                if ((ctrl & 0x80U) != 0) [[likely]]
                {
                    visitor(map->keys[i], map->counts[i].load(std::memory_order_relaxed));
                }
            }
            map = map->next_map.load(std::memory_order_acquire);
        }
    }

    static void set_memory_pool(ConcurrentMemoryPool* memory_pool)
    {
        pool = memory_pool;
        map_mem_block = reinterpret_cast<char*>(pool->allocate());
    }

    static char* get_map_metadata_mem()
    {
        if (map_num_in_block >= MAP_NUM_PER_BLOCK)
        {
            map_mem_block = reinterpret_cast<char*>(pool->allocate());
            map_num_in_block = 0;
        }
        char* mem_area = map_mem_block + map_num_in_block * MAP_SIZE;
        ++map_num_in_block;
        return mem_area;
    }

    InsertResult try_increment(
        uint64_t index,
        const kmer<N>& key,
        const uint8_t fp,
        const uint32_t value,
        uint64_t& local_count)
    {
        const int64_t cur_max_size = static_cast<int64_t>(capacity * LOAD_FACTOR);
        const uint64_t mod = capacity - 1;
        uint64_t offset = 0;

        for (;;) {
            uint8_t ctrl = ctrls[index].load(std::memory_order_acquire);
            if (ctrl == EMPTY)
            {
                // Try to insert
#ifdef CONCURRENT_OPEN_ADDRESS_HASH_MAP_TESTING
                notify_debug_hook(DebugEvent::BEFORE_SEALED_CHECK, this, key);
#endif
                if (sealed.load(std::memory_order_acquire))
                {
#ifdef CONCURRENT_OPEN_ADDRESS_HASH_MAP_TESTING
                    notify_debug_hook(DebugEvent::AFTER_SEALED_OBSERVED, this, key);
#endif
                    ctrl = ctrls[index].load(std::memory_order_acquire);
                    if (ctrl == EMPTY)
                    {
                        return InsertResult::FULL; // Indicate that the map is full
                    }
                    else
                    {
                        cpu_relax();
                        continue; // Another thread inserted, retry
                    }
                }

#ifdef CONCURRENT_OPEN_ADDRESS_HASH_MAP_TESTING
                notify_debug_hook(DebugEvent::AFTER_OPEN_CHECK_BEFORE_CTRL_CAS, this, key);
#endif

                uint8_t expected = EMPTY;
                if (ctrls[index].compare_exchange_strong(expected,
                    INSERTING,
                    std::memory_order_relaxed,
                    std::memory_order_relaxed))
                {

                    if (sealed.load(std::memory_order_acquire)) [[unlikely]]
                    {
                        ctrls[index].store(EMPTY, std::memory_order_release); // Revert control to EMPTY
                        return InsertResult::FULL; // Indicate that the map is full
                    }
                    // Successfully marked as inserting
                    const uint64_t old_size = size.fetch_add(1, std::memory_order_acq_rel);

                    if (old_size >= static_cast<uint64_t>(cur_max_size)) [[unlikely]]
                    {
                        size.fetch_sub(1, std::memory_order_acq_rel); // Revert size increment
                        ctrls[index].store(EMPTY, std::memory_order_release); // Revert control to EMPTY
                        sealed.store(true, std::memory_order_release);
                        return InsertResult::FULL; // Indicate that the map is full
                    }
                    else
                    {
#ifdef CONCURRENT_OPEN_ADDRESS_HASH_MAP_TESTING
                        notify_debug_hook(DebugEvent::AFTER_SIZE_RESERVED_BEFORE_PUBLISH, this, key);
#endif
                        keys[index] = key;
                        const uint32_t corrected_value = std::min<uint32_t>(value, count_max);
                        counts[index].store(corrected_value, std::memory_order_release);
                        ctrls[index].store(fp, std::memory_order_release); // Set the fingerprint
                        ++local_count;
                        return InsertResult::INSERTED; // Return true if we need to build next map
                    }
                }
                else
                {
                    cpu_relax(); // Another thread is inserting, pause and retry
                    continue;
                }
            }
            else if (ctrl == fp)
            {
                if (keys[index] == key) [[likely]]
                {
                    const uint32_t cur = counts[index].load(std::memory_order_relaxed);
                    if (cur < count_max)
                    {
                        const uint32_t increment = static_cast<uint32_t>(
                            std::min<uint32_t>(value, count_max));
                        const uint32_t prev = counts[index].fetch_add(increment, std::memory_order_relaxed);
                        if (increment > count_max - std::min(prev, count_max)) [[unlikely]]
                        {
                            counts[index].store(count_max, std::memory_order_relaxed);
                        }
                    }
                    return InsertResult::UPDATED;
                }
            }
            else if (ctrl == INSERTING)
            {
                cpu_relax(); // Another thread is inserting, pause and retry
                continue;
            }

            // probe to next index
            ++offset;
            offset &= mod;
            index = (index + offset) & mod;
        }
    }

    void increment(const kmer<N>& key, const uint32_t& value, uint64_t& local_count) {
        const uint64_t h = hash_key(key);
        const uint8_t fp = fingerprint(h);
        const uint64_t mod = capacity - 1;
        uint64_t index = h & mod;

        ConcurrentOpenAddressHashMap<N>* map_ptr = this;

        for (;;)
        {
            InsertResult res = map_ptr->try_increment(index, key, fp, value, local_count);
            if (res == InsertResult::FULL)
            {
                map_ptr = map_ptr->ensure_next_map();
            }
            else
            {
                return;
            }
        }
    }

#ifdef CONCURRENT_OPEN_ADDRESS_HASH_MAP_TESTING
    template <typename Visitor>
    void debug_visit_entries(Visitor&& visitor) const
    {
        const ConcurrentOpenAddressHashMap<N>* map = this;
        uint64_t segment_index = 0;
        while (map != nullptr)
        {
            for (uint64_t i = 0; i < map->capacity; ++i)
            {
                const uint8_t ctrl = map->ctrls[i].load(std::memory_order_acquire);
                if ((ctrl & 0x80U) != 0)
                {
                    visitor(segment_index,
                        map->keys[i],
                        map->counts[i].load(std::memory_order_relaxed));
                }
            }
            map = map->next_map.load(std::memory_order_acquire);
            ++segment_index;
        }
    }

    template <typename Visitor>
    void debug_visit_segments(Visitor&& visitor) const
    {
        const ConcurrentOpenAddressHashMap<N>* map = this;
        uint64_t segment_index = 0;
        while (map != nullptr)
        {
            visitor(segment_index,
                map->size.load(std::memory_order_acquire),
                map->sealed.load(std::memory_order_acquire),
                map->capacity);
            map = map->next_map.load(std::memory_order_acquire);
            ++segment_index;
        }
    }

    uint64_t debug_inserting_slot_count() const
    {
        uint64_t inserting_count = 0;
        const ConcurrentOpenAddressHashMap<N>* map = this;
        while (map != nullptr)
        {
            for (uint64_t i = 0; i < map->capacity; ++i)
            {
                if (map->ctrls[i].load(std::memory_order_acquire) == INSERTING)
                {
                    ++inserting_count;
                }
            }
            map = map->next_map.load(std::memory_order_acquire);
        }
        return inserting_count;
    }
#endif

    ConcurrentOpenAddressHashMap<N>* ensure_next_map() {
        ConcurrentOpenAddressHashMap<N>* next_map_ptr = next_map.load(std::memory_order_acquire);
        if (next_map_ptr == nullptr) [[unlikely]]
        {
            if (building_next.load(std::memory_order_acquire))
            {
                // Another thread is already building the next map, wait for it to finish
                cpu_relax();
                while (building_next.load(std::memory_order_acquire)) {
                    cpu_relax();
                    cpu_relax();
                    cpu_relax();
                    cpu_relax();
                }
                return next_map.load(std::memory_order_acquire);
            }
            else
            {
                bool expected = false;
                if (building_next.compare_exchange_strong(expected,
                    true,
                    std::memory_order_release,
                    std::memory_order_acquire))
                {
                    next_map_ptr = next_map.load(std::memory_order_acquire);
                    if (next_map_ptr != nullptr)
                    {
                        building_next.store(false, std::memory_order_release);
                        return next_map_ptr;
                    }
                    else
                    {
                        // Build the next map
                        ConcurrentOpenAddressHashMap<N>* new_map = reinterpret_cast<ConcurrentOpenAddressHashMap<N>*>(get_map_metadata_mem());
                        new(new_map) ConcurrentOpenAddressHashMap<N>(capacity);
                        next_map.store(new_map, std::memory_order_release);
                        building_next.store(false, std::memory_order_release);
                        return new_map;
                    }
                }
                else
                {
                    cpu_relax();
                    while (building_next.load(std::memory_order_acquire)) {
                        cpu_relax();
                        cpu_relax();
                        cpu_relax();
                        cpu_relax();
                    }
                    return next_map.load(std::memory_order_acquire);
                }
            }
        }
        else
        {
            return next_map_ptr;
        }
    }
};

#endif
