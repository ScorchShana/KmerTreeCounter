#ifndef FLAT_CONCURRENT_HASH_MAP_HEADER
#define FLAT_CONCURRENT_HASH_MAP_HEADER

#include "../src/definition.h"
#include "../src/kmer.h"
#include "../include/xxh3.h"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <iostream>
#include <sys/mman.h>


#if defined(__AVX2__)
#define FLAT_HASHMAP_AVX2 1
#define FLAT_HASHMAP_GROUP_SIZE 32
#elif defined(__SSE4_2__)
#define FLAT_HASHMAP_SSE 1
#define FLAT_HASHMAP_GROUP_SIZE 16
#else
#define FLAT_HASHMAP_SCALAR 1
#define FLAT_HASHMAP_GROUP_SIZE 8
#endif

#if defined(FLAT_HASHMAP_AVX2) || defined(FLAT_HASHMAP_SSE)
#include <immintrin.h>
#endif

template <uint32_t N>
class FlatConcurrentHashMap
{
public:
    static constexpr double LOAD_FACTOR = 0.75;
    static constexpr uint64_t GROUP_SIZE = FLAT_HASHMAP_GROUP_SIZE;

    using Entry = ExportRecord<N>;

    [[nodiscard]] static uint64_t required_mmap_bytes(uint64_t expected_unique_key)
    {
        return calculate_layout(expected_unique_key).mmap_bytes;
    }

    explicit FlatConcurrentHashMap(uint64_t expected_unique)
        : FlatConcurrentHashMap(calculate_layout(expected_unique))
    {
    }

    ~FlatConcurrentHashMap() noexcept
    {
        if (mmap_base_ != nullptr)
        {
            munmap(mmap_base_, static_cast<size_t>(mmap_size_));
            mmap_base_ = nullptr;
            mmap_size_ = 0;
        }
    }

    FlatConcurrentHashMap(const FlatConcurrentHashMap&) = delete;
    FlatConcurrentHashMap& operator=(const FlatConcurrentHashMap&) = delete;
    FlatConcurrentHashMap(FlatConcurrentHashMap&&) = delete;
    FlatConcurrentHashMap& operator=(FlatConcurrentHashMap&&) = delete;

    bool insert_unique(const kmer<N>& key, uint32_t count) noexcept
    {
        if (sealed_.load(std::memory_order_acquire)) [[unlikely]]
        {
            return false;
        }

        const uint64_t hash = XXH3_64bits(&key, sizeof(key));
        const uint8_t fp = fingerprint(hash);
        uint64_t base = initial_slot(hash);

        for (uint64_t probe = 0; probe < group_count_; ++probe)
        {
            for (uint64_t offset = 0; offset < GROUP_SIZE; ++offset)
            {
                const uint64_t slot = wrap_slot(base + offset);
                std::atomic_ref<uint8_t> ctrl_ref(ctrl_[static_cast<size_t>(slot)]);
                uint8_t expected = 0;
                if (ctrl_ref.compare_exchange_strong(
                    expected,
                    fp,
                    std::memory_order_acq_rel,
                    std::memory_order_acquire))
                {
                    std::construct_at(entries_ + static_cast<size_t>(slot), Entry{ key, count });
                    size_.fetch_add(1, std::memory_order_relaxed);
                    return true;
                }
            }
            base = next_group_base(base);
        }

        return false;
    }

    void seal() noexcept
    {
        for (uint64_t i = 0; i < GROUP_SIZE - 1; ++i)
        {
            ctrl_[static_cast<size_t>(capacity_ + i)] = ctrl_[static_cast<size_t>(i)];
        }
        sealed_.store(true, std::memory_order_release);
    }

    bool find(const kmer<N>& key, uint32_t& out_count) const noexcept
    {
        if (!sealed_.load(std::memory_order_acquire)) [[unlikely]]
        {
            return false;
        }

        const uint64_t hash = XXH3_64bits(&key, sizeof(key));
        const uint8_t fp = fingerprint(hash);
        uint64_t base = initial_slot(hash);

        for (uint64_t probe = 0; probe < group_count_; ++probe)
        {
            const auto masks = match_and_empty(base, fp);
            uint32_t match_mask = masks.first;
            const uint32_t empty_mask = masks.second;

            while (match_mask != 0)
            {
                const uint32_t bit = static_cast<uint32_t>(__builtin_ctz(match_mask));
                const uint64_t slot = wrap_slot(base + bit);
                const Entry& entry = entries_[static_cast<size_t>(slot)];
                if (entry.key == key) [[likely]]
                {
                    out_count = entry.count;
                    return true;
                }
                match_mask &= match_mask - 1;
            }

            if (empty_mask != 0)
            {
                return false;
            }
            base = next_group_base(base);
        }

        return false;
    }

    bool contains(const kmer<N>& key) const noexcept
    {
        uint32_t ignored = 0;
        return find(key, ignored);
    }

    uint64_t size() const noexcept
    {
        return size_.load(std::memory_order_relaxed);
    }

    uint64_t capacity() const noexcept
    {
        return capacity_;
    }

private:
    static constexpr uint64_t HUGE_PAGE_BYTES = 2ULL * 1024ULL * 1024ULL;
    static_assert(std::is_trivially_destructible_v<Entry>, "Entry must be trivially destructible");

    struct Layout
    {
        uint64_t capacity;
        uint64_t ctrl_bytes;
        uint64_t entries_offset;
        uint64_t entries_bytes;
        uint64_t payload_bytes;
        uint64_t mmap_bytes;
    };

    explicit FlatConcurrentHashMap(const Layout& layout)
        : capacity_(layout.capacity),
        group_count_(capacity_ / GROUP_SIZE),
        mmap_base_(map_memory(layout.mmap_bytes)),
        mmap_size_(layout.mmap_bytes),
        ctrl_(static_cast<uint8_t*>(mmap_base_)),
        entries_(reinterpret_cast<Entry*>(static_cast<uint8_t*>(mmap_base_) + layout.entries_offset))
    {
        std::memset(ctrl_, 0, static_cast<size_t>(layout.ctrl_bytes));
    }

    static uint64_t round_up(uint64_t value, uint64_t alignment) noexcept
    {
        const uint64_t remainder = value % alignment;
        return remainder == 0 ? value : value + alignment - remainder;
    }

    static uint64_t calculate_capacity(uint64_t expected_unique) noexcept
    {
        if (expected_unique == 0)
        {
            return GROUP_SIZE;
        }

        const long double raw =
            static_cast<long double>(expected_unique) / static_cast<long double>(LOAD_FACTOR);
        uint64_t requested = static_cast<uint64_t>(raw);
        if (static_cast<long double>(requested) < raw)
        {
            ++requested;
        }
        if (requested < GROUP_SIZE)
        {
            requested = GROUP_SIZE;
        }
        if (requested > std::numeric_limits<uint64_t>::max() - GROUP_SIZE) [[unlikely]]
        {
            return std::numeric_limits<uint64_t>::max() / GROUP_SIZE * GROUP_SIZE;
        }
        return round_up(requested, GROUP_SIZE);
    }

    static Layout calculate_layout(uint64_t expected_unique_key)
    {
        const uint64_t capacity = calculate_capacity(expected_unique_key);

        const uint64_t ctrl_bytes = capacity + GROUP_SIZE - 1;
        const uint64_t entries_offset = (CACHE_LINE_SIZE % alignof(Entry) == 0) ? round_up(ctrl_bytes, CACHE_LINE_SIZE) : round_up(ctrl_bytes, alignof(Entry));

        const uint64_t entries_bytes = capacity * sizeof(Entry);

        const uint64_t payload_bytes = entries_offset + entries_bytes;
        const uint64_t mmap_bytes = round_up(payload_bytes, HUGE_PAGE_BYTES);

        return Layout{ capacity, ctrl_bytes, entries_offset, entries_bytes, payload_bytes, mmap_bytes };
    }

    static void* map_memory(uint64_t mmap_bytes)
    {
        if (mmap_bytes == 0 || mmap_bytes > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) [[unlikely]]
        {
            throw std::bad_alloc();
        }

        const size_t length = static_cast<size_t>(mmap_bytes);
        void* base = mmap(nullptr,
            length,
            PROT_READ | PROT_WRITE,
            MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB,
            -1,
            0);

        if (base != MAP_FAILED)
        {
            return base;
        }

        base = mmap(nullptr,
            length,
            PROT_READ | PROT_WRITE,
            MAP_PRIVATE | MAP_ANONYMOUS,
            -1,
            0);

        if (base == MAP_FAILED)
        {
            throw std::bad_alloc();
        }

        madvise(base, length, MADV_HUGEPAGE);
        return base;
    }

    static uint8_t fingerprint(uint64_t hash) noexcept
    {
        return static_cast<uint8_t>((hash & 0x7FULL) | 0x80U);
    }

    uint64_t initial_slot(uint64_t hash) const noexcept
    {
        return (hash >> 7) % capacity_;
    }

    uint64_t wrap_slot(uint64_t slot) const noexcept
    {
        return slot >= capacity_ ? slot - capacity_ : slot;
    }

    uint64_t next_group_base(uint64_t base) const noexcept
    {
        base += GROUP_SIZE;
        return base >= capacity_ ? base - capacity_ : base;
    }

    std::pair<uint32_t, uint32_t> match_and_empty(uint64_t base, uint8_t fp) const noexcept
    {
        const uint8_t* ctrl = ctrl_ + base;

#if defined(FLAT_HASHMAP_AVX2)
        const __m256i ctrl_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(ctrl));
        const __m256i fp_vec = _mm256_set1_epi8(static_cast<char>(fp));
        const __m256i zero_vec = _mm256_setzero_si256();
        const uint32_t match_mask =
            static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(ctrl_vec, fp_vec)));
        const uint32_t empty_mask =
            static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(ctrl_vec, zero_vec)));
        return { match_mask, empty_mask };
#elif defined(FLAT_HASHMAP_SSE)
        const __m128i ctrl_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(ctrl));
        const __m128i fp_vec = _mm_set1_epi8(static_cast<char>(fp));
        const __m128i zero_vec = _mm_setzero_si128();
        const uint32_t match_mask =
            static_cast<uint32_t>(_mm_movemask_epi8(_mm_cmpeq_epi8(ctrl_vec, fp_vec)));
        const uint32_t empty_mask =
            static_cast<uint32_t>(_mm_movemask_epi8(_mm_cmpeq_epi8(ctrl_vec, zero_vec)));
        return { match_mask, empty_mask };
#else
        uint32_t match_mask = 0;
        uint32_t empty_mask = 0;
        for (uint64_t i = 0; i < GROUP_SIZE; ++i)
        {
            const uint8_t c = ctrl[i];
            if (c == fp)
            {
                match_mask |= (1U << i);
            }
            if (c == 0)
            {
                empty_mask |= (1U << i);
            }
    }
        return { match_mask, empty_mask };
#endif
}

    const uint64_t capacity_;
    const uint64_t group_count_;
    void* mmap_base_ = nullptr;
    uint64_t mmap_size_ = 0;
    uint8_t* ctrl_ = nullptr;
    Entry* entries_ = nullptr;
    alignas(64) std::atomic<uint64_t> size_{ 0 };
    std::atomic<bool> sealed_{ false };
};

#endif
