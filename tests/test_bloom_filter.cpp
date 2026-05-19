#include "../src/BloomFilter.h"

#include <algorithm>
#include <atomic>
#include <barrier>
#include <chrono>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <thread>
#include <unordered_map>
#include <vector>

namespace
{

    constexpr size_t kGlobalHotKeyCount = 1;
    constexpr size_t kLocalHotKeyCount = 1;

    struct Config
    {
        uint64_t total_inserts = 0;
        uint64_t inserts_per_thread = 0;
        uint64_t num_threads = 0;
        double duplicate_prob = 0.0;
    };

    void print_usage(const char *argv0)
    {
        std::cerr << "Usage: " << argv0
                  << " <inserts_per_thread> <num_threads> [duplicate_prob]\n";
    }

    bool parse_u64(const char *text, uint64_t &value)
    {
        if (!text || *text == '\0')
        {
            return false;
        }
        char *end = nullptr;
        errno = 0;
        unsigned long long parsed = std::strtoull(text, &end, 10);
        if (errno != 0 || end == text || *end != '\0')
        {
            return false;
        }
        value = static_cast<uint64_t>(parsed);
        return true;
    }

    bool parse_double(const char *text, double &value)
    {
        if (!text || *text == '\0')
        {
            return false;
        }
        char *end = nullptr;
        errno = 0;
        double parsed = std::strtod(text, &end);
        if (errno != 0 || end == text || *end != '\0')
        {
            return false;
        }
        value = parsed;
        return true;
    }

    uint64_t next_power_of_two(uint64_t value)
    {
        if (value == 0)
        {
            return 1;
        }
        --value;
        value |= value >> 1;
        value |= value >> 2;
        value |= value >> 4;
        value |= value >> 8;
        value |= value >> 16;
        value |= value >> 32;
        return value + 1;
    }

    template <uint32_t N>
    struct KmerHash
    {
        size_t operator()(const kmer<N> &value) const noexcept
        {
            return static_cast<size_t>(hash_func(value));
        }
    };

    template <uint32_t N>
    kmer<N> make_random_kmer(std::mt19937_64 &rng, std::uniform_int_distribution<uint32_t> &base_dist)
    {
        kmer<N> key{};
        constexpr uint32_t kBasesPerWord = 32;
        for (uint32_t word = 0; word < N; ++word)
        {
            uint64_t packed = 0;
            for (uint32_t i = 0; i < kBasesPerWord; ++i)
            {
                packed = (packed << 2) | static_cast<uint64_t>(base_dist(rng));
            }
            key.data[word] = packed;
        }
        return key;
    }

    std::vector<uint8_t> make_base_sequence(size_t length, std::mt19937_64 &rng)
    {
        std::vector<uint8_t> bases(length);
        std::uniform_int_distribution<uint32_t> base_dist(0, 3);
        for (auto &base : bases)
        {
            base = static_cast<uint8_t>(base_dist(rng));
        }
        return bases;
    }

    template <uint32_t N>
    kmer<N> make_kmer_from_bases(const uint8_t *bases)
    {
        kmer<N> key{};
        constexpr uint32_t kBasesPerWord = 32;
        for (uint32_t word = 0; word < N; ++word)
        {
            uint64_t packed = 0;
            const uint8_t *word_bases = bases + (word * kBasesPerWord);
            for (uint32_t i = 0; i < kBasesPerWord; ++i)
            {
                packed = (packed << 2) | static_cast<uint64_t>(word_bases[i] & 0x3);
            }
            key.data[word] = packed;
        }
        return key;
    }

    template <uint32_t N>
    struct ExpectedStats
    {
        uint64_t unique_count = 0;
        uint64_t expected_true = 0;
    };

    template <uint32_t N>
    ExpectedStats<N> compute_expected_stats(const std::vector<std::vector<kmer<N>>> &thread_keys, uint64_t total_inserts)
    {
        std::unordered_map<kmer<N>, uint32_t, KmerHash<N>> counts;
        counts.reserve(static_cast<size_t>(total_inserts));

        for (const auto &keys : thread_keys)
        {
            for (const auto &key : keys)
            {
                ++counts[key];
            }
        }

        ExpectedStats<N> stats{};
        stats.unique_count = static_cast<uint64_t>(counts.size());
        for (const auto &entry : counts)
        {
            stats.expected_true += std::min<uint32_t>(entry.second, 2u);
        }
        return stats;
    }

    template <uint32_t N>
    void generate_thread_keys(const Config &cfg,
                              const std::vector<kmer<N>> &global_hot_keys,
                              std::vector<std::vector<kmer<N>>> &thread_keys)
    {
        constexpr uint32_t kBasesPerWord = 32;
        constexpr uint32_t kKmerBases = kBasesPerWord * N;

        thread_keys.resize(cfg.num_threads);

        for (uint64_t thread_index = 0; thread_index < cfg.num_threads; ++thread_index)
        {
            std::mt19937_64 rng(0xC0FFEEULL + thread_index);
            std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
            std::uniform_int_distribution<uint32_t> base_dist(0, 3);

            std::vector<kmer<N>> local_hot_keys;
            local_hot_keys.reserve(kLocalHotKeyCount);
            for (size_t i = 0; i < kLocalHotKeyCount; ++i)
            {
                local_hot_keys.push_back(make_random_kmer<N>(rng, base_dist));
            }

            std::vector<kmer<N>> sequence_kmers;
            if (cfg.duplicate_prob < 1.0)
            {
                // Build overlapping k-mers from a base stream to mimic real k-mer locality.
                const size_t base_count = cfg.inserts_per_thread + kKmerBases - 1;
                std::vector<uint8_t> bases = make_base_sequence(base_count, rng);
                sequence_kmers.reserve(cfg.inserts_per_thread);
                for (uint64_t i = 0; i < cfg.inserts_per_thread; ++i)
                {
                    sequence_kmers.push_back(make_kmer_from_bases<N>(bases.data() + i));
                }
            }

            std::vector<kmer<N>> keys;
            keys.reserve(cfg.inserts_per_thread);

            std::uniform_int_distribution<size_t> global_index(0, global_hot_keys.size() - 1);
            std::uniform_int_distribution<size_t> local_index(0, local_hot_keys.size() - 1);

            for (uint64_t i = 0; i < cfg.inserts_per_thread; ++i)
            {
                const double roll = prob_dist(rng);
                if (roll < cfg.duplicate_prob)
                {
                    const bool use_global = prob_dist(rng) < 0.5;
                    if (use_global)
                    {
                        keys.push_back(global_hot_keys[global_index(rng)]);
                    }
                    else
                    {
                        keys.push_back(local_hot_keys[local_index(rng)]);
                    }
                }
                else
                {
                    keys.push_back(sequence_kmers[i]);
                }
            }

            thread_keys[thread_index] = std::move(keys);
        }
    }

    template <uint32_t N>
    int run_bloom_filter_test(const Config &cfg)
    {
        using Kmer = kmer<N>;

        std::mt19937_64 global_rng(0xABCDEFULL + N);
        std::uniform_int_distribution<uint32_t> base_dist(0, 3);
        std::vector<Kmer> global_hot_keys;
        global_hot_keys.reserve(kGlobalHotKeyCount);
        for (size_t i = 0; i < kGlobalHotKeyCount; ++i)
        {
            global_hot_keys.push_back(make_random_kmer<N>(global_rng, base_dist));
        }

        std::vector<std::vector<Kmer>> thread_keys;
        generate_thread_keys<N>(cfg, global_hot_keys, thread_keys);

        const ExpectedStats<N> expected = compute_expected_stats<N>(thread_keys, cfg.total_inserts);
        const uint64_t unique_count = expected.unique_count;
        const uint64_t expected_true = expected.expected_true;

        // 计算 0.1% (0.001) 误判率所需的 capacity:
        // 因为是两层独立的 Bloom Filter，每层 P_layer = sqrt(0.001) ≈ 0.03162
        // 且每层 k = 3，根据公式 P = (1 - e^(-k * n / m))^k 
        // 求解得到 m/n ≈ -3 / ln(1 - 0.03162^(1/3)) ≈ 7.893 bits/element
        // 每层容量为 capacity 个 uint64_t (即 capacity * 64 bits)
        // 所以 capacity * 64 = 7.893 * n => capacity = n * 7.893 / 64
        uint64_t target_capacity = static_cast<uint64_t>(unique_count * 7.893 / 64.0);
        const uint64_t capacity = next_power_of_two(target_capacity == 0 ? 1 : target_capacity);
        ConcurrentDoubleBloomFilter<N> filter(capacity);

        std::atomic<long long> max_insert_ns{0};
        std::barrier sync_point(static_cast<std::ptrdiff_t>(cfg.num_threads));
        std::vector<std::thread> threads;
        std::vector<uint64_t> true_counts(cfg.num_threads, 0);
        threads.reserve(cfg.num_threads);

        for (uint64_t thread_index = 0; thread_index < cfg.num_threads; ++thread_index)
        {
            threads.emplace_back([&, thread_index]()
                                 {
                                     const auto &keys = thread_keys[thread_index];
                                     sync_point.arrive_and_wait();

                                     const auto start_tp = std::chrono::steady_clock::now();
                                     uint64_t local_true = 0;
                                     for (const auto &key : keys)
                                     {
                                         if (filter.insert(key))
                                         {
                                             ++local_true;
                                         }
                                     }
                                     const auto end_tp = std::chrono::steady_clock::now();
                                     true_counts[thread_index] = local_true;

                                     const auto elapsed_ns =
                                         std::chrono::duration_cast<std::chrono::nanoseconds>(end_tp - start_tp).count();
                                     long long observed = max_insert_ns.load(std::memory_order_relaxed);
                                     while (elapsed_ns > observed &&
                                            !max_insert_ns.compare_exchange_weak(observed, elapsed_ns,
                                                                                std::memory_order_release,
                                                                                std::memory_order_relaxed))
                                     {
                                     }
                                 });
        }

        for (auto &thread : threads)
        {
            thread.join();
        }

        const long long max_ns = max_insert_ns.load(std::memory_order_acquire);
        if (max_ns <= 0)
        {
            std::cerr << "timing was not recorded\n";
            return 1;
        }

        uint64_t observed_true = 0;
        for (const auto count : true_counts)
        {
            observed_true += count;
        }

        const uint64_t false_positive = observed_true < expected_true ? (expected_true - observed_true) : 0;
        const uint64_t overcount = observed_true > expected_true ? (observed_true - expected_true) : 0;
        const double false_positive_rate = expected_true > 0
                                              ? static_cast<double>(false_positive) / static_cast<double>(expected_true)
                                              : 0.0;

        const double elapsed_seconds = static_cast<double>(max_ns) * 1e-9;
        const double ops_per_sec = elapsed_seconds > 0.0
                                       ? static_cast<double>(cfg.total_inserts) / elapsed_seconds
                                       : 0.0;

        std::cout << "kmer<" << N << ">\n"
                  << "capacity=" << capacity << "\n"
                  << "total_inserts=" << cfg.total_inserts << "\n"
                  << "unique_kmers=" << unique_count << "\n"
                  << "expected_true=" << expected_true << "\n"
                  << "observed_true=" << observed_true << "\n"
                  << "false_positive=" << false_positive << "\n"
                  << "overcount=" << overcount << "\n"
                  << std::fixed << std::setprecision(6)
                  << "false_positive_rate=" << false_positive_rate << "\n"
                  << "insert_seconds=" << elapsed_seconds << "\n"
                  << (ops_per_sec / 1000000.0) << "Mopt/s\n";

        return 0;
    }

} // namespace

int main(int argc, char **argv)
{
    Config cfg{};
    if ((argc != 3 && argc != 4) ||
        !parse_u64(argv[1], cfg.inserts_per_thread) ||
        !parse_u64(argv[2], cfg.num_threads) ||
        (argc == 4 && !parse_double(argv[3], cfg.duplicate_prob)))
    {
        print_usage(argv[0]);
        return 1;
    }

    if (cfg.inserts_per_thread == 0 || cfg.num_threads == 0)
    {
        std::cerr << "inserts_per_thread and num_threads must be > 0\n";
        return 1;
    }

    if (cfg.duplicate_prob < 0.0 || cfg.duplicate_prob > 1.0)
    {
        std::cerr << "duplicate_prob must be in [0.0, 1.0]\n";
        return 1;
    }

    if (cfg.inserts_per_thread > std::numeric_limits<uint64_t>::max() / cfg.num_threads)
    {
        std::cerr << "inserts_per_thread * num_threads overflow\n";
        return 1;
    }

    cfg.total_inserts = cfg.inserts_per_thread * cfg.num_threads;

    int rc = 0;
    rc |= run_bloom_filter_test<1>(cfg);
    rc |= run_bloom_filter_test<2>(cfg);
    rc |= run_bloom_filter_test<4>(cfg);

    return rc;
}
