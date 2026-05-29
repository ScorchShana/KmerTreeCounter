#ifndef CLASSIFIER_THREAD_POOL_HEADER
#define CLASSIFIER_THREAD_POOL_HEADER

#include "definition.h"
#include "MPMCRingQueue.h"
#include "NewKmerTree.h"
#include "ConcurrentMemoryPool.h"
#include "RingMemoryPool.h"
#include "FastqClassifier.h"
#include "SchedulerThreadPool.h"

#include <array>
#include <atomic>
#include <cstdint>
#include <memory>
#include <thread>
#include <vector>
#include <utility>
#include <chrono>
#include <bitset>


template <uint32_t N>
class ClassifierThreadPool
{

    int k_len;
    const uint32_t classifier_count;
    RingMemoryPool<PARSER_CLASSIFIER_RING_MEMORY_POOL_CAPACITY> *parser_classifier_ring_pool;
    SchedulerThreadPool<N> *task_thread_pool;
    KmerTree<N> *tree;
    std::vector<std::unique_ptr<std::thread>> threads_ptr;

public:
#ifdef TEST_MODE
    std::atomic<uint64_t> producer_enqueue_spin_time{0};
    std::atomic<uint64_t> producer_dequeue_spin_time{0};
    std::atomic<uint64_t> consumer_enqueue_spin_time{0};
    std::atomic<uint64_t> consumer_dequeue_spin_time{0};

    std::atomic<uint64_t> classifier_processed_kmers{0};
    std::atomic<uint64_t> prefix_count_cycles{0};
    std::atomic<uint64_t> prefix_reorder_cycles{0};
    std::atomic<uint64_t> classify_local_block_cycles{0};
    std::atomic<uint64_t> tree_add_cycles{0};

    std::atomic<uint64_t> bloom_calls{0};
    std::atomic<uint64_t> bloom_sampled_calls{0};
    std::atomic<uint64_t> bloom_true_count{0};
    std::atomic<uint64_t> bloom_second_layer_attempts{0};
    std::atomic<uint64_t> bloom_sampled_second_layer_attempts{0};
    std::atomic<uint64_t> bloom_first_fetch_or_count{0};
    std::atomic<uint64_t> bloom_second_fetch_or_count{0};
    std::atomic<uint64_t> bloom_sampled_first_fetch_or_count{0};
    std::atomic<uint64_t> bloom_sampled_second_fetch_or_count{0};

    std::atomic<uint64_t> bloom_total_cycles{0};
    std::atomic<uint64_t> bloom_hash_cycles{0};
    std::atomic<uint64_t> bloom_index_mask_cycles{0};
    std::atomic<uint64_t> bloom_first_load_cycles{0};
    std::atomic<uint64_t> bloom_first_filter_cycles{0};
    std::atomic<uint64_t> bloom_first_fetch_or_cycles{0};
    std::atomic<uint64_t> bloom_second_load_cycles{0};
    std::atomic<uint64_t> bloom_second_filter_cycles{0};
    std::atomic<uint64_t> bloom_second_fetch_or_cycles{0};
#endif

    explicit ClassifierThreadPool(const int in_k, KmerTree<N> *tree_ptr,
                                  RingMemoryPool<PARSER_CLASSIFIER_RING_MEMORY_POOL_CAPACITY> *in_parser_classifier_ring_pool_ptr,
                                  SchedulerThreadPool<N> *task_thread_pool_ptr,
                                  uint32_t in_classifier_count)
        : k_len(in_k),
          tree(tree_ptr),
          parser_classifier_ring_pool(in_parser_classifier_ring_pool_ptr),
          classifier_count(in_classifier_count),
          task_thread_pool(task_thread_pool_ptr)
    {
        threads_ptr.reserve(classifier_count);
    }

    void start()
    {
        for (uint32_t i = 0; i < classifier_count; ++i)
        {
            threads_ptr.push_back(std::make_unique<std::thread>([&]
                                                                {
                                                                    std::this_thread::sleep_for(std::chrono::nanoseconds(40));
                                                                    FastqClassifier<N> parser(k_len, parser_classifier_ring_pool, tree);
                                                                    parser.classify_and_push();
                                                                    task_thread_pool->mark_producer_done();
#ifdef TEST_MODE
                                                                    consumer_enqueue_spin_time.fetch_add(parser.consumer_enqueue_spin_time, std::memory_order_relaxed);
                                                                    consumer_dequeue_spin_time.fetch_add(parser.consumer_dequeue_spin_time, std::memory_order_relaxed);
                                                                    producer_enqueue_spin_time.fetch_add(parser.producer_enqueue_spin_time, std::memory_order_relaxed);
                                                                    producer_dequeue_spin_time.fetch_add(parser.producer_dequeue_spin_time, std::memory_order_relaxed);

                                                                    classifier_processed_kmers.fetch_add(parser.classifier_processed_kmers, std::memory_order_relaxed);
                                                                    prefix_count_cycles.fetch_add(parser.prefix_count_cycles, std::memory_order_relaxed);
                                                                    prefix_reorder_cycles.fetch_add(parser.prefix_reorder_cycles, std::memory_order_relaxed);
                                                                    classify_local_block_cycles.fetch_add(parser.classify_local_block_cycles, std::memory_order_relaxed);
                                                                    tree_add_cycles.fetch_add(parser.tree_add_cycles, std::memory_order_relaxed);

                                                                    bloom_calls.fetch_add(parser.bloom_timing_stats.calls, std::memory_order_relaxed);
                                                                    bloom_sampled_calls.fetch_add(parser.bloom_timing_stats.sampled_calls, std::memory_order_relaxed);
                                                                    bloom_true_count.fetch_add(parser.bloom_timing_stats.true_count, std::memory_order_relaxed);
                                                                    bloom_second_layer_attempts.fetch_add(parser.bloom_timing_stats.second_layer_attempts, std::memory_order_relaxed);
                                                                    bloom_sampled_second_layer_attempts.fetch_add(parser.bloom_timing_stats.sampled_second_layer_attempts, std::memory_order_relaxed);
                                                                    bloom_first_fetch_or_count.fetch_add(parser.bloom_timing_stats.first_fetch_or_count, std::memory_order_relaxed);
                                                                    bloom_second_fetch_or_count.fetch_add(parser.bloom_timing_stats.second_fetch_or_count, std::memory_order_relaxed);
                                                                    bloom_sampled_first_fetch_or_count.fetch_add(parser.bloom_timing_stats.sampled_first_fetch_or_count, std::memory_order_relaxed);
                                                                    bloom_sampled_second_fetch_or_count.fetch_add(parser.bloom_timing_stats.sampled_second_fetch_or_count, std::memory_order_relaxed);

                                                                    bloom_total_cycles.fetch_add(parser.bloom_timing_stats.total_cycles, std::memory_order_relaxed);
                                                                    bloom_hash_cycles.fetch_add(parser.bloom_timing_stats.hash_cycles, std::memory_order_relaxed);
                                                                    bloom_index_mask_cycles.fetch_add(parser.bloom_timing_stats.index_mask_cycles, std::memory_order_relaxed);
                                                                    bloom_first_load_cycles.fetch_add(parser.bloom_timing_stats.first_load_cycles, std::memory_order_relaxed);
                                                                    bloom_first_filter_cycles.fetch_add(parser.bloom_timing_stats.first_filter_cycles, std::memory_order_relaxed);
                                                                    bloom_first_fetch_or_cycles.fetch_add(parser.bloom_timing_stats.first_fetch_or_cycles, std::memory_order_relaxed);
                                                                    bloom_second_load_cycles.fetch_add(parser.bloom_timing_stats.second_load_cycles, std::memory_order_relaxed);
                                                                    bloom_second_filter_cycles.fetch_add(parser.bloom_timing_stats.second_filter_cycles, std::memory_order_relaxed);
                                                                    bloom_second_fetch_or_cycles.fetch_add(parser.bloom_timing_stats.second_fetch_or_cycles, std::memory_order_relaxed);
#endif
                                                                }));
        }
    }

    void join()
    {
        for (auto &t : threads_ptr)
        {
            if (t->joinable())
            {
                t->join();
            }
        }

    }
};

#endif
