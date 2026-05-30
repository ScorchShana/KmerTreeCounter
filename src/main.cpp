#include "NewKmerTree.h"
#include "SchedulerThreadPool.h"
#include "ConcurrentMemoryPool.h"
#include "MPMCRingQueue.h"
#include "GetKmer.h"
#include "FastqReader.h"
#include "ParserThreadPool.h"
#include "ExportWriter.h"
#include "ClassifierThreadPool.h"
#include "FastqPrefixCounter.h"
#include "FastqPreReader.h"
#include "FastqParser.h"

#include <chrono>
#include <iostream>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <memory>
#include <limits>
#include <format>
#include <vector>
#include <thread>
#include <algorithm>
#include <bit>

uint32_t k_len;        // k-mer 的长度 (例如: 31, 41)
uint32_t n_thread;     // 允许使用的总线程数
uint64_t memory_limit; // 全局最大内存限制，单位为 GB

int main(int argc, char *argv[])
{

    if (argc < 5 || argc > 9)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <fastq_file> <k_len> <n_thread> <memory_limit_gb> [map_capacity] [min_count] [max_count] [parser_threads]" << std::endl;
        return 1;
    }
    bool enable_prefault = false;
    uint32_t parser_threads_override = 0;
    min_count = 1;
    max_count = std::numeric_limits<uint32_t>::max();

    std::string filename;
    try
    {
        filename = argv[1];
        k_len = std::stoul(argv[2]);
        n_thread = std::stoul(argv[3]);
        memory_limit = std::stoull(argv[4]);

        if (argc >= 6)
        {
            kmer_concurrent_hash_map_capacity = std::stoul(argv[5]);
        }
        if (argc >= 7)
        {
            min_count = std::stoul(argv[6]);
        }
        if (argc >= 8)
        {
            max_count = std::stoul(argv[7]);
        }

        if (n_thread < 6)
        {
            std::cerr << "n_thread must be >= 6" << std::endl;
            std::exit(-1);
        }

        std::cout << "Input parameters:" << std::endl;
        std::cout << "  Fastq file: " << filename << std::endl;
        std::cout << "  k-mer length: " << k_len << std::endl;
        std::cout << "  Thread count: " << n_thread << std::endl;
        std::cout << "  Memory limit (GB): " << memory_limit << std::endl;
        std::cout << "  Map capacity: " << kmer_concurrent_hash_map_capacity << std::endl;
        std::cout << "  Min count: " << min_count << std::endl;
        std::cout << "  Max count: " << max_count << std::endl;
    }
    catch (const std::exception &)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <fastq_file> <k_len> <n_thread> <memory_limit_gb> [map_capacity] [min_count] [max_count] [parser_threads]" << std::endl;
        return 1;
    }

    if (kmer_concurrent_hash_map_capacity <= 1 || kmer_concurrent_hash_map_capacity >= 16ULL * 1024 * 1024 || max_count < min_count)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <fastq_file> <k_len> <n_thread> <memory_limit_gb> [map_capacity] [min_count] [max_count] [parser_threads]" << std::endl;
        return 1;
    }

    constexpr uint32_t N1 = 2;

    // 根据预算计算分给 Worker 的解析(Parser)线程和任务(Tasker)线程的数量
    // Worker 预算去除了主线程(Reader)和导出线程(ExportWriter)
    const uint32_t parser_num = (n_thread / 8 > 0) ? (n_thread / 8) : 1; // 预留至少 1 个线程给 Parser，剩余线程在 Parser 和 Tasker 之间分配
    const uint32_t worker_budget = n_thread - 2 - parser_num;
    const uint32_t classifier_num = (worker_budget / (1.0 + TASK_CLASSIFIER_RATIO) == 0) ? 1 : (uint32_t)(worker_budget / (1.0 + TASK_CLASSIFIER_RATIO));
    const uint32_t tasker_num = worker_budget - classifier_num;

    std::cout << "Thread split:" << std::endl;
    std::cout << "  parser threads: " << parser_num << std::endl;
    std::cout << "  classifier threads: " << classifier_num << std::endl;
    std::cout << "  task threads: " << tasker_num << std::endl;

    const auto init_start = std::chrono::steady_clock::now();

    // 初始化层级队列，用于在树的不同深度间传递任务
    auto layer_queues = std::make_shared<LayerQueues<N1>>();
    // 初始化解析器环形内存池，管理 Reader 读取后的碱基字符串数据块
    auto reader_parser_ring_pool = std::make_shared<SPMCRingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY>>(READER_PARSER_RING_MEMORY_POOL_BLOCK_SIZE, 1);
    // 初始化分类器环形内存池，管理 Parser 线程处理后的 k-mer 数据块
    auto parser_classifier_ring_pool = std::make_shared<RingMemoryPool<PARSER_CLASSIFIER_RING_MEMORY_POOL_CAPACITY>>(PARSER_CLASSIFIER_RING_MEMORY_POOL_BLOCK_SIZE, 1);
    // 初始化导出用的环形内存池，管理低频 k-mer 的导出数据块
    auto export_ring_pool = std::make_shared<RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY>>(EXPORT_RING_MEMORY_POOL_BLOCK_SIZE, 1);
    // 初始化全局并发内存池，用于节点分配、哈希表等
    auto pool = std::make_shared<ConcurrentMemoryPool>(memory_limit * 1024ULL * 1024ULL * 1024ULL);

    // PreRead阶段
    FastqPreReader<N1> pre_reader(filename, k_len, FASTQ_FILE_CHUNK_SIZE, reader_parser_ring_pool.get());
    auto pre_parser_thread = std::thread([&]
                                         {
        FastqParser<N1> parser(k_len, reader_parser_ring_pool.get(), parser_classifier_ring_pool.get());
        parser.parse_and_push();
        parser_classifier_ring_pool->producer_set_finished(); });

    std::vector<std::atomic<uint32_t>> prefix_counts(1U << (2 * ROOT_BASES)); // 256 个前缀的计数器
    for (auto &v : prefix_counts)
    {
        v.store(0); // 或 v.store(init_value);
    }
    std::vector<std::thread> prefix_counter_threads;
    uint32_t prefix_counter_thread_num = std::min(std::max(1u, n_thread - 2), 8u); // 预留至少 1 个线程给前缀计数器
    for (uint32_t i = 0; i < prefix_counter_thread_num; ++i)
    {
        prefix_counter_threads.emplace_back([&]
                                            {
            FastqPrefixCounter<N1> prefix_counter(k_len, parser_classifier_ring_pool.get());
            prefix_counter.count_prefixes();
            for(size_t j = 0; j < prefix_counts.size(); j++)
            {
                prefix_counts[j].fetch_add(prefix_counter.prefix_counts[j], std::memory_order_relaxed);
            } });
    }
    pre_reader.pre_read();
    pre_parser_thread.join();
    for (auto &t : prefix_counter_threads)
    {
        t.join();
    }

    uint64_t average_count = 0;

    for (uint64_t i = 0; i < prefix_counts.size(); i++)
    {
        if (prefix_counts[i].load(std::memory_order_relaxed) < 1000)
        {
            prefix_counts[i].store(1000, std::memory_order_relaxed);
        }
        average_count += prefix_counts[i].load(std::memory_order_relaxed);
        // std::cout << std::format("Prefix {:0>8b}: count = {}", i, prefix_counts[i].load()) << std::endl;
    }
    average_count /= prefix_counts.size();
    std::cout << "Average prefix count: " << average_count << std::endl;

    uint8_t bloom_index = 0;
    uint64_t partition_prefix_sum = 0;
    for (uint64_t i = 0; i < prefix_counts.size(); i++)
    {
        if (prefix_counts[i].load(std::memory_order_relaxed) > average_count * 4)
        {
            
            bloom_filter_capacity[bloom_index] = std::bit_ceil(std::max<uint64_t>(1, partition_prefix_sum / (average_count * 4)) * standard_bloom_filter_capacity);
            bloom_index++;
            prefix_hot[i] = 1;
            bloom_filter_capacity[bloom_index] = std::bit_ceil(prefix_counts[i].load(std::memory_order_relaxed) / (average_count * 4) * standard_bloom_filter_capacity);
            prefix_to_bloom_index[i] = bloom_index++;
            partition_prefix_sum = 0;
            std::cout << std::format("Prefix {:0>8b} is hot with count = {}", i, prefix_counts[i].load(std::memory_order_relaxed)) << std::endl;
            std::cout << std::bit_ceil(prefix_counts[i].load(std::memory_order_relaxed) / (average_count * 4)) << " times of standard bloom filter" << std::endl;
        }
        else
        {
            partition_prefix_sum += prefix_counts[i].load(std::memory_order_relaxed);
            if (partition_prefix_sum > average_count * 4)
            {
                bloom_filter_capacity[bloom_index] = std::bit_ceil(partition_prefix_sum / (average_count * 4) * standard_bloom_filter_capacity);
                prefix_to_bloom_index[i] = bloom_index++;
                partition_prefix_sum = 0;
            }
            else
            {
                prefix_to_bloom_index[i] = bloom_index;
            }

            prefix_hot[i] = 0;
        }
    }

    if (partition_prefix_sum > 0)
    {
        const uint64_t threshold = average_count * 4;
        const uint64_t partition_scale = (partition_prefix_sum + threshold - 1) / threshold;
        bloom_filter_capacity[bloom_index] = std::bit_ceil(std::max<uint64_t>(1, partition_scale) * standard_bloom_filter_capacity);
    }

    std::cout << "Bloom filter partition count: " << bloom_index + 1 << std::endl;

    // 初始化 k-mer 字典树(KmerTree)的核心结构，整合前述多个组件
    auto tree = std::make_shared<KmerTree<N1>>(k_len, pool.get(), layer_queues.get(), export_ring_pool.get());

    // 初始化并构建 Tasker 线程池，负责消费层级队列并将 k-mer 路由到深层节点 / 哈希表
    auto task_thread_pool = std::make_shared<SchedulerThreadPool<N1>>(tasker_num, classifier_num, tree.get(), layer_queues.get());
    // 初始化并构建 Parser 线程池，负责消费 FASTQ 读取器产生的数据，提取 k-mer 进行初步布隆过滤
    auto parser_thread_pool = std::make_shared<ParserThreadPool<N1>>(k_len, pool.get(), reader_parser_ring_pool.get(), parser_classifier_ring_pool.get(), parser_num);
    // 初始化并且构建 Classifier 线程池，负责消费 Parser 线程产生的 k-mer 数据块，进行更精细的分类和路由
    auto classifier_thread_pool = std::make_shared<ClassifierThreadPool<N1>>(k_len, tree.get(), parser_classifier_ring_pool.get(), task_thread_pool.get(), classifier_num);

    // 初始化导出写入器，用于将低频 k-mer 单线程安全落盘
    auto export_writer = std::make_shared<ExportWriter<N1>>(export_ring_pool.get());

    // 正式计数阶段
    reader_parser_ring_pool->reset_producers(1); // Reader 线程是单生产者
    parser_classifier_ring_pool->reset_producers(parser_num);
    // 初始化最终排干(drain)阶段所需的初始任务
    layer_queues->initialize_final_drain_queue(tree->root_nodes);

    // 初始化 FASTQ 读取器，将大文件分块读取并送入 ring_pool 用作流水线起点
    FastqReader<N1> reader(filename, k_len, FASTQ_FILE_CHUNK_SIZE, reader_parser_ring_pool.get());
    // FastqClassifier<N1> classifier(k_len, parser_classifier_ring_pool.get(), tree.get());

    const auto init_end = std::chrono::steady_clock::now();

    const auto mid_start = std::chrono::steady_clock::now();

    // 启动多线程流水线的各个工作线程
    parser_thread_pool->start();
    classifier_thread_pool->start();
    task_thread_pool->start();
    export_writer->start();

    const auto read_start = std::chrono::steady_clock::now();

    // 在主线程中运行读取器，读取文件块并推送给 Parser
    reader.read();

    const auto read_end = std::chrono::steady_clock::now();

    // classifier.classify_and_push();
    // task_thread_pool->mark_producer_done(); // 显式标记 Classifier 的生产者角色完成，通知 Scheduler 不再有新数据

    // 阻塞等待 Parser 线程将所有文件块解析并生成 k-mer 完成
    parser_thread_pool->join();
    const auto parser_end = std::chrono::steady_clock::now();
    classifier_thread_pool->join();
    const auto classifier_end = std::chrono::steady_clock::now();
    // 阻塞等待 Tasker 线程消费完所有的分发任务
    task_thread_pool->join();

    const auto task_end = std::chrono::steady_clock::now();

    const auto mid_end = std::chrono::steady_clock::now();

    // 标记所有的缓存已经发送完毕，向导出环形队列发送 nullptr 或退出标志
    const auto export_join_start = std::chrono::steady_clock::now();
    tree->mark_finish_export();

    // 阻塞等待导出器将所有低频 k-mer 都安全写入磁盘完成
    export_writer->join();
    const auto export_join_end = std::chrono::steady_clock::now();

    const auto final_start = std::chrono::steady_clock::now();

    // Final drain 阶段：多线程并行遍历整个字典树，将在节点中暂存但未下发的 k-mers 全部合并到全局哈希表中
    tree->final_drain_parallel(n_thread);

    std::cout << "Total read k-mer count: " << parser_thread_pool->get_total_read_kmer() << std::endl;

    const auto final_end = std::chrono::steady_clock::now();

    const auto init_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(init_end - init_start).count();
    const auto read_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(read_end - read_start).count();
    const auto mid_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(mid_end - mid_start).count();
    const auto export_join_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(export_join_end - export_join_start).count();
    const auto final_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(final_end - final_start).count();
    const auto total_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(final_end - init_start).count();
    const auto classifier_task_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(task_end - classifier_end).count();
    const auto parse_classifier_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(classifier_end - parser_end).count();
    const auto read_parse_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(parser_end - read_end).count();

    std::cout << "Init elapsed us: " << init_elapsed_us << std::endl;
    std::cout << "Read elapsed us: " << read_elapsed_us << std::endl;
    std::cout << "Mid elapsed us: " << mid_elapsed_us << std::endl;
    std::cout << "Between Read and Parse elapsed us: " << read_parse_elapsed_us << std::endl;
    std::cout << "Between Parse and Classifier elapsed us: " << parse_classifier_elapsed_us << std::endl;
    std::cout << "Between Classifier and Task elapsed us: " << classifier_task_elapsed_us << std::endl;
    std::cout << "Export join elapsed us: " << export_join_elapsed_us << std::endl;
    std::cout << "Final elapsed us: " << final_elapsed_us << std::endl;
    std::cout << "Total elapsed us: " << total_elapsed_us << std::endl;

#ifdef TEST_MODE
    SpinLock::flush_spin_loops_for_current_thread();
    std::cout << "SpinLock spin_loops: " << SpinLock::spin_loops() << std::endl;
    std::cout << "Parser producer enqueue total spin time: " << parser_thread_pool->producer_enqueue_spin_time.load() << std::endl;
    std::cout << "Parser producer dequeue total spin time: " << parser_thread_pool->producer_dequeue_spin_time.load() << std::endl;
    std::cout << "Parser consumer dequeue total spin time: " << parser_thread_pool->consumer_dequeue_spin_time.load() << std::endl;
    std::cout << "Parser consumer enqueue total spin time: " << parser_thread_pool->consumer_enqueue_spin_time.load() << std::endl;
    std::cout << "Classifier consumer dequeue total spin time: " << classifier_thread_pool->consumer_dequeue_spin_time.load() << std::endl;
    std::cout << "Classifier consumer enqueue total spin time: " << classifier_thread_pool->consumer_enqueue_spin_time.load() << std::endl;
    std::cout << "Classifier producer enqueue total spin time: " << classifier_thread_pool->producer_enqueue_spin_time.load() << std::endl;
    std::cout << "Classifier producer dequeue total spin time: " << classifier_thread_pool->producer_dequeue_spin_time.load() << std::endl;

    auto per = [](uint64_t value, uint64_t count) -> double
    {
        return count == 0 ? 0.0 : static_cast<double>(value) / static_cast<double>(count);
    };

    const uint64_t classifier_kmers = classifier_thread_pool->classifier_processed_kmers.load();
    const uint64_t bloom_calls = classifier_thread_pool->bloom_calls.load();
    const uint64_t sampled_calls = classifier_thread_pool->bloom_sampled_calls.load();
    const uint64_t second_attempts = classifier_thread_pool->bloom_second_layer_attempts.load();
    const uint64_t sampled_second_attempts = classifier_thread_pool->bloom_sampled_second_layer_attempts.load();
    const uint64_t true_count = classifier_thread_pool->bloom_true_count.load();
    const uint64_t first_fetch_or_count = classifier_thread_pool->bloom_first_fetch_or_count.load();
    const uint64_t second_fetch_or_count = classifier_thread_pool->bloom_second_fetch_or_count.load();
    const uint64_t sampled_first_fetch_or_count = classifier_thread_pool->bloom_sampled_first_fetch_or_count.load();
    const uint64_t sampled_second_fetch_or_count = classifier_thread_pool->bloom_sampled_second_fetch_or_count.load();

    std::cout << "Classifier processed k-mers: " << classifier_kmers << std::endl;
    std::cout << "Prefix count cycles/kmer: "
              << per(classifier_thread_pool->prefix_count_cycles.load(), classifier_kmers) << std::endl;
    std::cout << "Prefix reorder cycles/kmer: "
              << per(classifier_thread_pool->prefix_reorder_cycles.load(), classifier_kmers) << std::endl;
    std::cout << "Classify local block cycles/kmer: "
              << per(classifier_thread_pool->classify_local_block_cycles.load(), classifier_kmers) << std::endl;
    std::cout << "Tree add cycles/kmer: "
              << per(classifier_thread_pool->tree_add_cycles.load(), classifier_kmers) << std::endl;

    std::cout << "Bloom calls: " << bloom_calls << std::endl;
    std::cout << "Bloom sampled calls: " << sampled_calls << std::endl;
    std::cout << "Bloom true ratio: " << per(true_count, bloom_calls) << std::endl;
    std::cout << "Bloom second-layer attempt ratio: "
              << per(second_attempts, bloom_calls) << std::endl;
    std::cout << "Bloom first fetch_or ratio: "
              << per(first_fetch_or_count, bloom_calls) << std::endl;
    std::cout << "Bloom second fetch_or ratio: "
              << per(second_fetch_or_count, second_attempts) << std::endl;

    std::cout << "Bloom total cycles/sample: "
              << per(classifier_thread_pool->bloom_total_cycles.load(), sampled_calls) << std::endl;
    std::cout << "Bloom hash cycles/sample: "
              << per(classifier_thread_pool->bloom_hash_cycles.load(), sampled_calls) << std::endl;
    std::cout << "Bloom index+mask cycles/sample: "
              << per(classifier_thread_pool->bloom_index_mask_cycles.load(), sampled_calls) << std::endl;
    std::cout << "Bloom first load cycles/sample: "
              << per(classifier_thread_pool->bloom_first_load_cycles.load(), sampled_calls) << std::endl;
    std::cout << "Bloom first filter cycles/sample: "
              << per(classifier_thread_pool->bloom_first_filter_cycles.load(), sampled_calls) << std::endl;
    std::cout << "Bloom first fetch_or cycles/first-fetch-sample: "
              << per(classifier_thread_pool->bloom_first_fetch_or_cycles.load(), sampled_first_fetch_or_count) << std::endl;
    std::cout << "Bloom second load cycles/second-attempt-sample: "
              << per(classifier_thread_pool->bloom_second_load_cycles.load(), sampled_second_attempts) << std::endl;
    std::cout << "Bloom second filter cycles/second-attempt-sample: "
              << per(classifier_thread_pool->bloom_second_filter_cycles.load(), sampled_second_attempts) << std::endl;
    std::cout << "Bloom second fetch_or cycles/second-fetch-sample: "
              << per(classifier_thread_pool->bloom_second_fetch_or_cycles.load(), sampled_second_fetch_or_count) << std::endl;

    /*for (uint64_t i = 0; i < (1ULL << (2 * ROOT_BASES)); ++i)
    {
        std::cout << "Root prefix " << std::format("{:0{}d}", i, 3) << " count: " << root_counts[i].load(std::memory_order_relaxed) << std::endl;
    }*/
#endif

    return 0;
}
