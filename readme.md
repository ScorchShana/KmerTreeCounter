# tree_v4 — High-Performance K-mer Counter

A header-only C++20, template-driven, lock-free k-mer counter for FASTQ files. Uses x86-64 SIMD (AVX2/SSE4.2), NUMA-aware memory allocation with huge pages, and an adaptive multi-threaded pipeline.

## Pipeline Architecture

```
FASTQ file
    │
    ├── [Stage 0: Pre-Read Sampling] ────────────────────────────────────┐
    │   FastqPreReader (main thread) + FastqParser (1 thread)            │
    │   + FastqPrefixCounter (up to 8 threads)                           │
    │   Reads ~128MB, counts per-root-prefix frequencies (256 buckets),  │
    │   determines hot/cold prefixes for bloom filter sizing             │
    └────────────────────────────────────────────────────────────────────┘
    │
    ▼ [FastqReader] (1 thread — main thread)
    reads FASTQ in chunks, parses Header→Sequence→Plus→Quality records,
    pushes raw data blocks
    ▼ SPMCRingMemoryPool (reader_parser_ring_pool)
[FastqParser × P] (P = max(1, n_thread/8) threads via ParserThreadPool)
    extracts canonical k-mers via SIMD batch GetKmer (no bloom check),
    pushes raw k-mer arrays
    ▼ RingMemoryPool (parser_classifier_ring_pool)
[FastqClassifier × C] (C threads via ClassifierThreadPool)
    groups k-mers by root prefix (256 buckets, first 4 bases),
    checks each k-mer against per-root ConcurrentDoubleBloomFilter
    ├─ First occurrence (bloom = true) → export_block_ptr → export_ring_pool
    └─ Repeated (bloom = false) → local_root_nodes → periodically flushed to KmerTree
    │   (overflow tasks enqueued to LayerQueues depth 0)
    ▼ LayerQueues — 4 × MPMCRingQueue<Task, 65536>
[KmerTree — 4-level radix tree]  +  [SchedulerThreadPool] (1 scheduler + W workers)
    ├─ Root layer (256 nodes, indexed by first 4 bases, 2-bit/base)
    ├─ Level 1   (16 children/node, indexed by next 2 bases)
    ├─ Level 2   (16 children/node, indexed by next 2 bases)
    └─ Level 3   (16 children/node, leaf) → ConcurrentCountingHashMap per node
    │
    │   Scheduler: monitors queue pressure across 4 depth queues,
    │   adaptively assigns/migrates worker threads via pressure scores
    │   Workers: dequeue tasks from assigned depth, route k-mers to child nodes,
    │   lazily create CountingHashMaps at leaf level
    │
    ├──→ [export_ring_pool → ExportWriter] (1 thread)
    │       writes first-occurrence k-mers to tmp/low.bin
    │
    └──→ [Final Drain] (parallel phase, n_thread threads)
            DFS traversal of all 256 root subtrees
            ├─ Leaf with hash map → insert remaining blocks
            └─ Leaf without hash map → sort + merge + write to tmp/root_{id}.bin
```

### Thread Allocation

```
parser_num    = max(1, n_thread / 8)
worker_budget = n_thread - 2 - parser_num          // minus reader + export
classifier_num  = worker_budget / (1 + TASK_CLASSIFIER_RATIO)   // RATIO = 1
tasker_num    = worker_budget - classifier_num
```

Total = parser_num + classifier_num + tasker_num + 1 (reader) + 1 (export writer) = n_thread.

## Component Details

| File | Role |
|------|------|
| `main.cpp` | Entry point, CLI argument parsing, pipeline orchestration, pre-read sampling, timing measurement |
| `definition.h` | All compile-time constants (pool sizes, queue capacities, k-mer params) and shared type definitions (`Task`, `KmerBatch`, `ExportBlock`, etc.) |
| `kmer.h` | `kmer<N>` — bit-packed 2-bit/base k-mer stored in N `uint64_t` words; `kmer_block<N>` — 4KB-aligned fixed-capacity batch of k-mers |
| `GetKmer.h` | `GetKmer<N>` — sliding-window k-mer extraction from raw DNA stream. Maintains both forward and reverse complement simultaneously using bit manipulation, with AVX2/SSE4.2 batch processing |
| `FastqReader.h` | `FastqReader<N>` — single-threaded FASTQ state-machine parser (Header → Sequence → Plus → Quality). Uses `posix_fadvise` for sequential readahead, handles cross-block k-mer overlap |
| `FastqPreReader.h` | `FastqPreReader<N>` — pre-read phase FASTQ reader, identical state machine to `FastqReader`, reads first ~128MB to gather prefix frequency statistics |
| `FastqParser.h` | `FastqParser<N>` — pure k-mer extractor. Consumes raw sequence blocks from the reader, extracts canonical k-mers via `GetKmer`, pushes raw k-mer arrays to the classifier (no bloom filter check) |
| `ParserThreadPool.h` | `ParserThreadPool<N>` — manages P `FastqParser` threads, each consuming from the shared `SPMCRingMemoryPool` |
| `FastqClassifier.h` | `FastqClassifier<N>` — bloom filter classifier. Groups k-mers by root prefix (256 buckets), checks each against per-root `ConcurrentDoubleBloomFilter`. First-occurrence k-mers are enqueued to `export_ring_pool`; repeated k-mers accumulate in per-thread local root nodes and flush to `KmerTree` |
| `ClassifierThreadPool.h` | `ClassifierThreadPool<N>` — manages C `FastqClassifier` threads, coordinates flush and producer-done signaling to the scheduler |
| `FastqPrefixCounter.h` | `FastqPrefixCounter<N>` — consumes k-mer arrays from the pre-read phase, counts how many k-mers fall into each of the 256 root prefixes. Results determine hot/cold bloom filter partitioning |
| `NewKmerTree.h` | `KmerTree<N>` — the core 4-level radix tree. Provides `main_add_kmer_block_with_local_root_nodes()` for classifier insertion via local root nodes, `thread_add_kmer()` for worker routing, `insert_kmer_in_task_to_node_hash_map()` for leaf counting, `final_drain_parallel()` for cleanup, and per-root `ConcurrentDoubleBloomFilter` instances for first-occurrence detection |
| `LayerQueues.h` | `LayerQueues<N>` — array of 4 `MPMCRingQueue<Task>` instances (one per tree depth 0–3) plus a final-drain queue, with atomic size tracking for the scheduler |
| `SchedulerThreadPool.h` | `SchedulerThreadPool<N>` — adaptive scheduling engine. 1 scheduler thread computes score-based pressure metrics (fill ratio, burst, trend, upstream pressure, empty penalty) to dynamically assign/migrate worker threads across queue depths; workers dequeue and process tasks, with thread-local task stack spillover |
| `SpinLock.h` | `SpinLock` — TATAS (test-and-test-and-set) spinlock with exponential backoff, yield-after-spin, and optional spin-loop counting for testing |
| `RingMemoryPool.h` | `RingMemoryPool<C>` — bipartite lock-free ring buffer with two MPMC queues (producer→consumer + consumer→producer), enabling fixed-size block reuse |
| `SPMCRingMemoryPool.h` | (in `RingMemoryPool.h`) — single-producer variant using SPSC + MPSC internal queues for the reader→parser path |
| `ConcurrentMemoryPool.h` | `ConcurrentMemoryPool` — NUMA-aware 4KB-block allocator. Per-node arenas, thread-local caching (local free stack + remote lists), 2MB huge page support with transparent huge page fallback |
| `MPMCRingQueue.h` | `MPMCRingQueue<T,C>` — lock-free multi-producer multi-consumer ring queue with a 4-state slot machine (EMPTY → STORING → STORED → LOADING) |
| `BloomFilter.h` | `ConcurrentDoubleBloomFilter<N>` — thread-safe double-layer bloom filter using atomic `fetch_or` for concurrent insertion and duplicate detection |
| `ConcurrentCountingHashMap.h` | `ConcurrentCountingHashMap<N>` — SIMD-accelerated (AVX2/SSE4.2) concurrent hash map with 4KB page-aligned NodeBlocks, lock-free probing with lock-based chain extension |
| `CountingHashTable.h` | `CountingHashTable<N,B,T>` — single-threaded open-addressing hash table with SIMD-accelerated control-byte matching, used as a thread-local accumulator before flushing to the global map |
| `ConcurrentMap.h` | `ConcurrentMap<N>` — lock-free concurrent hash map using linked-list buckets with CAS-based insertion |
| `ExportWriter.h` | `ExportWriter<N>` — single-threaded writer consuming from the export ring pool, writes first-occurrence k-mers to `tmp/low.bin` |
| `FinalDrainWriter.h` | `FinalDrainWriter` — per-root file writer for the final drain phase, writes sorted k-mer records to `tmp/root_{id}.bin` |
| `ExportReader.h` | `ExportReader<N>` — reads back exported `low.bin` files (for post-processing or analysis) |
| `HashFunction.h` | `hash_func<N>()` — simple hash utility for k-mers |
| `SplitMix.h` | `SplitMix64` — thread-safe split-mix pseudo-random number generator |
| `FixedStack.h` | `FixedStack<T,N>` — fixed-capacity stack (used for thread-local task stacks) |
| `FixedMinHeap.h` | `FixedMinHeap<T,C>` — fixed-capacity min-heap (available utility) |

## Key Design Highlights

- **Lock-free data flow**: MPMC ring queues for all cross-thread communication; atomic bloom filter and CAS-based hash map insertion; no global pipeline locks.
- **NUMA-aware memory**: Per-NUMA-node arenas with first-touch initialization and thread-local caching; automatic remote-list batching for cross-node frees.
- **Cache-line alignment**: All hot structures (`node`, `Task`, `MPMCRingQueue` head/tail, `Arena`, `SpinLock` counters) padded to 64-byte cache lines to prevent false sharing.
- **Adaptive scheduling**: Scheduler computes per-depth pressure scores (fill ratio, burst, trend, upstream pressure) and dynamically reassigns worker threads to balance load across tree depths.
- **Two-stage classification**: `FastqParser` is a pure k-mer extractor; a separate `FastqClassifier` stage checks bloom filters and splits first-occurrence vs. repeated k-mers, reducing per-parser overhead and improving cache locality.
- **Pre-read sampling**: A warmup phase reads ~128MB to profile root-prefix frequencies, enabling hot/cold bloom filter partitioning (hot prefixes get 2x-capacity dedicated filters; cool prefixes are grouped into shared partitions).
- **Early duplicate detection with two-path export**: Classifiers check each k-mer against a per-root `ConcurrentDoubleBloomFilter`. First occurrences are immediately exported as singletons; repeated k-mers flow into the radix tree and are counted in leaf-level `ConcurrentCountingHashMap` instances.
- **Thread-local aggregation**: Classifiers accumulate k-mers in per-thread local root node arrays (reducing tree lock contention); workers use local hash tables and local task stacks (reducing queue pressure).

## Build & Usage

### Dependencies
- C++20 compiler (GCC 11+ / Clang 14+)
- pthread, libnuma, libaio

### Build
```bash
cmake -B build
cmake --build build
```

### Run
```bash
./Tree <fastq_file> <k_len> <n_thread> <memory_limit_gb> [map_capacity] [min_count] [max_count] [parser_threads]
```

### Parameters
| Parameter | Description |
|-----------|-------------|
| `fastq_file` | Input FASTQ file path |
| `k_len` | K-mer length (≤ 128) |
| `n_thread` | Total thread count (≥ 6) |
| `memory_limit_gb` | Memory budget in GB |
| `map_capacity` | Hash map bucket capacity (default: 1024) |
| `min_count` | Minimum count threshold for export (default: 1) |
| `max_count` | Maximum count threshold for export (default: unlimited) |
| `parser_threads` | Override parser thread count (not yet implemented) |
