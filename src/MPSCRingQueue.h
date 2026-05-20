#ifndef MPSC_RING_QUEUE_HEADER
#define MPSC_RING_QUEUE_HEADER

#include "definition.h"

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <utility>

template <typename T, std::size_t Capacity>
class MPSCRingQueue
{

	static_assert(Capacity > 1, "Capacity must be greater than 1");
	static_assert((Capacity & (Capacity - 1)) == 0, "Capacity must be a power of 2");
	static_assert(std::is_trivially_copyable_v<T>, "T must be trivially copyable");

	struct alignas(CACHE_LINE_SIZE) Cell
	{
		std::atomic<std::size_t> sequence;
		T data;
	};

	constexpr static std::size_t mask = Capacity - 1;

	constexpr static int SPIN_LIMIT = 256;
	constexpr static int BACKOFF_LIMIT = 64;

	std::atomic<std::size_t> enqueue_pos;
	std::atomic<std::size_t> dequeue_pos;
	std::array<Cell, Capacity> buffer;

public:
	MPSCRingQueue(const MPSCRingQueue &) = delete;
	MPSCRingQueue &operator=(const MPSCRingQueue &) = delete;

	MPSCRingQueue() : enqueue_pos(0), dequeue_pos(0)
	{
		for (size_t i = 0; i < Capacity; ++i)
		{
			buffer[i].sequence.store(i, std::memory_order_relaxed);
		}
	}

	bool try_enqueue(const T &item)
	{
		Cell *cell;
		std::size_t pos = enqueue_pos.load(std::memory_order_relaxed);

		cell = &buffer[pos & mask];

		std::size_t seq = cell->sequence.load(std::memory_order_acquire);
		int64_t diff = static_cast<int64_t>(seq - pos);

		if (diff == 0)
		{
			if (enqueue_pos.compare_exchange_weak(
					pos,
					pos + 1,
					std::memory_order_relaxed,
					std::memory_order_relaxed))
			{
				cell->data = item;
				// 发布数据
				cell->sequence.store(pos + 1, std::memory_order_release);

				return true;
			}
			// CAS 失败时，pos 会被 compare_exchange_weak 更新为当前值
		}
		else
		{
			// diff < 0 当前槽位还没被消费者释放，队列满
			// diff > 0 其他生产者已经推进 enqueue_pos
			return false;
		}

		return false;
	}

	void enqueue(const T &item)
	{
		int spin_count = 1;
		int backoff = 1;
		while (!try_enqueue(item))
		{
			if (spin_count >= BACKOFF_LIMIT) [[unlikely]]
			{
				std::this_thread::yield();
				spin_count = 0;
				backoff = 1;
			}
			else
			{
				for (int i = 0; i < backoff; ++i)
				{
					cpu_relax();
				}
				spin_count++;
				backoff = std::min(backoff * 2, BACKOFF_LIMIT);
			}
		}
	}

	bool try_dequeue(T &item)
	{
		Cell *cell;
		std::size_t pos = dequeue_pos.load(std::memory_order_relaxed);

		cell = &buffer[pos & mask];

		std::size_t seq = cell->sequence.load(std::memory_order_acquire);
		int64_t diff = static_cast<int64_t>(seq - (pos + 1));

		if (diff == 0)
		{
			dequeue_pos.store(pos + 1, std::memory_order_relaxed);
			item = cell->data;
			// 回收槽位，发布给生产者
			cell->sequence.store(pos + Capacity, std::memory_order_release);
			return true;

			// CAS 失败时，pos 会被 compare_exchange_weak 更新为当前值
		}
		else
		{
			// diff < 0 当前槽位还没有生产者发布数据，队列空
			// diff > 0 其他消费者已经推进 dequeue_pos
			return false;
		}

		return false;
	}

	void dequeue(T &item)
	{
		int spin_count = 1;
		int backoff = 1;
		while (!try_dequeue(item))
		{
			if (spin_count >= BACKOFF_LIMIT) [[unlikely]]
			{
				std::this_thread::yield();
				spin_count = 0;
				backoff = 1;
			}
			else
			{
				for (int i = 0; i < backoff; ++i)
				{
					cpu_relax();
				}
				spin_count++;
				backoff = std::min(backoff * 2, BACKOFF_LIMIT);
			}
		}
	}

	std::uint64_t size() const
	{
		std::size_t cur_enqueue_pos = enqueue_pos.load(std::memory_order_relaxed);
		std::size_t cur_dequeue_pos = dequeue_pos.load(std::memory_order_relaxed);

		return static_cast<std::uint64_t>(cur_enqueue_pos - cur_dequeue_pos);
	}

	bool empty() const
	{
		std::size_t cur_enqueue_pos = enqueue_pos.load(std::memory_order_relaxed);
		std::size_t cur_dequeue_pos = dequeue_pos.load(std::memory_order_relaxed);

		return cur_enqueue_pos == cur_dequeue_pos;
	}
};

#endif // MPMC_RING_QUEUE_HEADER
