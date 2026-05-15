#ifndef FINAL_DRAIN_WRITER_HEADER
#define FINAL_DRAIN_WRITER_HEADER

#include "definition.h"

#include <cstdio>
#include <string>
#include <stdexcept>
#include <type_traits>

class FinalDrainWriter
{
private:
    std::FILE *out_file = nullptr;

public:
    FinalDrainWriter(const FinalDrainWriter &) = delete;
    FinalDrainWriter &operator=(const FinalDrainWriter &) = delete;
    FinalDrainWriter(FinalDrainWriter &&) = delete;
    FinalDrainWriter &operator=(FinalDrainWriter &&) = delete;

    explicit FinalDrainWriter() : out_file(nullptr) {}

    void open(const uint32_t root_id)
    {
        constexpr uint64_t kFileBufferSize = DRAIN_EXPORT_BUFFER_SIZE * 4;
        std::string file_name = temp_dir + "root_" + std::to_string(root_id) + ".bin";
        out_file = std::fopen(file_name.c_str(), "wb+");
        if (out_file == nullptr)
        {
            throw std::runtime_error("failed to open final drain output file");
        }
        if (std::setvbuf(out_file, nullptr, _IOFBF, kFileBufferSize) != 0)
        {
            std::fclose(out_file);
            throw std::runtime_error("failed to set buffer size for final drain output file");
        }
    }

    ~FinalDrainWriter()
    {
        if (out_file != nullptr)
        {
            std::fclose(out_file);
            out_file = nullptr;
        }
    }

    void close()
    {
        if (out_file != nullptr)
        {
            std::fclose(out_file);
            out_file = nullptr;
        }
    }

    template <typename T>
    std::size_t write(const T &kmer_info)
    {
        sorted_kmer_count.fetch_add(1, std::memory_order_relaxed);
        return std::fwrite(reinterpret_cast<const void *>(&kmer_info), sizeof(T), 1, out_file);
    }
};

#endif