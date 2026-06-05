#ifndef FINAL_DRAIN_WRITER_HEADER
#define FINAL_DRAIN_WRITER_HEADER

#include "definition.h"

#include <string>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>

class FinalDrainWriter
{
private:
    int out_file = -1;
    std::vector<char> buffer;
    int buffer_offset = 0;
    uint64_t local_sorted_kmer_count = 0;

public:
    FinalDrainWriter(const FinalDrainWriter &) = delete;
    FinalDrainWriter &operator=(const FinalDrainWriter &) = delete;
    FinalDrainWriter(FinalDrainWriter &&) = delete;
    FinalDrainWriter &operator=(FinalDrainWriter &&) = delete;

    explicit FinalDrainWriter() : out_file(-1), buffer(DRAIN_EXPORT_BUFFER_SIZE), buffer_offset(0), local_sorted_kmer_count(0) {}

    void open(const uint32_t root_id)
    {
        buffer_offset = 0;
        std::string file_name = temp_dir + "root_" + std::to_string(root_id) + ".bin";
        out_file = ::open(file_name.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (out_file < 0)
        {
            std::cerr << "failed to open final drain output file" << std::endl;
            std::exit(-1);
        }
    }

    ~FinalDrainWriter()
    {
        if (out_file != -1)
        {
            ::close(out_file);
            out_file = -1;
        }
        if(local_sorted_kmer_count > 0)
        {
            sorted_kmer_count.fetch_add(local_sorted_kmer_count, std::memory_order_relaxed);
            local_sorted_kmer_count = 0;
        }
    }

    void close()
    {
        if (out_file != -1)
        {
            if (buffer_offset > 0)
            {
                ::write(out_file, buffer.data(), buffer_offset);
                buffer_offset = 0;
            }
            ::close(out_file);
            out_file = -1;
        }
        sorted_kmer_count.fetch_add(local_sorted_kmer_count, std::memory_order_relaxed);
        local_sorted_kmer_count = 0;
    }

    template <typename T>
    std::size_t write(const T &kmer_info)
    {
        //sorted_kmer_count.fetch_add(1, std::memory_order_relaxed);
        local_sorted_kmer_count++;
        if (sizeof(T) + buffer_offset > DRAIN_EXPORT_BUFFER_SIZE) [[unlikely]]
        {
            ::write(out_file, buffer.data(), buffer_offset);
            buffer_offset = 0;
        }
        std::memcpy(buffer.data() + buffer_offset, &kmer_info, sizeof(T));
        buffer_offset += sizeof(T);
        return sizeof(T);
    }
};

#endif