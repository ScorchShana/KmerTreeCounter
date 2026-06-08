#ifndef EXPORT_READER_HEADER
#define EXPORT_READER_HEADER

#include "definition.h"
#include "kmer.h"

#include <cstdio>
#include <array>
#include <string>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <iostream>

template <uint32_t N>
class ExportReader
{
    std::FILE* file;
    uint64_t kmer_amount;

public:
    ExportReader(const ExportReader&) = delete;
    ExportReader& operator=(const ExportReader&) = delete;
    ExportReader(ExportReader&&) = delete;
    ExportReader& operator=(ExportReader&&) = delete;

    explicit ExportReader() : file(nullptr), kmer_amount(0)
    {
    }

    ~ExportReader()
    {
        if (file)
        {
            std::fclose(file);
            file = nullptr;
        }
    }

    void open(const uint64_t prefix)
    {
        const std::string filename = temp_dir + "low_" + std::to_string(prefix) + ".bin";
        file = std::fopen(filename.c_str(), "rb");
        if (!file) [[unlikely]]
        {
            std::cerr << "Failed to open file: " << filename << std::endl;
            std::exit(-1);
        }

        // 读取 4KB header，获取实际 kmer 数量
        alignas(PAGE_SIZE) char header[PAGE_SIZE];
        size_t readn = std::fread(header, 1, PAGE_SIZE, file);
        if (readn == PAGE_SIZE)
        {
            kmer_amount = *reinterpret_cast<uint64_t*>(header);
        }
        else
        {
            kmer_amount = 0;
        }
    }

    void close()
    {
        if (file)
        {
            std::fclose(file);
            file = nullptr;
        }
    }

    uint64_t get_kmer_amount() const
    {
        return kmer_amount;
    }

    uint64_t read_kmers(kmer<N>* buffer, uint64_t max_kmers_to_read)
    {
        uint64_t read_count = std::fread(buffer, sizeof(kmer<N>), max_kmers_to_read, file);
        return read_count;
    }
};

#endif
