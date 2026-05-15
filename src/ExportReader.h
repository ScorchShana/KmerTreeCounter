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

template <uint32_t N>
class ExportReader
{

    std::FILE *file;
    uint64_t kmer_amount;

public:
    ExportReader(const ExportReader &) = delete;
    ExportReader &operator=(const ExportReader &) = delete;
    ExportReader(ExportReader &&) = delete;
    ExportReader &operator=(ExportReader &&) = delete;

    explicit ExportReader() : file(nullptr),kmer_amount(0)
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
        std::string filename = temp_dir + "low_" + std::to_string(prefix) + ".bin";
        file = std::fopen(filename.c_str(), "rb");
        if (!file) [[unlikely]]
        {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        struct stat st;
        if (fstat(fileno(file), &st) != 0) [[unlikely]]
        {
            std::fclose(file);
            file = nullptr;
            throw std::runtime_error("Failed to stat file: " + filename);
        }
        kmer_amount = st.st_size / sizeof(kmer<N>);
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

    uint64_t read_kmers(kmer<N> *buffer, uint64_t max_kmers_to_read)
    {
        uint64_t read_size = std::fread(buffer, sizeof(kmer<N>), max_kmers_to_read, file);
        uint64_t read_count = read_size / sizeof(kmer<N>);
        return read_count;
    }

};

#endif