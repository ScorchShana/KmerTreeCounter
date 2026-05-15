#ifndef HASH_FUNCTION_HEADER
#define HASH_FUNCTION_HEADER

#include <cstdint>

#include "kmer.h"


template <uint32_t N>
inline uint64_t hash_func(const kmer<N> &key)
{
    uint64_t h = key.data[0];
    for (uint32_t i = 1; i < N; ++i)
    {
        h ^= key.data[i];
        h *= 0x9e3779b97f4a7c15ULL;
    }
    return h;
}

/*
template <uint32_t N>
inline uint64_t hash_func(const kmer<N> &key)
{
    uint64_t h = key.data[0];
    for (uint32_t i = 1; i < N; ++i)
    {
        h += key.data[i];
    }
    return h;
}
*/

/*
template <uint32_t N>
inline uint64_t hash_func(const kmer<N> &key)
{
    uint64_t h = XXH3_64bits(key.data.data(), sizeof(uint64_t) * N);
    return h;
}
*/
/*
inline uint64_t mum64(uint64_t A, uint64_t B)
{
    __uint128_t r = static_cast<__uint128_t>(A) * static_cast<__uint128_t>(B | 1ULL);
    return static_cast<uint64_t>(r ^ (r >> 64));
}

template <uint32_t N>
inline uint64_t hash_func(const kmer<N> &key)
{
    if constexpr (N == 1)
    {
        uint64_t h = key.data[0];
        h += 0x9e3779b97f4a7c15ULL;
        h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
        h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
        return h ^ (h >> 31);
    }
    if constexpr (N == 2)
    {
        uint64_t h = mum64(key.data[0] ^ 0x9e3779b97f4a7c15ULL, key.data[1] ^ 0x9e3779b97f4a7c15ULL);
        return h;
    }
    else
    {
        uint64_t h1 = mum64(key.data[0], key.data[1]);
        uint64_t h2 = mum64(key.data[2], key.data[3]);
        uint64_t h = mum64(h1 ^ 0x9e3779b97f4a7c15ULL, h2);
        return h;
    }
}
*/
#endif