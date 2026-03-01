#pragma once
#include <cstddef>

namespace flatmem
{
    // Returns aligned, uninitialized memory. Throws nothing; returns 0 on failure.
    void* aligned_malloc(std::size_t bytes, std::size_t alignment);
    void  aligned_free(void* p);
}