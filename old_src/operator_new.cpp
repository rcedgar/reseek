#include "myutils.h"

#if 0
static atomic<size_t> allocation_count = 0;

void* operator new(std::size_t size) {
    allocation_count++;
    
    void* p = std::malloc(size);
    
    if (!p) throw std::bad_alloc();
    return p;
}

void operator delete(void* p) noexcept {
    if (p) {
        allocation_count--;
        std::free(p);
    }
}
#endif