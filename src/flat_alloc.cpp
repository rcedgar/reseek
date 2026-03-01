#include "flat_alloc.h"

#if defined(_MSC_VER)
    #include <malloc.h>   // _aligned_malloc/_aligned_free
#else
    #include <stdlib.h>   // posix_memalign, free
#endif

namespace flatmem
{
    void* aligned_malloc(std::size_t bytes, std::size_t alignment)
    {
        if (bytes == 0) bytes = 1;

#if defined(_MSC_VER)
        return _aligned_malloc(bytes, alignment);
#else
        void* p = 0;
        // posix_memalign requires alignment to be power-of-two and multiple of sizeof(void*)
        if (alignment < sizeof(void*)) alignment = sizeof(void*);
        int rc = posix_memalign(&p, alignment, bytes);
        return (rc == 0) ? p : 0;
#endif
    }

    void aligned_free(void* p)
    {
        if (!p) return;
#if defined(_MSC_VER)
        _aligned_free(p);
#else
        free(p);
#endif
    }
}