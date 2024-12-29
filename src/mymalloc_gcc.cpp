#if defined(__GNUC__) && !defined(__APPLE__)
#include <dlfcn.h>
#include <atomic>
#include <stdio.h>
#include <malloc.h>

#if 0

std::atomic<int> g_nrmallocs;
std::atomic<int> g_nrfrees;

typedef void (*ptr_freefunc)(void *ptr);
typedef void *(*ptr_mallocfunc)(size_t);

void *malloc(size_t size)
	{
	printf("malloc\n");
	ptr_mallocfunc mf = (ptr_mallocfunc) dlsym(RTLD_NEXT, "malloc");
	++g_nrmallocs;
	return mf(size);
	}

void free(void *ptr)
	{
	printf("free\n");
	ptr_freefunc ff = (ptr_freefunc) dlsym(RTLD_NEXT, "free");
	++g_nrfrees;
	ff(ptr);
	}
#endif

#if 0
// https://stackoverflow.com/questions/262439/create-a-wrapper-function-for-malloc-and-free-in-c
/* 
 * Link-time interposition of malloc and free using the static
 * linker's (ld) "--wrap symbol" flag.
 * 
 * Compile the executable using "-Wl,--wrap,malloc -Wl,--wrap,free".
 * This tells the linker to resolve references to malloc as
 * __wrap_malloc, free as __wrap_free, __real_malloc as malloc, and
 * __real_free as free.
 */
#include <stdio.h>

extern "C" {
void *__real_malloc(size_t size);
void __real_free(void *ptr);


/* 
 * __wrap_malloc - malloc wrapper function 
 */
void *__wrap_malloc(size_t size)
{
    void *ptr = __real_malloc(size);
    printf("malloc(%d) = %p\n", int(size), ptr);
    return ptr;
}

/* 
 * __wrap_free - free wrapper function 
 */
void __wrap_free(void *ptr)
{
    __real_free(ptr);
    printf("free(%p)\n", ptr);
}
};
#endif

#if 0
// https://www.gnu.org/software/libc/manual/html_node/Hooks-for-Malloc.html
// *****************************************
// HOOKS NOT SUPPORTED IN NEWER GCC RELEASES
// *****************************************
void *my_malloc_hook(size_t size, const void *caller)
	{
	printf("my_malloc_hook()\n");
	}

void my_free_hook(void *ptr, const void *caller)
	{
	printf("my_free_hook()\n");
	}

void my_hooks_init()
	{
	printf("my_hooks_init\n");
	__malloc_hook = my_malloc_hook;
	__free_hook = my_free_hook;
	}

#endif 

#endif