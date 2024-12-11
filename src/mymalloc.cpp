#include "myutils.h"

#if !DEBUG

void *mymalloc(unsigned n, unsigned bytes)
	{
	size_t B = size_t(n)*size_t(bytes);
	void *p = malloc(B);
	if (0 == p)
		{
		double b = GetMemUseBytes();
		fprintf(stderr, "\n\nOut of memory mymalloc(%u, %u), curr %.3g bytes\n\n",
		  n, bytes, b);
		exit(1);
		}
	return p;
	}

void myfree(void *p)
	{
	if (p == 0)
		return;
	free(p);
	}

#endif
