#include "myutils.h"

void *mymalloc(unsigned n, unsigned bytes1)
	{
	size_t tot = size_t(n)*size_t(bytes1);
	if (size_t(uint32(tot)) != tot)
		{
		fprintf(stderr, "\n\nmymalloc(%u, %u) overflow\n\n", n, bytes1);
		exit(1);
		}
	void *p = malloc(tot);
	if (0 == p)
		{
		double b = GetMemUseBytes();
		fprintf(stderr, "\n\nOut of memory mymalloc(%u, %u), curr %.3g bytes\n\n",
		  n, bytes1, b);
		exit(1);
		}
	return p;
	}

#if !MYMALLOC_DBG

void myfree(void *p)
	{
	if (p == 0)
		return;
	free(p);
	}

#endif // !MYMALLOC_DBG
