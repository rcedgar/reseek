#include "myutils.h"

void *mymalloc(unsigned n, unsigned bytes1)
	{
	size_t tot = size_t(n)*size_t(bytes1);
	if (size_t(uint32(tot)) != tot)
		{
		fprintf(stderr, "\n\nmymalloc(%u, %u) overflow\n\n", n, bytes1);
		dbrk(1);
		exit(1);
		}
	void *p = malloc(tot);
	if (0 == p)
		{
		double b = GetMemUseBytes();
		fprintf(stderr, "\n\nOut of memory mymalloc(%u, %u), curr %.3g bytes\n\n",
		  n, bytes1, b);
		dbrk(1);
		exit(1);
		}
	return p;
	}

void *mymalloc64(uint64_t n, uint64_t bytes1)
	{
	size_t tot = size_t(n)*size_t(bytes1);
	if (tot/bytes1 != n || tot/n != bytes1)
		{
		fprintf(stderr, "\n\nmymalloc64(%.3g, %.3g) overflow\n\n",
				double(n), double(bytes1));
		dbrk(1);
		exit(1);
		}
	void *p = malloc(tot);
	if (0 == p)
		{
		double b = GetMemUseBytes();
		fprintf(stderr, "\n\nOut of memory mymalloc64(%.3g, %.3g), curr %.3g bytes\n\n",
		  double(n), double(bytes1), b);
		dbrk(1);
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
