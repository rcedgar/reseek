#include "myutils.h"

#if !DEBUG

#if MYMALLOC_TRACK_SIZE

static mutex s_Lock;
static uint64 s_Total;

double GetMymallocTotal()
	{
	return double(s_Total);
	}

void *mymalloc(unsigned n, unsigned bytes)
	{
	size_t B = size_t(n)*size_t(bytes) + sizeof(size_t);
	uint Bu = uint(B);
	if (size_t(Bu) != B)
		{
		fprintf(stderr, "\n\nmymalloc(%u, %u) overflow\n\n", n, bytes);
		exit(1);
		}
	void *p = malloc(B);
	if (0 == p)
		{
		double b = GetMemUseBytes();
		fprintf(stderr, "\n\nOut of memory mymalloc(%u, %u), curr %.3g bytes\n\n",
		  n, bytes, b);
		exit(1);
		}
	size_t *sp = (size_t *) p;
	*sp = B;
	void *effp = (void *) (sp + 1);
	s_Lock.lock();
	s_Total += B;
	s_Lock.unlock();
	return (void *) effp;
	}

void myfree(void *p)
	{
	if (p == 0)
		return;
	size_t *sp = (size_t *) p;
	size_t B = sp[-1];
	s_Lock.lock();
	s_Total -= B;
	s_Lock.unlock();
	void *origp = (void *) (sp - 1);
	free(origp);
	}

void *operator new(size_t n)
	{
	void *p = mymalloc(uint(n), 1);
	return p;
	}

void operator delete(void *p)
	{
	myfree(p);
	}


#else

void *mymalloc(unsigned n, unsigned bytes)
	{
	size_t B = size_t(n)*size_t(bytes);
	uint Bu = uint(B);
	if (size_t(Bu) != B)
		{
		fprintf(stderr, "\n\nmymalloc(%u, %u) overflow\n\n", n, bytes);
		exit(1);
		}
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
#endif
