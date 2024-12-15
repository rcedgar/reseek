#include "myutils.h"

#if defined(__GNUC__) && DEBUG

static mutex s_Lock;
bool s_trace = false;

static FILE *f;
static time_t t0;

void MemTrace(bool On)
	{
	s_trace = On;
	}

static void init()
	{
	if (f != 0)
		return;
	f = fopen("myalloc.trace", "w");
	t0 = time(0);
	setbuf(f, 0);
	}

void *operator new(size_t n)
	{
	void *p = malloc(n);
	if (s_trace)
		{
		s_Lock.lock();
		init();
		uint secs = uint(time(0) - t0);
		fprintf(f, "%u\tnew\t%u\t%p\n", secs, uint(n), p);
		s_Lock.unlock();
		}
	return p;
	}

void operator delete(void *p)
	{
	if (s_trace)
		{
		s_Lock.lock();
		init();
		uint secs = uint(time(0) - t0);
		fprintf(f, "%u\tdelete\t%p\n", secs, p);
		s_Lock.unlock();
		}
	free(p);
	}

void *mymalloc(unsigned n, unsigned bytes, const char *fn, int linenr)
	{
	if (bytes == 0)
		return 0;
	size_t B = size_t(n)*size_t(bytes);
	void *p = malloc(B);
	if (0 == p)
		{
		double b = GetMemUseBytes();
		fprintf(stderr, "\n\n%s:%d Out of memory mymalloc(%u, %u), curr %.3g bytes\n\n",
		  fn, linenr, n, bytes, b);
		exit(1);
		}

	if (s_trace)
		{
		s_Lock.lock();
		init();
		uint secs = uint(time(0) - t0);
		fprintf(f, "%u\talloc\t%s\t%d\t%d\t%u\t%p\n", secs, fn, linenr, n, bytes, p);
		s_Lock.unlock();
		}

	return p;
	}

void myfree(void *p)
	{
	if (p == 0)
		return;
	if (s_trace)
		{
		s_Lock.lock();
		init();
		uint secs = uint(time(0) - t0);
		fprintf(f, "%u\tfree\t%p\n", secs, p);
		s_Lock.unlock();
		}
	free(p);
	}

void LogHeapSummary(const char *s)
	{
	}

void LogHeapBlocks(int aType, const char *s)
	{
	}

#endif // DEBUG
