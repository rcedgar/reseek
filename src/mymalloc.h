#pragma once

void myfree(void *p);

#define MYMALLOC_TRACK_SIZE	0
#define MYMALLOC_TRACE		0

#if MYMALLOC_TRACK_SIZE
double GetMymallocTotal();
#endif

#if DEBUG

void *mymalloc(unsigned n, unsigned bytes, const char *fn, int linenr);
#define myalloc(t, n)	(t *) mymalloc((n), sizeof(t), __FILE__, __LINE__)
void LogHeapSummary(const char *s);
void LogHeapBlocks(int Type, const char *fn);
void MemTrace(bool on);

#else

void *mymalloc(unsigned n, unsigned bytes);
#define myalloc(t, n)	(t *) mymalloc((n), sizeof(t))
#define LogHeapSummary(s)	/* empty */
#define LogHeapBlocks(Type, fn)	/* empty */
#define MemTrace(on) /* empty */
#endif
