#pragma once

void myfree(void *p);

#if _MSC_VER && DEBUG

void *mymalloc(unsigned n, unsigned bytes, const char *fn, int linenr);
#define myalloc(t, n)	(t *) mymalloc((n), sizeof(t), __FILE__, __LINE__)
void LogHeapSummary(const char *s);
void LogHeapBlocks(int Type, const char *s);

#else

void *mymalloc(unsigned n, unsigned bytes);
#define myalloc(t, n)	(t *) mymalloc((n), sizeof(t))
#define LogHeapSummary(s)	/* empty */
#define LogHeapBlocks(Type, s)	/* empty */

#endif
