#pragma once

void myfree(void *p);

#define MYMALLOC_DBG		1

void myfree(void *p);

#if MYMALLOC_DBG

void mymalloc_print_state(FILE *f, const char *Msg);
void mymalloc_trace(bool on);
void *mymalloc(unsigned n, unsigned bytes, const char *fn, int linenr);
void *mymalloc(unsigned n, unsigned bytes);
#define myalloc(t, n)	(t *) mymalloc((n), sizeof(t), __FILE__, __LINE__)

#else
#define mymalloc_print_state(f, Msg)	/* empty */
void *mymalloc(unsigned n, unsigned bytes);
#define myalloc(t, n)	(t *) mymalloc((n), sizeof(t))

#endif
