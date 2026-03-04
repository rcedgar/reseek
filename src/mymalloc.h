#pragma once

void myfree(void *p);
void myfreea(void *p);

#define MYMALLOC_DBG		0

#if MYMALLOC_DBG

uint mymalloc_set_summary_secs(uint secs);
void mymalloc_set_dump_secs(uint secs);
void mymalloc_print_summary(const char *Msg = "");
void mymalloc_print_summary(int fd, const char *Msg);
void mymalloc_print_state(int fd, const char *Msg);
void mymalloc_write_state(const char *fn);
void mymalloc_write_map(const char *fn);
void mymalloc_trace(bool on);
void *mymalloc(unsigned n, unsigned bytes, const char *fn, int linenr);
void *mymalloc(unsigned n, unsigned bytes);
#define myalloc(t, n)	(t *) mymalloc((n), sizeof(t), __FILE__, __LINE__)

#else

#define mymalloc_print_summary(x)	0
#define mymalloc_set_summary_secs(secs) 0
#define mymalloc_set_dump_secs(secs) 0
#define mymalloc_print_state(f, Msg)	0
#define mymalloc_write_map(fn)	0
#define mymalloc_write_state(fn)	0
#define mymalloc_trace(on)	0
void *mymalloc(unsigned n, unsigned bytes);
void *mymalloca(unsigned n, unsigned bytes);
void *mymalloc64(uint64_t n, uint64_t bytes);
void *mymalloc64a(uint64_t n, uint64_t bytes);
#define myalloc(t, n)	(t *) mymalloc((n), sizeof(t))
#define myalloca(t, n)	(t *) mymalloca((n), sizeof(t))
#define myalloc64(t, n)	(t *) mymalloc64((n), sizeof(t))
#define myalloc64a(t, n)	(t *) mymalloc64a((n), sizeof(t))

#endif
