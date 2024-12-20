#include "myutils.h"

#define TEST_MYMALLOC_DBG	1

#if MYMALLOC_DBG

struct MemData
	{
	MemData *fwd;
	MemData *bwd;
	uint32_t magic;
	const char *fn;
	int linenr;
	uint bytes;
	uint idx;
	uint secs;
	};

static mutex s_Lock;
static bool s_trace = false;
static bool s_use_malloc = false;
static FILE *s_ftrace;
static time_t s_t0 = time(0);
static const uint32_t MEM_MAGIC = 0x3e3e3e3e;
uint64 g_myalloc_netbytes;
atomic<int> g_mymalloc_netnewcount;
atomic<int> g_mymalloc_netallocount;
static uint s_nrblocks;

static MemData *s_BlockList;

static void InsertBlock(MemData *Block)
	{
	if (s_BlockList != 0)
		s_BlockList->bwd = Block;
	Block->bwd = 0;
	Block->fwd = s_BlockList;
	s_BlockList = Block;
	}

static void DeleteBlock(MemData *Block)
	{
	if (Block->bwd == 0)
		s_BlockList = Block->fwd;
	else
		Block->bwd->fwd = Block->fwd;

	if (Block->fwd != 0)
		Block->fwd->bwd = Block->bwd;
	}

static void MemTraceInit()
	{
	s_Lock.lock();
	if (s_ftrace == 0)
		{
		s_ftrace = fopen("mem.trace", "w");
		setbuf(s_ftrace, 0);
		}
	s_Lock.unlock();
	}

void MemTrace(bool On)
	{
	s_trace = On;
	if (On)
		MemTraceInit();
	}

void *mymalloc(uint n, uint bytes1, const char *fn, int linenr)
	{
	if (s_use_malloc)
		return mymalloc(n, bytes1);

	size_t tot = size_t(n)*size_t(bytes1) + sizeof(MemData);
	if (size_t(uint32(tot)) != tot)
		{
		fprintf(stderr, "\n\nmymalloc(%u, %u) overflow\n\n", n, bytes1);
		exit(1);
		}
	uint32_t requested_bytes = uint32_t(n*bytes1);
	void *p = malloc(tot);
	if (0 == p)
		{
		double b = GetMemUseBytes();
		fprintf(stderr, "\n\nOut of memory mymalloc(%u, %u), curr %.3g bytes\n\n",
		  n, bytes1, b);
		exit(1);
		}

	void *userp = (void *) (((byte *) p) + sizeof(MemData));
	s_Lock.lock();
	++g_mymalloc_netallocount;
	g_myalloc_netbytes += requested_bytes;
	if (s_trace)
		fprintf(s_ftrace, "%s:%d\t%p\t%u\n", fn, linenr, userp, requested_bytes);
	s_Lock.unlock();

	MemData *Block = (MemData *) p;

	++s_nrblocks;
	Block->magic = MEM_MAGIC;
	Block->fn = fn;
	Block->linenr = linenr;
	Block->bytes = requested_bytes;
	Block->idx = s_nrblocks;
	Block->secs = uint(time(0) - s_t0);

	InsertBlock(Block);

	return userp;
	}

void myfree(void *userp)
	{
	if (userp == 0)
		return;
	if (s_use_malloc)
		{
		free(userp);
		return;
		}

	MemData *p = (MemData *) (((byte *) userp) - sizeof(MemData));

	if (p->magic != MEM_MAGIC)
		{
		fprintf(stderr, "\nmyfree MAGIC ERROR\n");
		exit(1);
		}

	DeleteBlock(p);

	s_Lock.lock();
	--g_mymalloc_netallocount;
	g_myalloc_netbytes -= p->bytes;
	if (s_trace)
		{
		uint secs = uint(time(0) - s_t0);
		fprintf(s_ftrace, "%u\tfree\t%p\n", secs, userp);
		}
	s_Lock.unlock();
	free(p);
	}

void *operator new(size_t n)
	{
	void *p = mymalloc(uint(n), 1, "operator_new", 0);
	++g_mymalloc_netnewcount;
	return p;
	}

void operator delete(void *p)
	{
	--g_mymalloc_netnewcount;
	myfree(p);
	}

void mymalloc_print_state(FILE *f, const char *Msg)
	{
	if (f == 0)
		return;
	const bool saved_use_malloc = s_use_malloc;
	s_use_malloc = true;

	fprintf(f, "\n");
	fprintf(f, "mymalloc_print_state(%s)\n", Msg);
	fprintf(f, "    idx     bytes   secs  src\n");
	//          1234567  12345678  12345
	uint64 sumbytes = 0;
	for (MemData *Block = s_BlockList; Block != 0; Block = Block->fwd)
		{
		sumbytes += Block->bytes;
		fprintf(f, "%7u  %8u  %5u  %s:%d\n",
				Block->idx, 
				Block->bytes,
				Block->secs,
				GetBaseName(Block->fn),
				Block->linenr);
		}

	fprintf(f, "%10.0f  mem bytes\n", GetMemUseBytes());
	fprintf(f, "%10.0f  net bytes\n", double(g_myalloc_netbytes));
	fprintf(f, "%10.0f  sum bytes\n", double(sumbytes));
	fprintf(f, "%10u  net allocs\n", g_mymalloc_netallocount.load());
	fprintf(f, "%10u  net news\n", g_mymalloc_netnewcount.load());

	s_use_malloc = saved_use_malloc;
	}

#endif // MYMALLOC_DBG

#if TEST_MYMALLOC_DBG

void cmd_test()
	{
	vector<byte *> ps;
	for (int i = 0; i < 100; ++i)
		{
		byte *p = myalloc(byte, randu32()%100 + 3);
		ps.push_back(p);
		}
	mymalloc_print_state(g_fLog, "");
	}

#endif