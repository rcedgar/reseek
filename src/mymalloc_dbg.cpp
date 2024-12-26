#include "myutils.h"

#define TEST_MYMALLOC_DBG	0

#if MYMALLOC_DBG

struct mem_data
	{
	mem_data *fwd;
	mem_data *bwd;
	uint32_t magic;
	const char *fn;
	int linenr;
	uint bytes;
	uint idx;
	uint secs;
	};

static mutex s_Lock;
static bool s_trace = false;
static FILE *s_ftrace;
static time_t s_t0 = time(0);
static const uint32_t MEM_MAGIC = 0x3e3e3e3e;
uint64 g_myalloc_netbytes;
atomic<int> g_mymalloc_netnewcount;
atomic<uint64> g_mymalloc_netnewbytes;
atomic<int> g_mymalloc_netallocount;
static uint s_nrblocks;
static uint s_summary_secs = 0;
static uint s_dump_secs = 0;
static time_t s_last_summary_secs = 0;
static FILE *s_summary_file = 0;
static FILE *s_state_file = 0;
static FILE *s_proc_maps_copy_file = 0;
static FILE *s_proc_maps_file = 0;
static bool s_write_state_done = false;

void mymalloc_set_dump_secs(uint secs)
	{
	s_dump_secs = secs;
	if (secs > 0 && s_state_file == 0)
		{
		s_state_file = fopen("myalloc.state", "w");
		s_proc_maps_file = fopen("/proc/self/maps", "r");
		s_proc_maps_copy_file = fopen("myalloc.procmaps", "w");
		setbuf(s_state_file, 0);
		setbuf(s_proc_maps_file, 0);
		}
	}

uint mymalloc_set_summary_secs(uint secs)
	{
	uint old_secs = s_summary_secs;
	if (secs > 0 && s_summary_file == 0)
		s_summary_file = fopen("myalloc.summary", "w");
	s_summary_secs = secs;
	return old_secs;
	}

static mem_data *s_BlockList;

static void InsertBlock(mem_data *Block)
	{
	if (s_BlockList != 0)
		s_BlockList->bwd = Block;
	Block->bwd = 0;
	Block->fwd = s_BlockList;
	s_BlockList = Block;
	}

static void DeleteBlock(mem_data *Block)
	{
	if (Block->bwd == 0)
		s_BlockList = Block->fwd;
	else
		Block->bwd->fwd = Block->fwd;

	if (Block->fwd != 0)
		Block->fwd->bwd = Block->bwd;
	}

static void mymalloc_trace_init()
	{
	s_Lock.lock();
	if (s_ftrace == 0)
		{
		s_ftrace = fopen("mem.trace", "w");
		setbuf(s_ftrace, 0);
		}
	s_Lock.unlock();
	}

void mymalloc_trace(bool On)
	{
	s_trace = On;
	if (On)
		mymalloc_trace_init();
	}

void *mymalloc_lo(uint n, uint bytes1, const char *fn, int linenr)
	{
	size_t tot = size_t(n)*size_t(bytes1) + sizeof(mem_data);

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

	void *userp = (void *) (((byte *) p) + sizeof(mem_data));
	s_Lock.lock();
	++g_mymalloc_netallocount;
	g_myalloc_netbytes += requested_bytes;
	if (s_trace)
		fprintf(s_ftrace, "%s:%d\t%p\t%u\n", fn, linenr, userp, requested_bytes);

	mem_data *Block = (mem_data *) p;

	++s_nrblocks;
	time_t t = time(0);
	uint secs = uint(t - s_t0);
	Block->magic = MEM_MAGIC;
	Block->fn = fn;
	Block->linenr = linenr;
	Block->bytes = requested_bytes;
	Block->idx = s_nrblocks;
	Block->secs = secs;

	InsertBlock(Block);

	bool do_summary =
		(s_summary_secs > 0 && secs - s_last_summary_secs >= s_summary_secs);
	if (do_summary)
		s_last_summary_secs = secs;
	bool do_write_state = (s_dump_secs > 0 && secs >= s_dump_secs);
	s_Lock.unlock();

	if (do_summary)
		mymalloc_print_summary(s_summary_file, "");

	if (do_write_state)
		mymalloc_write_state();

	return userp;
	}

void *mymalloc(uint n, uint bytes1, const char *fn, int linenr)
	{
	void *p = mymalloc_lo(n, bytes1, fn, linenr);
	return p;
	}

uint myfree_lo(void *userp)
	{
	if (userp == 0)
		return 0;

	mem_data *p = (mem_data *) (((byte *) userp) - sizeof(mem_data));

	if (p->magic != MEM_MAGIC)
		{
		fprintf(stderr, "\nmyfree MAGIC ERROR\n");
		exit(1);
		}

	s_Lock.lock();
	DeleteBlock(p);
	--g_mymalloc_netallocount;
	uint bytes = p->bytes;
	g_myalloc_netbytes -= bytes;
	if (s_trace)
		{
		uint secs = uint(time(0) - s_t0);
		fprintf(s_ftrace, "%u\tfree\t%p\n", secs, userp);
		}
	s_Lock.unlock();
	free(p);
	return bytes;
	}

void myfree(void *userp)
	{
	myfree_lo(userp);
	}

void *operator new(size_t n)
	{
	void *p = mymalloc_lo(uint(n), 1, "operator_new", 0);
	++g_mymalloc_netnewcount;
	g_mymalloc_netnewbytes += n;
	return p;
	}

void operator delete(void *p)
	{
	--g_mymalloc_netnewcount;
	uint bytes = myfree_lo(p);
	g_mymalloc_netnewbytes -= bytes;
	}

void mymalloc_write_state()
	{
	s_Lock.lock();
	if (s_write_state_done)
		{
		s_Lock.unlock();
		return;
		}
	fprintf(stderr, "\nmymalloc_write_state()\n");
	for (mem_data *Block = s_BlockList; Block != 0; Block = Block->fwd)
		fprintf(s_state_file, "%p\t%u\t%u\t%u\t%s:%d\n",
				Block,
				Block->idx, 
				Block->bytes,
				Block->secs,
				GetBaseName(Block->fn),
				Block->linenr);
	fflush(s_state_file);

	for (;;)
		{
		int c = fgetc(s_proc_maps_file);
		if (c == -1)
			break;
		fputc(c, s_proc_maps_copy_file);
		}
	fflush(s_proc_maps_copy_file);

	s_write_state_done = true;
	std::quick_exit(0);
	s_Lock.unlock();
	}

void mymalloc_print_summary(FILE *f, const char *Msg)
	{
	if (f == 0)
		return;
	s_Lock.lock();

	fprintf(f, "\n");
	if (strcmp(Msg, "") == 0)
		fprintf(f, "mymalloc_print_summary(%u)\n", uint(time(0)-s_t0));
	else
		fprintf(f, "mymalloc_print_summary(%s)\n", Msg);

	uint64 sumbytes = 0;
	uint nrblocks = 0;
	for (mem_data *Block = s_BlockList; Block != 0; Block = Block->fwd)
		{
		++nrblocks;
		sumbytes += Block->bytes;
		}

	double overhead = double(nrblocks)*sizeof(mem_data);

	double mem = GetMemUseBytes();
	double net = double(g_myalloc_netbytes);
	double sum = double(sumbytes);
	double newbytes = double(g_mymalloc_netnewbytes.load());

	fprintf(f, "%10.0f  mem_bytes (%.3g)\n", mem, mem);
	fprintf(f, "%10.0f  net_bytes (%.3g)\n", net, net);
	fprintf(f, "%10.0f  sum_block_bytes (%.3g)\n", sum, sum);
	fprintf(f, "%10.0f  net_new_bytes (%.3g)\n", newbytes, newbytes);
	fprintf(f, "%10.0f  overhead (%.3g)\n", overhead, overhead);
	fprintf(f, "%10u  net_allocs\n", g_mymalloc_netallocount.load());
	fprintf(f, "%10u  net_news\n", g_mymalloc_netnewcount.load());
	fprintf(f, "%10u  nrblocks\n", nrblocks);
	fflush(f);
	s_Lock.unlock();
	}

void mymalloc_print_state(FILE *f, const char *Msg)
	{
	if (f == 0)
		return;

	fprintf(f, "\n");
	fprintf(f, "mymalloc_print_state(%s)\n", Msg);
	fprintf(f, "    idx     bytes   secs  src\n");
	//          1234567  12345678  12345
	uint64 sumbytes = 0;
	s_Lock.lock();
	for (mem_data *Block = s_BlockList; Block != 0; Block = Block->fwd)
		fprintf(f, "%7u  %8u  %5u  %s:%d\n",
				Block->idx, 
				Block->bytes,
				Block->secs,
				GetBaseName(Block->fn),
				Block->linenr);

	s_Lock.unlock();
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
	mymalloc_write_state("mymalloc.state");
	}

#endif