#include "myutils.h"
#include <fcntl.h>
#include <io.h>

#define TEST_MYMALLOC_DBG	1

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
static int s_trace_fd = -1;
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
static int s_summary_fd = -1;
static bool s_write_state_done = false;

static const size_t s_printf_buffer_sz = 1023;
static char s_printf_buffer[s_printf_buffer_sz+1];

static const char *ps(const char *Format, ...)
	{
	va_list ArgList;
	va_start(ArgList, Format);
	vsnprintf(s_printf_buffer, s_printf_buffer_sz, Format, ArgList);
	s_printf_buffer[s_printf_buffer_sz] = 0;
	va_end(ArgList);
	}

static void writes(int fd, const char *Format, ...)
	{
	va_list ArgList;
	va_start(ArgList, Format);
	vsnprintf(s_printf_buffer, s_printf_buffer_sz, Format, ArgList);
	s_printf_buffer[s_printf_buffer_sz] = 0;
	va_end(ArgList);
	size_t n = strlen(s_printf_buffer);
	write(fd, s_printf_buffer, uint(n));
	}

void mymalloc_set_dump_secs(uint secs)
	{
	s_dump_secs = secs;
	}

uint mymalloc_set_summary_secs(uint secs)
	{
	uint old_secs = s_summary_secs;
	if (secs > 0 && s_summary_fd < 0)
		s_summary_fd = open("myalloc.summary", O_WRONLY | O_CREAT, S_IREAD | S_IWRITE);
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
	if (s_trace_fd < 0)
		s_trace_fd = open("mem.trace", O_WRONLY | O_CREAT, S_IREAD | S_IWRITE);
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
		writes(2, "\n\nmymalloc(%u, %u) overflow\n\n", n, bytes1);
		exit(1);
		}
	uint32_t requested_bytes = uint32_t(n*bytes1);
	void *p = malloc(tot);
	if (0 == p)
		{
		double b = GetMemUseBytes();
		writes(2, "\n\nOut of memory mymalloc(%u, %u), curr %.3g bytes\n\n",
		  n, bytes1, b);
		exit(1);
		}

	void *userp = (void *) (((byte *) p) + sizeof(mem_data));
	s_Lock.lock();
	++g_mymalloc_netallocount;
	g_myalloc_netbytes += requested_bytes;
	if (s_trace)
		writes(s_trace_fd, "%s:%d\t%p\t%u\n", fn, linenr, userp, requested_bytes);

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
	bool do_write_state = (s_dump_secs > 0 && secs >= s_dump_secs && !s_write_state_done);
	s_write_state_done = true;
	s_Lock.unlock();

	if (do_summary)
		mymalloc_print_summary(s_summary_fd, "");

	if (do_write_state)
		mymalloc_write_state("myalloc.state");

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
		writes(2, "\nmyfree MAGIC ERROR\n");
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
		writes(s_trace_fd, "%u\tfree\t%p\n", secs, userp);
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

void mymalloc_write_state(const char *fn)
	{
	if (fn[0] == 0)
		return;
	int fd = open(fn, O_WRONLY | O_CREAT, S_IREAD | S_IWRITE);
	if (fd < 0)
		{
		writes(2, "\nmymalloc_write_state(%s) open error %d\n", fn, errno);
		return;
		}
	writes(2, "\nmymalloc_write_state(%s)\n", fn);
	s_Lock.lock();
	for (mem_data *Block = s_BlockList; Block != 0; Block = Block->fwd)
		writes(fd, "%p\t%u\t%u\t%u\t%s:%d\n",
				Block,
				Block->idx, 
				Block->bytes,
				Block->secs,
				GetBaseName(Block->fn),
				Block->linenr);

#if defined(__GNUC__)
	{
	int proc_maps_fd = open("/proc/self/maps", O_RDONLY);
	if (proc_maps_fd < 0)
		{
		int e = errno;
		writes(2, "\nopen(/proc/self/maps) failed errno=%d\n", e);
		}
	else
		{
		int proc_maps_copy_fd = open("myalloc.maps", O_WRONLY | O_CREAT, S_IREAD | S_IWRITE);
		if (proc_maps_copy_fd < 0)
			{
			int e = errno;
			writes(2, "\nopen(/proc/self/maps) failed errno=%d\n", e);
			}
		else
			{
			for (;;)
				{
				byte data;
				int c = read(proc_maps_fd, &data, 1);
				if (c == -1)
					break;
				write(proc_maps_copy_fd, &data, 1);
				}
			close(proc_maps_copy_fd);
			}
		close(proc_maps_fd);
		}
	}
#endif // __GNUC__
	close(fd);
	s_Lock.unlock();
	}

void mymalloc_print_summary(int fd, const char *Msg)
	{
	if (fd < 0)
		return;
	s_Lock.lock();

	writes(fd, "\n");
	if (strcmp(Msg, "") == 0)
		writes(fd, "mymalloc_print_summary(%u)\n", uint(time(0)-s_t0));
	else
		writes(fd, "mymalloc_print_summary(%s)\n", Msg);

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

	writes(fd, "%10.0f  mem_bytes (%.3g)\n", mem, mem);
	writes(fd, "%10.0f  net_bytes (%.3g)\n", net, net);
	writes(fd, "%10.0f  sum_block_bytes (%.3g)\n", sum, sum);
	writes(fd, "%10.0f  net_new_bytes (%.3g)\n", newbytes, newbytes);
	writes(fd, "%10.0f  overhead (%.3g)\n", overhead, overhead);
	writes(fd, "%10u  net_allocs\n", g_mymalloc_netallocount.load());
	writes(fd, "%10u  net_news\n", g_mymalloc_netnewcount.load());
	writes(fd, "%10u  nrblocks\n", nrblocks);
	s_Lock.unlock();
	}

void mymalloc_print_state(int fd, const char *Msg)
	{
	if (fd < 0)
		return;

	writes(fd, "\n");
	writes(fd, "mymalloc_print_state(%s)\n", Msg);
	writes(fd, "    idx     bytes   secs  src\n");
	//          1234567  12345678  12345
	uint64 sumbytes = 0;
	s_Lock.lock();
	for (mem_data *Block = s_BlockList; Block != 0; Block = Block->fwd)
		writes(fd, "%7u  %8u  %5u  %s:%d\n",
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
	if (g_fLog != 0)
		{
		int fd = fileno(g_fLog);
		mymalloc_print_state(fd, "");
		}
	mymalloc_write_state("test.state");
	}

#endif