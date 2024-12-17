#include "myutils.h"

#if _MSC_VER && DEBUG

#if MYMALLOC_TRACE
static mutex s_Lock;
static bool s_trace;
#endif

void MemTrace(bool On)
	{
#if MYMALLOC_TRACE
	s_trace = On;
#endif
	}

struct HEAP_BLOCK
	{
// Pointer to the block allocated just before this one:
	HEAP_BLOCK *_block_header_next;

// Pointer to the block allocated just after this one:
	HEAP_BLOCK *_block_header_prev;
	char const *_file_name;
	int                 _line_number;

	int                 _block_use;      // Type of block
	size_t              _data_size;      // Size of user block

	long                _request_number; // Allocation number

// Buffer just before (lower than) the user's memory:
	unsigned char       _gap[4];

// Followed by:
// unsigned char    _data[_data_size];
// unsigned char    _another_gap[no_mans_land_size];
	};

void *operator new(size_t n)
	{
	void *p = _malloc_dbg(n, _CLIENT_BLOCK, __FILE__, __LINE__);
	return p;
	}

void operator delete(void *p)
	{
	_free_dbg(p, _CLIENT_BLOCK);
	}

static FILE *f;
static time_t t0;

void *mymalloc(unsigned n, unsigned bytes, const char *fn, int linenr)
	{
	if (bytes == 0)
		return 0;
	size_t B = size_t(n)*size_t(bytes);
	void *p = _malloc_dbg(B, _NORMAL_BLOCK, fn, linenr);
	if (0 == p)
		{
		double b = GetMemUseBytes();
		fprintf(stderr, "\n\n%s:%d Out of memory mymalloc(%u, %u), curr %.3g bytes\n\n",
		  fn, linenr, n, bytes, b);
		exit(1);
		}

#if MYMALLOC_TRACE
	s_Lock.lock();
	if (f == 0)
		{
		f = fopen("myalloc.trace", "w");
		t0 = time(0);
		setbuf(f, 0);
		}
	uint secs = uint(time(0) - t0);
	fprintf(f, "%u\talloc\t%s\t%d\t%d\t%u\t%p\n", secs, fn, linenr, n, bytes, p);
	s_Lock.unlock();
#endif
	return p;
	}

void myfree(void *p)
	{
	if (p == 0)
		return;
#if MYMALLOC_TRACE
	s_Lock.lock();
	uint secs = uint(time(0) - t0);
	fprintf(f, "%u\tfree\t%p\n", secs, p);
	s_Lock.unlock();
#endif
	_free_dbg(p, _NORMAL_BLOCK);
	}

//#define _FREE_BLOCK      0
//#define _NORMAL_BLOCK    1
//#define _CRT_BLOCK       2
//#define _IGNORE_BLOCK    3
//#define _CLIENT_BLOCK    4
//#define _MAX_BLOCKS      5

void LogHeapSummary(const char *s)
	{
	static const char *TypeNames[6] =
		{
		"Free",
		"Normal",
		"CRT",
		"Ignore",
		"Client",
		"ERROR"
		};
	Log("\n");
	Log("LogHeapSummary(%s)\n", s);
	_CrtMemState MS;
	_CrtMemCheckpoint(&MS);
	uint n = 0;
	uint Counts[6];
	uint64_t Totals[6];
	uint64_t Totals2[6];
	zero_array(Counts, 6);
	zero_array(Totals, 6);
	zero_array(Totals2, 6);
	int n1 = 0;
	for (HEAP_BLOCK *b = (HEAP_BLOCK *) MS.pBlockHeader;
		 b != 0; b = b->_block_header_next)
		{
		int Type = b->_block_use & 0xffff;
		int SubType = (b->_block_use >> 16);
		if (SubType == 1)
			++n1;
		if (Type >= 0 && Type < 5)
			{
			++Counts[Type];
			Totals[Type] += b->_data_size;
			Totals2[Type] += b->_data_size + sizeof(HEAP_BLOCK);
			}
		else
			{
			++Counts[5];
			Totals[5] += b->_data_size;
			}
		}

	double Total = 0;
	double Total2 = 0;
	for (int i = 0; i < 6; ++i)
		{
		double t = double(Totals[i]);
		double t2 = double(Totals2[i]);
		Total += t;
		Total2 += t2;
		Log("%8.8s", TypeNames[i]);
		Log("  %7u", Counts[i]);
		Log("  %10.0f", double(t));
		Log("  %10.0f", double(t2));
		Log("  %s", MemBytesToStr(double(t)));
		Log("\n");
		}
	double MemUse = GetMemUseBytes();
	Log("MemUse = %.0f (%s)\n", MemUse, MemBytesToStr(MemUse));
	Log("Total = %.0f (%s)\n", Total, MemBytesToStr(double(Total)));
	Log("Total2 = %.0f (%s)\n", Total2, MemBytesToStr(double(Total2)));
	}

void LogHeapBlocks(int aType, const char *s)
	{
	_CrtCheckMemory();
	_CrtMemState MS;
	_CrtMemCheckpoint(&MS);
	//Log("LogHeapBlocks(%d, %s)\n", aType, s);

	uint Count = 0;
	double Total = 0;
	for (HEAP_BLOCK *b = (HEAP_BLOCK *) MS.pBlockHeader;
		 b != 0; b = b->_block_header_next)
		{
		int Type = b->_block_use & 0xffff;
		if (b->_file_name == 0)
			Log("@HB_%s\tNULL:0\t%.3g\n", s, double(b->_data_size));
		else
			Log("@HB_%s\t%s:%d\t%.3g\n", s, b->_file_name, b->_line_number, double(b->_data_size));
		if (Type == aType)
			{
			++Count;
			++Total;
			}
		}

	Log("Blocks %u, total = %.0f (%s)\n",
		Count, Total, MemBytesToStr(double(Total)));
	}

#endif // DEBUG
