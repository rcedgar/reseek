#ifndef myutils_h
#define myutils_h

#define RCE_MALLOC		0
#define TRACK_ALLOC		0
#define USE_OMP			0

#if defined(__x86_64__) || defined(_M_X64)
#define	BITS			64
#else
#define	BITS			32
#endif

#include <stdio.h>
#include <sys/types.h>
#include <string>
#include <string.h>
#include <memory.h>
#include <vector>
#include <math.h>
#include <stdarg.h>
#include <cstdlib>
#include <climits>
#include <float.h>
#include <algorithm>
#include <mutex>
#include <filesystem>
#include <thread>
#include <set>

#ifndef _MSC_VER
#include <inttypes.h>
#endif

#ifndef _MSC_VER
#define _stricmp	strcasecmp
#define stricmp	strcasecmp
#endif

using namespace std;

#ifdef _MSC_VER
#include <crtdbg.h>
#pragma warning(disable: 4996)	// deprecated functions
#define _CRT_SECURE_NO_DEPRECATE	1
#endif

#if defined(_DEBUG) && !defined(DEBUG)
#define DEBUG	1
#endif

#if defined(DEBUG) && !defined(_DEBUG)
#define _DEBUG	1
#endif

#ifndef NDEBUG
#define	DEBUG	1
#define	_DEBUG	1
#endif

#define byte rce__byte
typedef unsigned char byte;
typedef unsigned short uint16;
typedef unsigned uint32;
typedef unsigned uint;

// typedefs for int64 and uint64
#if		defined(_MSC_VER)
typedef __int64 int64;
typedef unsigned __int64 uint64;
#elif defined(__GNUC__)
typedef long long int64;
typedef unsigned long long uint64;
#else	
#error	"int64 typedefs"
#endif

#if	BITS==32
typedef uint32 uintb;
#else
typedef uint64 uintb;
#endif

#ifndef UINT32_MAX
const uint32 UINT32_MAX = (~(uint32(1)));
#endif

#ifndef UINT64_MAX
const uint64 UINT64_MAX = (~(uint64(1)));
#endif

#ifndef SIZE_T_MAX
const size_t SIZE_T_MAX = (~(size_t(1)));
#endif

void myassertfail(const char *Exp, const char *File, unsigned Line);
#undef  assert
#ifdef  NDEBUG
#define assert(exp)     ((void)0)
#define myassert(exp)     ((void)0)
#else
#define assert(exp) (void)( (exp) || (myassertfail(#exp, __FILE__, __LINE__), 0) )
#define myassert(exp) (void)( (exp) || (myassertfail(#exp, __FILE__, __LINE__), 0) )
#endif
#define asserta(exp) (void)( (exp) || (myassertfail(#exp, __FILE__, __LINE__), 0) )

#define ureturn(x)	return (x)

#define NotUsed(v)	((void *) &v)

// pom=plus or minus, tof=true or false, yon=yes or no
static inline char pom(bool Plus)	{ return Plus ? '+' : '-'; }
static inline char tof(bool x)		{ return x ? 'T' : 'F';	}
static inline char yon(bool x)		{ return x ? 'Y' : 'N';	}
static inline const char *YesOrNo(bool x)	{ return x ? "Yes" : "No"; }
static inline const char *plurals(unsigned n) { return n == 1 ? "" : "s"; }

extern const char *g_GitVer;
extern vector<string> g_Argv;
extern string g_Arg1;
extern mutex g_DieLock;

const char *GetPlatform();
unsigned GetElapsedSecs();
void mysleep(unsigned ms);
void mylistdir(const string &DirName, vector<string> &FileNames,
  vector<bool> &IsSubDirs);

#if	RCE_MALLOC

void *rce_malloc(size_t bytes, const char *FileName, int Line);
void rce_free(void *p, const char *FileName, int LineNr);
void rce_chkmem();

void rce_assertvalidptr_(void *p, const char *FileName, int LineNr);
#define rce_assertvalidptr(p)	rce_assertvalidptr_(p, __FILE__, __LINE__)

void rce_dumpptr_(void *p, const char *FileName, int LineNr);
#define rce_dumpptr(p)	rce_dumpptr_(p, __FILE__, __LINE__)

#define mymalloc(n)		rce_malloc((n), __FILE__, __LINE__)
#define myfree(p)		rce_free(p, __FILE__, __LINE__)
#define myalloc(t, n)	(t *) rce_malloc((n)*sizeof(t), __FILE__, __LINE__)

#else // RCE_MALLOC

void Version(FILE *f);
void *mymalloc(size_t bytes);
void myfree(void *p);
#define rce_chkmem()	/* empty */
#define myalloc(t, n)	(t *) mymalloc((n)*sizeof(t))

#endif // RCE_MALLOC

#if TRACK_ALLOC

#undef myalloc
#undef myfree

void *myalloc_track(unsigned Bytes, unsigned AllocId, const char *FileName, int LineNr);
void myfree_track(void *p, const char *FileName, int LineNr);
void myalloc_trace(bool On);

#define myalloc(t, n)	(t *) myalloc_track(sizeof(t)*(n), AllocId, __FILE__, __LINE__)
#define myfree(p)		myfree_track(p, __FILE__, __LINE__)

#endif // TRACK_ALLOC

#define SIZE(c)	unsigned((c).size())
#define RoundUp(Bytes, BlockSize)	((Bytes) + ((BlockSize) - (Bytes)%(BlockSize)))

bool myisatty(int fd);

#ifdef _MSC_VER
#define off_t	__int64
#endif

// Stdio functions without "nr of bytes" arg.
FILE *OpenStdioFile(const string &FileName);
FILE *CreateStdioFile(const string &FileName);
void CloseStdioFile(FILE *f);
bool ReadLineStdioFile(FILE *f, string &Line);
void FlushStdioFile(FILE *f);
bool StdioFileExists(const string &FileName);
void LogStdioFileState(FILE *f);
void RenameStdioFile(const string &FileNameFrom, const string &FileNameTo);
void MoveStdioFile(const string &FileName1, const string &FileName2);
void DeleteStdioFile(const string &FileName);
void WriteStdioFileStr(FILE *f, const char *s);
void Pr(FILE *f, const char *Format, ...);

// Stdio functions with size args:
byte *ReadAllStdioFile32(FILE *f, uint32 &FileSize);
byte *ReadAllStdioFile64(FILE *f, uint64 &FileSize);

byte *ReadAllStdioFile(FILE *f, uint32 &FileSize);
byte *ReadAllStdioFile64(FILE *f, uint64 &FileSize);

byte *ReadAllStdioFile32(const string &FileName, uint32 &FileSize);
byte *ReadAllStdioFile64(const string &FileName, uint64 &FileSize);

bool ReadLineStdioFile(FILE *f, char *Line, uint32 Bytes);
bool ReadLineStdioFile64(FILE *f, char *Line, uint64 Bytes);

void SetStdioFilePos(FILE *f, uint32 Pos);
void SetStdioFilePos64(FILE *f, uint64 Pos);

uint32 GetStdioFilePos32(FILE *f);
uint64 GetStdioFilePos64(FILE *f);

uint32 GetStdioFileSize32(FILE *f);
uint64 GetStdioFileSize64(FILE *f);

#if	BITS==32
#define uintB	uint32
#define GetStdioFilePosB	GetStdioFilePos32
#define GetStdioFileSizeB	GetStdioFileSize32
#define SetStdioFilePosB	SetStdioFilePos
#else
#define uintB	uint64
#define GetStdioFilePosB	GetStdioFilePos64
#define GetStdioFileSizeB	GetStdioFileSize64
#define SetStdioFilePosB	SetStdioFilePos64
#endif

uint32 ReadStdioFile_NoFail(FILE *f, void *Buffer, uint32 Bytes);
void ReadStdioFile(FILE *f, uint32 Pos, void *Buffer, uint32 Bytes);
void ReadStdioFile64(FILE *f, uint64 Pos, void *Buffer, uint64 Bytes);

void ReadStdioFile(FILE *f, void *Buffer, uint32 Bytes);
void ReadStdioFile64(FILE *f, void *Buffer, uint64 Bytes);

void WriteStdioFile(FILE *f, uint32 Pos, const void *Buffer, uint32 Bytes);
void WriteStdioFile64(FILE *f, uint64 Pos, const void *Buffer, uint64 Bytes);

void WriteStdioFile(FILE *f, const void *Buffer, uint32 Bytes);
void WriteStdioFile64(FILE *f, const void *Buffer, uint64 Bytes);

void Ps(string &Str, const char *Format, ...);
void Psa(string &Str, const char *Format, ...);
void Psasc(string &Str, const char *Format, ...);
void Pf(FILE *f, const char *Format, ...);

#define MAGIC(a, b, c, d)	uint32(uint32(a)<<24 | uint32(b)<<16 | uint32(c)<<8 | (d))

void myvstrprintf(string &Str, const char *szFormat, va_list ArgList);
void myvstrprintf(string &Str, const char *szFormat, ...);

void SetLogFileName(const string &FileName);
void Log(const char *szFormat, ...);

void Die_(const char *szFormat, ...);
void Warning_(const char *szFormat, ...);

typedef void (*PTR_PRINTFLIKE_FN)(const char *Format, ...);

static inline PTR_PRINTFLIKE_FN DiePtr(const char *FileName, unsigned LineNr)
	{
	fprintf(stderr, "\n\n%s(%u): ", FileName, LineNr);
	Log("\n\n%s(%u): ", FileName, LineNr);
	return Die_;
	}

static inline PTR_PRINTFLIKE_FN WarningPtr(const char *FileName, unsigned LineNr)
	{
	fprintf(stderr, "\n\n%s(%u): ", FileName, LineNr);
	Log("\n\n%s(%u): ", FileName, LineNr);
	return Warning_;
	}

#define Die		(g_DieLock.lock(), *DiePtr(__FILE__, __LINE__))
#define Warning	(*WarningPtr(__FILE__, __LINE__))

typedef const char *(*FN_PROGRESS_CALLBACK)();
void SetPCB(FN_PROGRESS_CALLBACK PCB);
void ProgressCallback(unsigned i, unsigned N);

bool ProgressPrefix(bool On);
string &GetProgressPrefixStr(string &s);
void ProgressStep(size_t i, size_t N, const char *Format, ...);
void Progress(const char *szFormat, ...);
void Progress(const string &Str);
void ProgressLog(const char *szFormat, ...);
void ProgressLogPrefix(const char *Format, ...);

void ProgressFileInit(FILE *f, const char *Format = 0, ...);
void ProgressFileStep(const char *Format = 0, ...);
void ProgressFileDone(const char *Format = 0, ...);

void LogElapsedTimeAndRAM();
void LogProgramInfoAndCmdLine();
void Help();

inline unsigned ustrlen(const char *s) { return (unsigned) strlen(s); }
inline unsigned ustrlen(const string &s) { return SIZE(s); }
char *mystrsave(const char *s);
unsigned myipow(unsigned x, unsigned y);
static inline unsigned atou(const char *s) { return unsigned(atoi(s)); }
static inline unsigned atou(const string &s) { return unsigned(atoi(s.c_str())); }

unsigned GetThreadIndex();

double GetMemUseBytes();
double GetPeakMemUseBytes();
double GetPhysMemBytes();
double GetUsableMemBytes();

// Are two floats equal to within epsilon?
const double epsilon = 0.01;
inline bool feq(double x, double y, double epsilon)
	{
	if (fabs(x) > 10000)
		epsilon = fabs(x)/10000;
	if (fabs(x - y) > epsilon)
		return false;
	return true;
	}

inline bool feq(double x, double y)
	{
	if (x < -1e6 && y < -1e6)
		return true;
	double e = epsilon;
	if (fabs(x) > 10000)
		e = fabs(x)/10000;
	if (fabs(x - y) > e)
		return false;
	return true;
	}

#define asserteq(x, y)	assert(feq(x, y))
#define assertaeq(x, y)	asserta(feq(x, y))

#define	zero_array(a, n)	memset(a, 0, (n)*sizeof(a[0]))
#define	memset_zero(a, n)	zero_array((a), (n))

void InitRand();
unsigned randu32();
void SplitWhite(const string &Str, vector<string> &Fields);
void Split(const string &Str, vector<string> &Fields, char Sep = '\t');
bool EndsWith(const string &s, const string &t);
bool StartsWith(const string &s, const string &t);
bool StartsWith(const char *s, const char *t);
static inline double GetRatio(double x, double y) { if (y == 0) { asserta(x == 0); return 0; } return x/y; }
static inline double GetPct(double x, double y) { return 100.0*GetRatio(x, y); }
const char *GetPctStr(double x, double y, string &s);
double GetMemUseBytes();
void PrintCmdLine(FILE *f);
void PrintProgramInfo(FILE *f);
void PrintCopyright(FILE *f);
extern string g_ShortCmdLine;

const char *MemBytesToStr(double Bytes);
static inline const char *MemBytesToStr(uint64 Bytes) { return MemBytesToStr((double) Bytes); }
static inline const char *MemBytesToStr(uint32 Bytes) { return MemBytesToStr((double) Bytes); }
unsigned StrToUint(const string &s);
int StrToInt(const string &s);
double StrToMemBytes(const string &s);
double StrToFloat(const string &s);
double StrToFloat(const char *s);
const char *GetElapsedTimeStr(string &s);
const char *GetMaxRAMStr(string &s);

const char *GetBaseName(const char *PathName);
void GetBaseName(const string &PathName, string &Base);
void GetStemName(const string &PathName, string &Stem);
const char *GetExtFromPathName(const char *PathName);
void GetExtFromPathName(const string &PathName, string &Ext);

const char *IntToStr(unsigned i);
const char *IntToStr2(unsigned i);
const char *Int64ToStr(uint64 i);
const char *FloatToStr(double d);
const char *IntFloatToStr(double d);
const char *SecsToStr(double Secs);

void LogInt(unsigned i, unsigned w = UINT_MAX);
void Logu(unsigned u, unsigned w, unsigned prefixspaces = 2);
void Logf(float x, unsigned w, unsigned prefixspaces = 2);
const char *SecsToHHMMSS(int Secs);
unsigned GetCPUCoreCount();

void MyCmdLine(int argc, char **argv);
//void CmdLineErr(const char *Format, ...);
void GetCmdLine(string &s);

#define oc(optname)	(opt_##optname)
#define os(optname)	(string(opt_##optname))

#define FLAG_OPT(Name)						extern bool opt_##Name; extern bool optset_##Name;
#define UNS_OPT(Name, Default, Min, Max)	extern unsigned opt_##Name; extern bool optset_##Name;
#define FLT_OPT(Name, Default, Min, Max)	extern double opt_##Name; extern bool optset_##Name;
#define STR_OPT(Name)						extern const char *opt_##Name; extern bool optset_##Name;
#include "myopts.h"

extern FILE *g_fLog;

void LogAllocs();

unsigned GetRequestedThreadCount();

void SeqToFasta(FILE *f, const char *Label, const char *Seq, unsigned L,
  uint ROWLEN = 80);
void SeqToFasta(FILE *f, const string &Label, const string &Seq,
  uint ROWLEN = 80);
void RevCompSeq(string &Seq);
void StripGaps(string &Seq);
void ToLower(string &Str);
void StripWhiteSpace(string &Str);
char GetOneFromThree(const string &AAA);
void ReadLinesFromFile(const string &FileName, vector<string> &Lines);
void Shuffle(vector<unsigned> &v);
void ToUpper(string &s);
void ToLower(string &s);
bool IsDirectory(const string &PathName);
bool IsRegularFile(const string &PathName);

#ifdef _MSC_VER
#define brk(x)       if (x) __debugbreak()
#else
#define brk(x)		(0)
#endif

#endif	// myutils_h
