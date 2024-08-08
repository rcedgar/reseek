#include <time.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <signal.h>
#include <float.h>
#include <mutex>

#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include <crtdbg.h>
#include <process.h>
#include <windows.h>
#include <psapi.h>
#include <io.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <dirent.h>
#endif

#include "myutils.h"

bool IsDirectory(const string& PathName)
{
	std::error_code ec;
	bool IsDir = std::filesystem::is_directory(PathName);
	return IsDir;
}

bool IsRegularFile(const string& PathName)
{
	std::error_code ec;
	bool IsFile = std::filesystem::is_regular_file(PathName);
	return IsFile;
}

const int SIZE_16 = 16;
const int SIZE_32 = 32;
const int SIZE_64 = 64;

mutex g_DieLock;

static void CompilerInfo();

static map<FILE *, string> FileToFileName;

const unsigned MY_IO_BUFSIZ = 32000;
const unsigned MAX_FORMATTED_STRING_LENGTH = 64000;


static char *g_IOBuffers[256];
static time_t g_StartTime = time(0);
extern vector<string> g_Argv;
static double g_PeakMemUseBytes;

static char **g_ThreadStrs;
static unsigned g_ThreadStrCount;

const char *GetPlatform()
	{
#if	BITS==32
	asserta(sizeof(void *) == 4);
#ifdef _MSC_VER
	return "win32";
#elif defined(__APPLE__)
	return "i86osx32";
#elif defined(__GNUC__)
	return "i86linux32";
#else
#error "Unknown compiler"
#endif
#elif BITS==64
	asserta(sizeof(void *) == 8);
#ifdef _MSC_VER
	return "win64";
#elif defined(__APPLE__)
	return "i86osx64";
#elif defined(__GNUC__)
	return "i86linux64";
#else
#error "Unknown compiler"
#endif
#else
#error "Bad BITS"
#endif
	}

const char *GetBaseName(const char *PathName)
	{
	const char *q = 0;
	for (const char *p = PathName; *p; ++p)
		{
		char c = *p;
		if (c == '/' || c == '\\')
			q = p + 1;
		}
	if (q != 0)
		return q;
	return PathName;
	}

void GetBaseName(const string &PathName, string &Base)
	{
	Base = string(GetBaseName(PathName.c_str()));
	}

void GetStemName(const string &PathName, string &Stem)
	{
	string Base;
	GetBaseName(PathName, Base);
	vector<string> Fields;
	Split(Base, Fields, '.');
	Stem = Fields[0];
	}

void GetExtFromPathName(const string &PathName, string &Ext)
	{
	string Base;
	GetBaseName(PathName, Base);
	vector<string> Fields;
	Split(Base, Fields, '.');
	uint n = SIZE(Fields);
	if (n == 1)
		{
		Ext = "";
		return;
		}
	const string &LastField = Fields[n-1];

// Special case for .ext.gz, this is considered to be 
//   the filename extension
	if (LastField == "gz" && n > 2)
		{
		Ext = Fields[n-2] + ".gz";
		return;
		}
	Ext = LastField;
	}

static void AllocBuffer(FILE *f)
	{
#if	DEBUG
	setbuf(f, 0);
#else
	int fd = fileno(f);
	if (fd < 0 || fd >= 256)
		return;
	if (g_IOBuffers[fd] == 0)
		g_IOBuffers[fd] = myalloc(char, MY_IO_BUFSIZ);
	setvbuf(f, g_IOBuffers[fd], _IOFBF, MY_IO_BUFSIZ);
#endif
	}

static void FreeBuffer(FILE *f)
	{
#if	0
	int fd = fileno(f);
	if (fd < 0 || fd >= 256)
		return;
	if (g_IOBuffers[fd] == 0)
		return;
	myfree(g_IOBuffers[fd]);
	g_IOBuffers[fd] = 0;
#endif
	}

unsigned GetElapsedSecs()
	{
	return (unsigned) (time(0) - g_StartTime);
	}

bool StdioFileExists(const string &FileName)
	{
	struct stat SD;
	int i = stat(FileName.c_str(), &SD);
	return i == 0;
	return true;
	}

static mutex g_Lock;
void myassertfail(const char *Exp, const char *File, unsigned Line)
	{
	g_Lock.lock();
	Die("%s(%u) assert failed: %s", File, Line, Exp);
	g_Lock.unlock();
	}

bool myisatty(int fd)
	{
	return isatty(fd) != 0;
	}

#ifdef _MSC_VER
#include <io.h>
int fseeko(FILE *stream, off_t offset, int whence)
	{
	off_t FilePos = _fseeki64(stream, offset, whence);
	return (FilePos == -1L) ? -1 : 0;
	}
#define ftello(fm) (off_t) _ftelli64(fm)
#endif

void LogStdioFileState(FILE *f)
	{
	unsigned long tellpos = (unsigned long) ftello(f);
	long fseek_pos = fseek(f, 0, SEEK_CUR);
	int fd = fileno(f);
	Log("FILE *     %p\n", f);
	Log("fileno     %d\n", fd);
	Log("feof       %d\n", feof(f));
	Log("ferror     %d\n", ferror(f));
	Log("ftell      %ld\n", tellpos);
	Log("fseek      %ld\n", fseek_pos);
#if	!defined(_GNU_SOURCE) && !defined(__APPLE_CC__)
	fpos_t fpos;
	int fgetpos_retval = fgetpos(f, &fpos);
	Log("fpos       %ld (retval %d)\n", (long) fpos, fgetpos_retval);
//	Log("eof        %d\n", _eof(fd));
#endif
#ifdef _MSC_VER
	__int64 pos64 = _ftelli64(f);
	Log("_ftelli64  %lld\n", pos64);
#endif
	if (FileToFileName.find(f) == FileToFileName.end())
		Log("Not found in FileToFileName\n");
	else
		Log("Name       %s\n", FileToFileName[f].c_str());
	}

FILE *OpenStdioFile(const string &FileName)
	{
	if (FileName == "")
		Die("Missing input file name");

	const char *Mode = "rb";
	FILE *f = fopen(FileName.c_str(), Mode);
	if (f == 0)
		{
		if (errno == EFBIG)
			{
			if (sizeof(off_t) == 4)
				Die("File too big for 32-bit version (sizeof(off_t)=%d): %s", sizeof(off_t), FileName.c_str());
			else
				Die("Cannot open '%s', file too big (off_t=%u bits)",
				  FileName.c_str(), sizeof(off_t)*8);
			}
		int Err = errno;
		const char *StrErr = strerror(errno);
		Die("Cannot open %s, errno=%d %s",
		  FileName.c_str(), Err, StrErr);
		}
	AllocBuffer(f);
	FileToFileName[f] = FileName;
	return f;
	}

FILE *CreateStdioFile(const string &FileName)
	{
	if (FileName == "")
		// Die("Missing output file name");
		return 0;
	FILE *f = fopen(FileName.c_str(), "wb+");
	if (0 == f)
		Die("Cannot create %s, errno=%d %s",
		  FileName.c_str(), errno, strerror(errno));
	AllocBuffer(f);
	FileToFileName[f] = FileName;
	return f;
	}

void SetStdioFilePos(FILE *f, uint32 Pos)
	{
	if (0 == f)
		Die("SetStdioFilePos failed, f=NULL");
	int Ok = fseeko(f, Pos, SEEK_SET);
	off_t NewPos = ftello(f);
	if (Ok != 0 || Pos != NewPos)
		{
		LogStdioFileState(f);
		Die("SetStdioFilePos(%d) failed, Ok=%d NewPos=%d",
		  (int) Pos, Ok, (int) NewPos);
		}
	}

void SetStdioFilePos64(FILE *f, uint64 Pos)
	{
	if (0 == f)
		Die("SetStdioFilePos failed, f=NULL");
	int Ok = fseeko(f, Pos, SEEK_SET);
	off_t NewPos = ftello(f);
	if (Ok != 0 || Pos != NewPos)
		{
		LogStdioFileState(f);
		Die("SetStdioFilePos64(%ul) failed, Ok=%d NewPos=%ul",
		  (unsigned long) Pos, Ok, (unsigned long)  NewPos);
		}
	}

uint32 ReadStdioFile_NoFail(FILE *f, void *Buffer, uint32 Bytes)
	{
	asserta(f != 0);
	off_t PosBefore = ftello(f);
	size_t ElementsRead = fread(Buffer, Bytes, 1, f);
	off_t PosAfter = ftello(f);
	if (ElementsRead == 1)
		return Bytes;
	uint32 BytesRead = uint32(PosAfter - PosBefore);
	return BytesRead;
	}

void ReadStdioFile(FILE *f, uint32 Pos, void *Buffer, uint32 Bytes)
	{
	asserta(f != 0);
	SetStdioFilePos(f, Pos);
	uint32 BytesRead = (uint32) fread(Buffer, 1, Bytes, f);
	if (BytesRead != Bytes)
		{
		LogStdioFileState(f);
		Die("ReadStdioFile failed, attempted %lu bytes, read %lu bytes, errno=%d",
		  (unsigned long) Bytes, (unsigned long) BytesRead, errno);
		}
	}

void ReadStdioFile64(FILE *f, uint64 Pos, void *Buffer, uint64 Bytes)
	{
	asserta(f != 0);
	uint32 Bytes32 = (uint32) Bytes;
	asserta(Bytes32 == Bytes);
	SetStdioFilePos64(f, Pos);
	uint64 BytesRead = (uint64) fread(Buffer, 1, Bytes32, f);
	if (BytesRead != Bytes)
		{
		LogStdioFileState(f);
		Die("ReadStdioFile64 failed, attempted %lu bytes, read %lu bytes, errno=%d",
		  (unsigned long) Bytes, (unsigned long) BytesRead, errno);
		}
	}

void ReadStdioFile(FILE *f, void *Buffer, uint32 Bytes)
	{
	asserta(f != 0);
	uint32 BytesRead = (uint32) fread(Buffer, 1, Bytes, f);
	if (BytesRead != Bytes)
		{
		LogStdioFileState(f);
		Die("ReadStdioFile failed, attempted %u bytes, read %u bytes, errno=%d",
		  Bytes, BytesRead, errno);
		}
	}

void ReadStdioFile64(FILE *f, void *Buffer, uint64 Bytes)
	{
	asserta(f != 0);
	uint32 Bytes32 = (uint32) Bytes;
	asserta(Bytes32 == Bytes);
	uint64 BytesRead = (uint32) fread(Buffer, 1, Bytes32, f);
	if (BytesRead != Bytes)
		{
		LogStdioFileState(f);
		Die("ReadStdioFile failed, attempted %u bytes, read %u bytes, errno=%d",
		  Bytes, BytesRead, errno);
		}
	}

byte *ReadAllStdioFile(FILE *f, uint32 &FileSize)
	{
	uint64 Pos = GetStdioFilePos64(f);
	uint64 FileSize64 = GetStdioFileSize64(f);
#if	BITS == 32
	if (FileSize > UINT_MAX)
		Die("ReadAllStdioFile (32-bit): file too big");
#endif
	FileSize = uint32(FileSize64);
	SetStdioFilePos(f, 0);
	byte *Buffer = myalloc(byte, FileSize);
	ReadStdioFile(f, Buffer, FileSize);
	SetStdioFilePos64(f, Pos);
	return Buffer;
	}

byte *ReadAllStdioFile64(const string &FileName, uint64 &FileSize)
	{
	FILE *f = OpenStdioFile(FileName);
	FileSize = GetStdioFileSize64(f);
#if	BITS==32
	if (FileSize > UINT32_MAX)
		Die("File too big, requires 64-bit version: %s", FileName.c_str());
#endif
	byte *Buffer = ReadAllStdioFile64(f, FileSize);
	CloseStdioFile(f);
	return Buffer;
	}

byte *ReadAllStdioFile64(FILE *f, uint64 &FileSize)
	{
	uint64 SavedPos = GetStdioFilePos64(f);
	FileSize = GetStdioFileSize64(f);

#if BITS==32
	if (FileSize > UINT32_MAX)
		Die("File too big, requires 64-bit version");
	byte *Buffer = myalloc(byte, (uint32) FileSize);
#else
	if (FileSize > UINT_MAX)
		Die("ReadAllStdioFile64, file too big %s", MemBytesToStr((double) FileSize));
	unsigned uFileSize = (unsigned) FileSize;
	byte *Buffer = myalloc(byte, uFileSize);
#endif

	uint64 Pos = 0;
	uint64 BytesLeft = FileSize;

	const uint64 ChunkSize = 0x40000000; // 1Gb
	for (;;)
		{
		if (BytesLeft == 0)
			break;
		uint64 BytesToRead = BytesLeft;
		if (BytesToRead > ChunkSize)
			BytesToRead = ChunkSize;
		ReadStdioFile64(f, Pos, Buffer + Pos, BytesToRead);
		BytesLeft -= BytesToRead;
		}

	SetStdioFilePos64(f, SavedPos);
	return Buffer;
	}

byte *ReadAllStdioFile32(const std::string &FileName, uint32 &FileSize)
	{
#if	WIN32
	FILE *f = OpenStdioFile(FileName);
	FileSize = GetStdioFileSize32(f);
	CloseStdioFile(f);

	HANDLE h = CreateFile(FileName.c_str(), GENERIC_READ, FILE_SHARE_READ,
	  NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (h == INVALID_HANDLE_VALUE)
		Die("ReadAllStdioFile:Open(%s) failed", FileName.c_str());

	byte *Buffer = myalloc(byte, FileSize);
	DWORD BytesRead;
	ReadFile(h, Buffer, FileSize, &BytesRead, NULL);
	if (FileSize != BytesRead)
		Die("ReadAllStdioFile:Error reading %s, attempted %u got %u",
		  FileName.c_str(), FileSize, (unsigned) BytesRead);

	CloseHandle(h);
	return Buffer;
#else
	int h = open(FileName.c_str(), O_RDONLY);
	if (h < 0)
		Die("ReadAllStdioFile:Cannot open %s", FileName.c_str());
	FileSize = lseek(h, 0, SEEK_END);
#ifndef __APPLE__
	if (FileSize == (off_t) (-1))
		Die("ReadAllStdioFile:Error seeking %s", FileName.c_str());
#endif
	// byte *Buffer = myalloc<byte>(FileSize);
	size_t stBytes = (size_t) FileSize;
	if ((off_t) stBytes != FileSize)
		Die("ReadAllStdioFile: off_t overflow");
	byte *Buffer = (byte *) myalloc(byte, stBytes);
	if (Buffer == 0)
		Die("ReadAllStdioFile: failed to allocate %s", MemBytesToStr((double) stBytes));
	lseek(h, 0, SEEK_SET);
	size_t n = read(h, Buffer, stBytes);
	if (n != FileSize)
		Die("ReadAllStdioFile, Error reading %s, attempted %g got %g",
		  FileName.c_str(), (double) FileSize, (double) n);
	close(h);
	return Buffer;
#endif
#undef AllocId
#define AllocId ALLOCID_myutils
	}

void WriteStdioFile(FILE *f, uint32 Pos, const void *Buffer, uint32 Bytes)
	{
	if (0 == f)
		Die("WriteStdioFile failed, f=NULL");
	SetStdioFilePos(f, Pos);
	size_t BytesWritten = fwrite(Buffer, 1, Bytes, f);
	if (BytesWritten != Bytes)
		{
		LogStdioFileState(f);
		Die("WriteStdioFile failed, attempted %ul bytes, wrote %ul bytes, errno=%d",
		  (unsigned long) Bytes, (unsigned long) BytesWritten, errno);
		}
	}

void WriteStdioFileStr(FILE *f, const char *s)
	{
	uint32 Bytes = ustrlen(s);
	WriteStdioFile(f, s, Bytes);
	}

void WriteStdioFile(FILE *f, const void *Buffer, uint32 Bytes)
	{
	if (0 == f)
		Die("WriteStdioFile failed, f=NULL");
	size_t BytesWritten = fwrite(Buffer, 1, Bytes, f);
	if (BytesWritten != Bytes)
		{
		LogStdioFileState(f);
		Die("WriteStdioFile failed, attempted %ul bytes, wrote %ul bytes, errno=%d",
		  (unsigned long) Bytes, (unsigned long) BytesWritten, errno);
		}
	}

void WriteStdioFile64(FILE *f, const void *Buffer, uint64 Bytes)
	{
	if (0 == f)
		Die("WriteStdioFile failed, f=NULL");
	size_t BytesWritten = fwrite(Buffer, 1, Bytes, f);
	if (BytesWritten != Bytes)
		{
		LogStdioFileState(f);
		Die("WriteStdioFile64 failed, attempted %ul bytes, wrote %ul bytes, errno=%d",
		  (unsigned long) Bytes, (unsigned long) BytesWritten, errno);
		}
	}

// Return false on EOF, true if line successfully read.
bool ReadLineStdioFile(FILE *f, char *Line, uint32 Bytes)
	{
	if (feof(f))
		return false;
	if ((int) Bytes < 0)
		Die("ReadLineStdioFile: Bytes < 0");
	char *RetVal = fgets(Line, (int) Bytes, f);
	if (NULL == RetVal)
		{
		if (feof(f))
			return false;
		if (ferror(f))
			Die("ReadLineStdioFile: errno=%d", errno);
		Die("ReadLineStdioFile: fgets=0, feof=0, ferror=0");
		}

	if (RetVal != Line)
		Die("ReadLineStdioFile: fgets != Buffer");
	size_t n = strlen(Line);
	if (n < 1 || Line[n-1] != '\n')
		Die("ReadLineStdioFile: line too long or missing end-of-line");
	if (n > 0 && (Line[n-1] == '\r' || Line[n-1] == '\n'))
		Line[n-1] = 0;
	if (n > 1 && (Line[n-2] == '\r' || Line[n-2] == '\n'))
		Line[n-2] = 0;
	return true;
	}

// Return false on EOF, true if line successfully read.
bool ReadLineStdioFile(FILE *f, string &Line)
	{
	Line.clear();
	for (;;)
		{
		int c = fgetc(f);
		if (c == -1)
			{
			if (feof(f))
				{
				if (!Line.empty())
					return true;
				return false;
				}
			Die("ReadLineStdioFile, errno=%d", errno);
			}
		if (c == '\r')
			continue;
		if (c == '\n')
			return true;
		Line.push_back((char) c);
		}
	}

void RenameStdioFile(const string &FileNameFrom, const string &FileNameTo)
	{
	int Ok = rename(FileNameFrom.c_str(), FileNameTo.c_str());
	if (Ok != 0)
		Die("RenameStdioFile(%s,%s) failed, errno=%d %s",
		  FileNameFrom.c_str(), FileNameTo.c_str(), errno, strerror(errno));
	}

void FlushStdioFile(FILE *f)
	{
	int Ok = fflush(f);
	if (Ok != 0)
		Die("fflush(%p)=%d,", f, Ok);
	}

void CloseStdioFile(FILE *f)
	{
	if (f == 0)
		return;
	int Ok = fclose(f);
	if (Ok != 0)
		Die("fclose(%p)=%d", f, Ok);
	FreeBuffer(f);
	}

uint32 GetStdioFilePos32(FILE *f)
	{
	off_t FilePos = ftello(f);
	if (FilePos < 0)
		Die("ftello=%d", (int) FilePos);
	if (FilePos > UINT32_MAX)
		Die("File offset too big for 32-bit version (%s)", MemBytesToStr((double) FilePos));
	return (uint32) FilePos;
	}

uint64 GetStdioFilePos64(FILE *f)
	{
	off_t FilePos = ftello(f);
	if (FilePos < 0)
		Die("ftello=%d", (int) FilePos);
	return (uint64) FilePos;
	}

uint32 GetStdioFileSize32(FILE *f)
	{
	uint32 CurrentPos = GetStdioFilePos32(f);
	int Ok = fseeko(f, 0, SEEK_END);
	if (Ok < 0)
		Die("fseek in GetFileSize");

	off_t Length = ftello(f);
	SetStdioFilePos(f, CurrentPos);

	if (Length < 0)
		Die("ftello in GetFileSize");
#if	BITS == 32
	if (Length > UINT32_MAX)
		Die("File size too big for 32-bit version (%s)", MemBytesToStr((double) Length));
#endif
	return (uint32) Length;
	}

uint64 GetStdioFileSize64(FILE *f)
	{
	uint64 CurrentPos = GetStdioFilePos64(f);
	int Ok = fseeko(f, 0, SEEK_END);
	if (Ok < 0)
		Die("fseek in GetFileSize64");

	off_t Length = ftello(f);
	SetStdioFilePos64(f, CurrentPos);

	if (Length < 0)
		Die("ftello in GetFileSize");
	return (uint64) Length;
	}

void MoveStdioFile(const string &FileName1, const string &FileName2)
	{
	if (StdioFileExists(FileName2))
		DeleteStdioFile(FileName2);
	RenameStdioFile(FileName1, FileName2);
	}

void DeleteStdioFile(const string &FileName)
	{
	int Ok = remove(FileName.c_str());
	if (Ok != 0)
		Die("remove(%s) failed, errno=%d %s", FileName.c_str(), errno, strerror(errno));
	}

double GetUsableMemBytes()
	{
	double RAM = GetPhysMemBytes();
#if	BITS==32
#ifdef	_MSC_VER
	if (RAM > 2e9)
		return 2e9;
#else
	if (RAM > 4e9)
		return 4e9;
#endif
#endif
	return RAM;
	}

void myvstrprintf(string &Str, const char *Format, va_list ArgList)
	{
	char *szStr = myalloc(char, MAX_FORMATTED_STRING_LENGTH);
	vsnprintf(szStr, MAX_FORMATTED_STRING_LENGTH-1, Format, ArgList);
	szStr[MAX_FORMATTED_STRING_LENGTH - 1] = '\0';
	Str.assign(szStr);
	myfree(szStr);
	}

void myvstrprintf(string &Str, const char *Format, ...)
	{
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);
	}

FILE *g_fLog = 0;

void SetLogFileName(const string &FileName)
	{
	if (g_fLog != 0)
		CloseStdioFile(g_fLog);
	g_fLog = 0;
	if (FileName.empty())
		return;
	g_fLog = CreateStdioFile(FileName);
	}

void Log(const char *Format, ...)
	{
	if (g_fLog == 0)
		return;
	va_list ArgList;
	va_start(ArgList, Format);
	vfprintf(g_fLog, Format, ArgList);
	va_end(ArgList);
	fflush(g_fLog);
	}

void Die_(const char *Format, ...)
	{
	string Msg;

	if (g_fLog != 0)
		setbuf(g_fLog, 0);
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Msg, Format, ArgList);
	va_end(ArgList);

	fprintf(stderr, "\n\n");
	Log("\n");
	time_t t = time(0);
	Log("%s", asctime(localtime(&t)));
	for (unsigned i = 0; i < g_Argv.size(); i++)
		{
		fprintf(stderr, (i == 0) ? "%s" : " %s", g_Argv[i].c_str());
		Log((i == 0) ? "%s" : " %s", g_Argv[i].c_str());
		}
	fprintf(stderr, "\n");
	Log("\n");

	time_t CurrentTime = time(0);
	unsigned ElapsedSeconds = unsigned(CurrentTime - g_StartTime);
	const char *sstr = SecsToStr(ElapsedSeconds);
	Log("Elapsed time: %s\n", sstr);

	const char *szStr = Msg.c_str();
	fprintf(stderr, "\n---Fatal error---\n%s\n", szStr);
	Log("\n---Fatal error---\n%s\n", szStr);

#ifdef _MSC_VER
	if (IsDebuggerPresent())
		_CrtSetDbgFlag(0);
	__debugbreak();
#endif
	g_DieLock.unlock();
	exit(1);
	}

void Warning_(const char *Format, ...)
	{	
	string Msg;

	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Msg, Format, ArgList);
	va_end(ArgList);

	const char *szStr = Msg.c_str();

	fprintf(stderr, "\nWARNING: %s\n\n", szStr);
	if (g_fLog != stdout)
		{
		Log("\nWARNING: %s\n", szStr);
		fflush(g_fLog);
		}
	}

#ifdef _MSC_VER
void mysleep(unsigned ms)
	{
	Sleep(ms);
	}
#else
void mysleep(unsigned ms)
	{
	usleep(ms);
	}
#endif

#ifdef _MSC_VER
double GetMemUseBytes()
	{
	HANDLE hProc = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS PMC;
	BOOL bOk = GetProcessMemoryInfo(hProc, &PMC, sizeof(PMC));
	if (!bOk)
		return 1000000;
	double Bytes = (double) PMC.WorkingSetSize;
	if (Bytes > g_PeakMemUseBytes)
		g_PeakMemUseBytes = Bytes;
	return Bytes;
	}

double GetPhysMemBytes()
	{
	MEMORYSTATUSEX MS;
	MS.dwLength = sizeof(MS);
	BOOL Ok = GlobalMemoryStatusEx(&MS);
	if (!Ok)
		return 0.0;
	return double(MS.ullTotalPhys);
	}

#elif	linux || __linux__ || __CYGWIN__
double GetPhysMemBytes()
	{
	int fd = open("/proc/meminfo", O_RDONLY);
	if (fd < 0)
		return 0.0;
// MemTotal:       255908 kB
	char Line[128];
	int n = read(fd, Line, sizeof(Line));
	if (n < 0)
		return 0.0;
	Line[127] = 0;
	unsigned kb;
	n = sscanf(Line, "MemTotal: %u", &kb);
	if (n != 1)
		return 0.0;
	return double(kb)*1000.0;
	}

double GetMemUseBytes()
	{
	static char statm[SIZE_64];
	static int PageSize = 1;
	if (0 == statm[0])
		{
		PageSize = sysconf(_SC_PAGESIZE);
		pid_t pid = getpid();
		snprintf(statm, SIZE_64, "/proc/%d/statm", (int) pid);
		}

	int fd = open(statm, O_RDONLY);
	if (fd < 0)
		return 0.0;
	char Buffer[64];
	int n = read(fd, Buffer, sizeof(Buffer) - 1);
	close(fd);
	fd = -1;

	if (n <= 0)
		return 0.0;

	Buffer[n] = 0;
	const char *p = strchr(Buffer, ' ');
	if (p == 0)
		return 0.0;

	double Pages = atof(p);

	double Bytes = Pages*PageSize;
	if (Bytes > g_PeakMemUseBytes)
		g_PeakMemUseBytes = Bytes;
	return Bytes;
	}

#elif defined(__MACH__)
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/socket.h>
#include <sys/gmon.h>
#include <mach/vm_param.h>
#include <netinet/in.h>
#include <netinet/icmp6.h>
#include <sys/vmmeter.h>
#include <sys/proc.h>
// #include <mach/task_info.h>
#include <mach/task.h>
#include <mach/mach_init.h>
#include <mach/vm_statistics.h>

#define DEFAULT_MEM_USE	0.0

double GetMemUseBytes()
	{
	task_t mytask = mach_task_self();
	struct task_basic_info ti;
	memset((void *) &ti, 0, sizeof(ti));
	mach_msg_type_number_t count = TASK_BASIC_INFO_COUNT;
	kern_return_t ok = task_info(mytask, TASK_BASIC_INFO, (task_info_t) &ti, &count);
	if (ok == KERN_INVALID_ARGUMENT)
		return DEFAULT_MEM_USE;

	if (ok != KERN_SUCCESS)
		return DEFAULT_MEM_USE;

	double Bytes = (double ) ti.resident_size;
	if (Bytes > g_PeakMemUseBytes)
		g_PeakMemUseBytes = Bytes;
	return Bytes;
	}

double GetPhysMemBytes()
	{
	uint64_t mempages = 0;
	size_t len = sizeof(mempages);
	int rc = sysctlbyname("hw.memsize", &mempages, &len, NULL, 0);
	if (rc < 0)
		return 0.0;
	return double(mempages);
	}
#else
double GetMemUseBytes()
	{
	return 0.0;
	}
#endif

#ifdef _MSC_VER
// WARNING includes . and ..
void mylistdir(const string &DirName, vector<string> &FileNames,
  vector<bool> &IsSubDirs)
	{
	FileNames.clear();
	IsSubDirs.clear();
	bool First = true;
	HANDLE h = INVALID_HANDLE_VALUE;
	WIN32_FIND_DATA FFD;
	for (;;)
		{
		if (First)
			{
			string s = DirName + string("/*");
			h = FindFirstFile(s.c_str(), &FFD);
			if (h == INVALID_HANDLE_VALUE)
				return;
			First = false;
			}
		else
			{
			BOOL Ok = FindNextFile(h, &FFD);
			if (!Ok)
				return;
			}
		FileNames.push_back(string(FFD.cFileName));
		IsSubDirs.push_back(bool(FFD.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY));
		}
	}
#else
// WARNING includes . and ..
void mylistdir(const string &DirName, vector<string> &FileNames,
  vector<bool> &IsSubDirs)
	{
	FileNames.clear();
	IsSubDirs.clear();
	DIR *dir = opendir(DirName.c_str());
	if (dir == 0)
		Die("Directory not found: %s", DirName.c_str());
	for (;;)
		{
		struct dirent *dp = readdir(dir);
		if (dp == 0)
			break;
		FileNames.push_back(string(dp->d_name));
		IsSubDirs.push_back(dp->d_type == DT_DIR);
		}
	closedir(dir);
	}
#endif

double GetPeakMemUseBytes()
	{
	return g_PeakMemUseBytes;
	}

const char *SecsToHHMMSS(int Secs)
	{
	int HH = Secs/3600;
	int MM = (Secs - HH*3600)/60;
	int SS = Secs%60;
	static char Str[SIZE_32];
	if (HH == 0)
		snprintf(Str, SIZE_32, "%02d:%02d", MM, SS);
	else
		snprintf(Str, SIZE_32, "%02d:%02d:%02d", HH, MM, SS);
	return Str;
	}

const char *SecsToStr(double Secs)
	{
	if (Secs >= 60.0)
		return SecsToHHMMSS((int) Secs);

	static char Str[SIZE_16];
	if (Secs < 1e-6)
		snprintf(Str, SIZE_16, "%.2gs", Secs);
	else if (Secs < 1e-3)
		snprintf(Str, SIZE_16, "%.2fms", Secs*1e3);
	else if (Secs < 1.0)
		snprintf(Str, SIZE_16, "%.3fs", Secs);
	else if (Secs < 10.0)
		snprintf(Str, SIZE_16, "%.2fs", Secs);
	else
		snprintf(Str, SIZE_16, "%.1fs", Secs);
	return Str;
	}

const char *MemBytesToStr(double Bytes)
	{
	static char Str[SIZE_32];

	if (Bytes < 1e4)
		snprintf(Str, SIZE_32, "%.1fb", Bytes);
	else if (Bytes < 1e6)
		snprintf(Str, SIZE_32, "%.1fkb", Bytes/1e3);
	else if (Bytes < 10e6)
		snprintf(Str, SIZE_32, "%.1fMb", Bytes/1e6);
	else if (Bytes < 1e9)
		snprintf(Str, SIZE_32, "%.0fMb", Bytes/1e6);
	else if (Bytes < 100e9)
		snprintf(Str, SIZE_32, "%.1fGb", Bytes/1e9);
	else
		snprintf(Str, SIZE_32, "%.0fGb", Bytes/1e9);
	return Str;
	}

bool IsUintStr(const char *s)
	{
	if (!isdigit(*s++))
		return false;
	while (*s)
		if (!isdigit(*s++))
			return false;
	return true;
	}

int StrToInt(const char *s)
	{
	int i = atoi(s);
	return i;
	}

unsigned StrToUint(const char *s)
	{
	if (!IsUintStr(s))
		Die("Invalid integer '%s'", s);
	unsigned n = 0;
	while (char c = *s++)
		{
		if (!isdigit(c))
			return n;
		n = n*10 + (c - '0');
		}
	return n;
	}

unsigned StrToUint(const string &s)
	{
	return StrToUint(s.c_str());
	}

int StrToInt(const string &s)
	{
	asserta(s.size() > 0);
	return StrToInt(s.c_str());
	}

double StrToFloat(const char *s)
	{
	char *p;
	double d = strtod(s, &p);
	if (*p != 0)
		Die("Invalid floating-point number '%s'", s);
	return d;
	}

double StrToFloat(const string &s)
	{
	return StrToFloat(s.c_str());
	}

float StrToFloatf(const string &s)
	{
	return (float) StrToFloat(s.c_str());
	}

double StrToMemBytes(const string &s)
	{
	unsigned n = SIZE(s);
	if (n == 0)
		return 0.0;

	double d = StrToFloat(s.c_str());
	char c = toupper(s[n-1]);
	if (isdigit(c))
		return d;
	else if (c == 'K')
		return 1000.0*d;
	else if (c == 'M')
		return 1e6*d;
	else if (c == 'G')
		return 1e9*d;
	else
		Die("Invalid amount of memory '%s'", s.c_str());
	return 0.0;
	}

const char *IntToStr2(unsigned i)
	{
	static char Str[SIZE_64];
	if (i < 9999)
		snprintf(Str, SIZE_64, "%u", i);
	else
		snprintf(Str, SIZE_64, "%u (%s)", i, IntToStr(i));
	return Str;
	}

const char *IntToStr(unsigned i)
	{
	static char Str[SIZE_32];

	double d = (double) i;
	if (i < 10000)
		snprintf(Str, SIZE_32, "%u", i);
	else if (i < 1e6)
		snprintf(Str, SIZE_32, "%.1fk", d/1e3);
	else if (i < 100e6)
		snprintf(Str, SIZE_32, "%.1fM", d/1e6);
	else if (i < 1e9)
		snprintf(Str, SIZE_32, "%.0fM", d/1e6);
	else if (i < 10e9)
		snprintf(Str, SIZE_32, "%.1fG", d/1e9);
	else if (i < 100e9)
		snprintf(Str, SIZE_32, "%.0fG", d/1e9);
	else
		snprintf(Str, SIZE_32, "%.3g", d);
	return Str;
	}

const char *Int64ToStr(uint64 i)
	{
	static char Str[SIZE_64];

	double d = (double) i;
	if (i < 10000)
		snprintf(Str, SIZE_64, "%u", (unsigned) i);
	else if (i < 1e6)
		snprintf(Str, SIZE_64, "%.1fk", d/1e3);
	else if (i < 10e6)
		snprintf(Str, SIZE_64, "%.1fM", d/1e6);
	else if (i < 1e9)
		snprintf(Str, SIZE_64, "%.0fM", d/1e6);
	else if (i < 10e9)
		snprintf(Str, SIZE_64, "%.1fG", d/1e9);
	else if (i < 100e9)
		snprintf(Str, SIZE_64, "%.0fG", d/1e9);
	else
		snprintf(Str, SIZE_64, "%.3g", d);
	return Str;
	}

const char *FloatToStr(double d)
	{
	static char Str[SIZE_32];

	double a = fabs(d);
	if (a < 0.01)
		snprintf(Str, SIZE_32, "%.3g", a);
	else if (a >= 0.01 && a < 1)
		snprintf(Str, SIZE_32, "%.3f", a);
	else if (a <= 10 && a >= 1)
		{
		double intpart;
		if (modf(a, &intpart) < 0.05)
			snprintf(Str, SIZE_32, "%.0f", d);
		else
			snprintf(Str, SIZE_32, "%.1f", d);
		}
	else if (a > 10 && a < 10000)
		snprintf(Str, SIZE_32, "%.1f", d);
	else if (a < 1e6)
		snprintf(Str, SIZE_32, "%.1fk", d/1e3);
	else if (a < 10e6)
		snprintf(Str, SIZE_32, "%.1fM", d/1e6);
	else if (a < 1e9)
		snprintf(Str, SIZE_32, "%.1fM", d/1e6);
	else if (a < 10e9)
		snprintf(Str, SIZE_32, "%.1fG", d/1e9);
	else if (a < 100e9)
		snprintf(Str, SIZE_32, "%.1fG", d/1e9);
	else
		snprintf(Str, SIZE_32, "%.3g", d);
	return Str;
	}

const char *IntFloatToStr(double d)
	{
	static char Str[SIZE_32];

	double a = fabs(d);
	if (a < 1.0)
		snprintf(Str, SIZE_32, "%.3g", a);
	else if (a <= 10)
		snprintf(Str, SIZE_32, "%.0f", d);
	else if (a > 10 && a < 10000)
		snprintf(Str, SIZE_32, "%.0f", d);
	else if (a < 1e6)
		snprintf(Str, SIZE_32, "%.1fk", d/1e3);
	else if (a < 10e6)
		snprintf(Str, SIZE_32, "%.1fM", d/1e6);
	else if (a < 1e9)
		snprintf(Str, SIZE_32, "%.1fM", d/1e6);
	else if (a < 10e9)
		snprintf(Str, SIZE_32, "%.1fG", d/1e9);
	else if (a < 100e9)
		snprintf(Str, SIZE_32, "%.1fG", d/1e9);
	else
		snprintf(Str, SIZE_32, "%.3g", d);
	return Str;
	}

static string g_CurrentProgressLine;
static string g_ProgressDesc;
static size_t g_ProgressIndex;
static size_t g_ProgressCount;

static size_t g_CurrProgressLineLength;
static size_t g_LastProgressLineLength;
static size_t g_CountsInterval;
static size_t g_StepCalls;
static time_t g_TimeLastOutputStep;

string &GetProgressPrefixStr(string &s)
	{
	double Bytes = GetMemUseBytes();
	unsigned Secs = GetElapsedSecs();
	s = string(SecsToHHMMSS(Secs));
	if (Bytes > 0)
		{
		s.push_back(' ');
		char Str[SIZE_32];
		snprintf(Str, SIZE_32, "%5s", MemBytesToStr(Bytes));
		s += string(Str);
		}
	s.push_back(' ');
	return s;
	}

const char *GetElapsedTimeStr(string &s)
	{
	unsigned Secs = GetElapsedSecs();
	s = string(SecsToHHMMSS(Secs));
	return s.c_str();
	}

const char *GetMaxRAMStr(string &s)
	{
	char Str[32];
	snprintf(Str, SIZE_32, "%5s", MemBytesToStr(g_PeakMemUseBytes));
	s = string(Str);
	return s.c_str();
	}

static bool g_ProgressPrefixOn = true;

bool ProgressPrefix(bool On)
	{
	bool OldValue = g_ProgressPrefixOn;
	g_ProgressPrefixOn = On;
	return OldValue;
	}

void ProgressLog(const char *Format, ...)
	{
	string Str;
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);

	Log("%s", Str.c_str());
	bool SavedPrefix = g_ProgressPrefixOn;
	g_ProgressPrefixOn = false;
	Progress("%s", Str.c_str());
	g_ProgressPrefixOn = SavedPrefix;
	}

void ProgressLogPrefix(const char *Format, ...)
	{
	string Str;
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);

	Log("%s", Str.c_str());
	bool SavedPrefix = g_ProgressPrefixOn;
	g_ProgressPrefixOn = true;
	Progress("%s", Str.c_str());
	g_ProgressPrefixOn = SavedPrefix;
	}

void Pr(FILE *f, const char *Format, ...)
	{
	if (f == 0)
		return;

	va_list args;
	va_start(args, Format);
	vfprintf(f, Format, args);
	va_end(args);
	}

void Progress(const char *Format, ...)
	{
	if (opt_quiet)
		return;

	string Str;
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);

#if	0
	Log("Progress(");
	for (unsigned i = 0; i < Str.size(); ++i)
		{
		char c = Str[i];
		if (c == '\r')
			Log("\\r");
		else if (c == '\n')
			Log("\\n");
		else
			Log("%c", c);
		}
	Log(")\n");
#endif //0

	for (unsigned i = 0; i < Str.size(); ++i)
		{
		if (g_ProgressPrefixOn && g_CurrProgressLineLength == 0)
			{
			string s;
			GetProgressPrefixStr(s);
			for (unsigned j = 0; j < s.size(); ++j)
				{
				fputc(s[j], stderr);
				++g_CurrProgressLineLength;
				}
			}

		char c = Str[i];
		if (c == '\n' || c == '\r')
			{
			for (size_t j = g_CurrProgressLineLength; j < g_LastProgressLineLength; ++j)
				fputc(' ', stderr);
			if (c == '\n')
				g_LastProgressLineLength = 0;
			else
				g_LastProgressLineLength = g_CurrProgressLineLength;
			g_CurrProgressLineLength = 0;
			fputc(c, stderr);
			}
		else
			{
			fputc(c, stderr);
			++g_CurrProgressLineLength;
			}
		}
	}

void LogProgramInfoAndCmdLine()
	{
	PrintProgramInfo(g_fLog);
	PrintCmdLine(g_fLog);
#ifdef	_MSC_VER
	const char *e = getenv("CYGTZ");
	if (e != 0 && strcmp(e, "YES") == 0)
		putenv("TZ=");
#endif
	time_t Now = time(0);
	struct tm *t = localtime(&Now);
	const char *s = asctime(t);
	Log("Started %s", s); // there is a newline in s
	}

void LogElapsedTimeAndRAM()
	{
	time_t Now = time(0);
	struct tm *t = localtime(&Now);
	const char *s = asctime(t);
	unsigned Secs = GetElapsedSecs();

	Log("\n");
	Log("Finished %s", s); // there is a newline in s
	Log("Elapsed time %s\n", SecsToHHMMSS((int) Secs));
	Log("Max memory %s\n", MemBytesToStr(g_PeakMemUseBytes));
#if	WIN32 && DEBUG
// Skip exit(), which can be very slow in DEBUG build
// VERY DANGEROUS practice, because it skips global destructors.
// But if you know the rules, you can break 'em, right?
	ExitProcess(0);
#endif
	}

const char *PctStr(double x, double y)
	{
	if (y == 0)
		{
		if (x == 0)
			return "100%";
		else
			return "inf%";
		}
	static char Str[SIZE_16];
	double p = x*100.0/y;
	snprintf(Str, SIZE_16, "%5.1f%%", p);
	return Str;
	}

string &GetProgressLevelStr(string &s)
	{
	uint Index = (uint) g_ProgressIndex;
	uint Count = (uint) g_ProgressCount;
	if (Count == UINT_MAX)
		{
		if (Index == UINT_MAX)
			s = "100%";
		else
			{
			char Tmp[SIZE_16];
			snprintf(Tmp, SIZE_16, "%u", Index); 
			s = Tmp;
			}
		}
	else
		s = string(PctStr(Index+1, Count));
	s += string(" ") + g_ProgressDesc;
	return s;
	}

static const char *DefaultPCB()
	{
	return "Processing";
	}
static FN_PROGRESS_CALLBACK g_PCB = DefaultPCB;

void SetPCB(FN_PROGRESS_CALLBACK PCB)
	{
	g_PCB = PCB;
	}

static FILE *g_fProg;
static double g_ProgFileSize;
static unsigned g_ProgFileTick;
static const char *g_ProgFileMsg = "Processing";
static string g_ProgFileStr;

void ProgressFileInit(FILE *f, const char *Format, ...)
	{
	g_fProg = f;
	g_ProgFileSize = (double) GetStdioFileSize64(f);
	g_ProgFileTick = 0;
	if (Format == 0)
		g_ProgFileMsg = "Processing";
	else
		{
		va_list ArgList;
		va_start(ArgList, Format);
		myvstrprintf(g_ProgFileStr, Format, ArgList);
		va_end(ArgList);
		g_ProgFileMsg = g_ProgFileStr.c_str();
		}
	ProgressStep(0, 1002, "%s", g_ProgFileMsg);
	}

void ProgressFileStep(const char *Format, ...)
	{
	double Pos = (double) GetStdioFilePos64(g_fProg);
	unsigned Tick = (unsigned) ((Pos*1000.0)/g_ProgFileSize);
	if (Tick <= g_ProgFileTick)
		return;
	if (Format != 0)
		{
		va_list ArgList;
		va_start(ArgList, Format);
		myvstrprintf(g_ProgFileStr, Format, ArgList);
		va_end(ArgList);
		g_ProgFileMsg = g_ProgFileStr.c_str();
		}
	ProgressStep(Tick, 1002, "%s", g_ProgFileMsg);
	g_ProgFileTick = Tick;
	}

void ProgressFileDone(const char *Format, ...)
	{
	if (Format != 0)
		{
		va_list ArgList;
		va_start(ArgList, Format);
		myvstrprintf(g_ProgFileStr, Format, ArgList);
		va_end(ArgList);
		g_ProgFileMsg = g_ProgFileStr.c_str();
		}
	ProgressStep(1001, 1002, "%s", g_ProgFileMsg);
	}

#if TIMING
static time_t g_LastLogTimerSecs; 
#endif

static void ProgressStepHook(time_t Now)
	{
#if TIMING
	if (opt_log_timer_mins > 0)
		{
		if (g_LastLogTimerSecs == 0)
			g_LastLogTimerSecs = Now;
		else
			{
			//Log("Now %.0f Last %.0f Diff %.0f opt %.0f\n",
			//  float(g_TimeLastOutputStep),
			//  float(g_LastLogTimerSecs),
			//  float(g_TimeLastOutputStep - g_LastLogTimerSecs)/60,
			//  float(opt_log_timer_mins));
			if ((g_TimeLastOutputStep - g_LastLogTimerSecs)/60 >= opt_log_timer_mins)
				{
				void LogTiming();
				LogTiming();
				g_LastLogTimerSecs = Now;
				}
			}
		}
#endif
	}

void ProgressCallback(unsigned i, unsigned N)
	{
	if (opt_quiet)
		return;

	if (i == 0)
		{
		g_ProgressIndex = 0;
		g_ProgressCount = N;
		g_CountsInterval = 1;
		g_StepCalls = 0;
		g_TimeLastOutputStep = 0;
		if (g_CurrProgressLineLength > 0)
			Progress("\n");
		}

	bool IsLastStep = (i == UINT_MAX || i + 1 == N);
	if (!IsLastStep)
		{
		++g_StepCalls;
		if (g_StepCalls%g_CountsInterval != 0)
			return;

		time_t Now = time(0);
		if (Now == g_TimeLastOutputStep)
			{
			if (g_CountsInterval < 128)
				g_CountsInterval = (g_CountsInterval*3)/2;
			else
				g_CountsInterval += 64;
			return;
			}
		else
			{
			time_t Secs = Now - g_TimeLastOutputStep;
			if (Secs > 1)
				g_CountsInterval = unsigned(g_CountsInterval/(Secs*8));
			}

		if (g_CountsInterval < 1)
			g_CountsInterval = 1;

		g_TimeLastOutputStep = Now;
		}

	g_ProgressIndex = i;

	ProgressStepHook(g_TimeLastOutputStep);

	Progress(" %s", PctStr(i+1, N));
	Progress(" %s\r", (g_PCB)());

	if (IsLastStep)
		{
		g_CountsInterval = 1;
		fputc('\n', stderr);
		}
	}

void ProgressStep(size_t i, size_t N, const char *Format, ...)
	{
	if (opt_quiet)
		return;

	if (i == 0)
		{
		string Str;
		va_list ArgList;
		va_start(ArgList, Format);
		myvstrprintf(Str, Format, ArgList);
		va_end(ArgList);
		g_ProgressDesc = Str;
		g_ProgressIndex = 0;
		g_ProgressCount = N;
		g_CountsInterval = 1;
		g_StepCalls = 0;
		g_TimeLastOutputStep = 0;
		if (g_CurrProgressLineLength > 0)
			Progress("\n");
		}

	assert(N == g_ProgressCount);
	if (i >= N && i != UINT_MAX)
		Die("ProgressStep(%u,%u)", i, N);
	bool IsLastStep = (i == UINT_MAX || i + 1 == N);
	if (!IsLastStep)
		{
		++g_StepCalls;
		if (g_StepCalls%g_CountsInterval != 0)
			return;

		time_t Now = time(0);
		if (Now == g_TimeLastOutputStep)
			{
			if (g_CountsInterval < 128)
				g_CountsInterval = (g_CountsInterval*3)/2;
			else
				g_CountsInterval += 64;
			return;
			}
		else
			{
			time_t Secs = Now - g_TimeLastOutputStep;
			if (Secs > 1)
				g_CountsInterval = unsigned(g_CountsInterval/(Secs*8));
			}

		if (g_CountsInterval < 1)
			g_CountsInterval = 1;

		g_TimeLastOutputStep = Now;
		}

	g_ProgressIndex = i;

	if (i > 0)
		{
		va_list ArgList;
		va_start(ArgList, Format);
		myvstrprintf(g_ProgressDesc, Format, ArgList);
		}

	ProgressStepHook(g_TimeLastOutputStep);

	string LevelStr;
	GetProgressLevelStr(LevelStr);
	Progress(" %s\r", LevelStr.c_str());

	if (IsLastStep)
		{
		g_CountsInterval = 1;
		fputc('\n', stderr);
		}
	}

vector<string> g_Argv;

#define FLAG_OPT(Name)						bool opt_##Name; bool optset_##Name;
#define UNS_OPT(Name, Default, Min, Max)	unsigned opt_##Name; bool optset_##Name;
#define FLT_OPT(Name, Default, Min, Max)	double opt_##Name; bool optset_##Name;
#define STR_OPT(Name)						const char *opt_##Name = ""; bool optset_##Name;
#include "myopts.h"

static void CmdLineErr(const char *Format, ...)
	{
	va_list ArgList;
	va_start(ArgList, Format);
	string Str;
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);
	fprintf(stderr, "\n");
	fprintf(stderr, "Invalid command line\n");
	fprintf(stderr, "%s\n", Str.c_str());
	void PrintHelp();
//	PrintHelp();
	asserta(false);
	exit(1);
	}

static void GetArgsFromFile(const string &FileName, vector<string> &Args)
	{
	Args.clear();

	FILE *f = OpenStdioFile(FileName);
	string Line;
	while (ReadLineStdioFile(f, Line))
		{
		size_t n = Line.find('#');
		if (n != string::npos)
			Line = Line.substr(0, n);
		vector<string> Fields;
		Split(Line, Fields, 0);
		Args.insert(Args.end(), Fields.begin(), Fields.end());
		}
	CloseStdioFile(f);
	}

static bool TryFlagOpt(const char *OptName)
	{
#define FLAG_OPT(Name)	if (strcmp(OptName, #Name) == 0) { opt_##Name = true; optset_##Name = true; return true; }
#define UNS_OPT(Name, Default, Min, Max)	/* empty */
#define FLT_OPT(Name, Default, Min, Max)	/* empty */
#define STR_OPT(Name)						/* empty */

#include "myopts.h"
	return false;
	}

static bool TryUnsOpt(const char *OptName, const char *Value)
	{
#define UNS_OPT(Name, Default, Min, Max)	if (strcmp(OptName, #Name) == 0) { opt_##Name = atou(Value); optset_##Name = true; return true; }
#define FLAG_OPT(Name)						/* empty */
#define FLT_OPT(Name, Default, Min, Max)	/* empty */
#define STR_OPT(Name)						/* empty */
#include "myopts.h"
	return false;
	}

static bool TryFloatOpt(const char *OptName, const char *Value)
	{
#define FLT_OPT(Name, Default, Min, Max)	if (strcmp(OptName, #Name) == 0) { opt_##Name = StrToFloat(Value); optset_##Name = true; return true; }
#define UNS_OPT(Name, Default, Min, Max)	/* empty */
#define FLAG_OPT(Name)						/* empty */
#define STR_OPT(Name)						/* empty */
#include "myopts.h"
	return false;
	}

static bool TryStrOpt(const char *OptName, const char *Value)
	{
#define STR_OPT(Name)	if (strcmp(OptName, #Name) == 0) { opt_##Name = mystrsave(Value); optset_##Name = true; return true; }
#define UNS_OPT(Name, Default, Min, Max)	/* empty */
#define FLT_OPT(Name, Default, Min, Max)	/* empty */
#define FLAG_OPT(Name)						/* empty */
#include "myopts.h"
	return false;
	}

void MyCmdLine(int argc, char **argv)
	{
	setbuf(stdout, 0);
	setbuf(stderr, 0);

	if (argc == 1 ||
	  argc == 2 && !strcmp(argv[1], "-h") ||
	  argc == 2 && !strcmp(argv[1], "--help"))
		Help();

#define UNS_OPT(Name, Default, Min, Max)	opt_##Name = Default;
#define FLT_OPT(Name, Default, Min, Max)	opt_##Name = Default;
#define FLAG_OPT(Name)						/* empty */
#define STR_OPT(Name)						/* empty */
#include "myopts.h"

	for (unsigned i = 0; i < (unsigned) argc; )
		{
		const string &Arg = argv[i];
		if (Arg == "file:" && i + 1 < (unsigned) argc)
			{
			const string &FileName = argv[i+1];
			vector<string> Args;
			GetArgsFromFile(FileName, Args);
			for (unsigned k = 0; k < SIZE(Args); ++k)
				g_Argv.push_back(Args[k]);
			i += 2;
			}
		else
			{
			g_Argv.push_back(Arg);
			i += 1;
			}
		}

	const unsigned ArgCount = SIZE(g_Argv);
	unsigned ArgIndex = 1;
	for (;;)
		{
		if (ArgIndex >= ArgCount)
			break;
		const string &Arg = g_Argv[ArgIndex];
		if (Arg.size() > 1 && Arg[0] == '-')
			{
			string LongName = (Arg.size() > 2 && Arg[1] == '-' ? Arg.substr(2) : Arg.substr(1));
			if (LongName == "help")
				Help();
			else if (LongName == "version")
				{
				Version(stdout);
				printf("\n");
				exit(0);
				}

			bool IsFlag = TryFlagOpt(LongName.c_str());
			if (IsFlag)
				{
				++ArgIndex;
				continue;
				}

			++ArgIndex;
			if (ArgIndex >= ArgCount)
				CmdLineErr("Missing value or invalid option -%s", LongName.c_str());

			const char *Value = g_Argv[ArgIndex].c_str();

			bool IsUns = TryUnsOpt(LongName.c_str(), Value);
			if (IsUns)
				{
				++ArgIndex;
				continue;
				}

			bool IsFloat = TryFloatOpt(LongName.c_str(), Value);
			if (IsFloat)
				{
				++ArgIndex;
				continue;
				}

			bool IsStr = TryStrOpt(LongName.c_str(), Value);
			if (IsStr)
				{
				++ArgIndex;
				continue;
				}
			
			CmdLineErr("Unknown option %s", LongName.c_str());
			}
		else
			CmdLineErr("Expected -option_name or --option_name, got '%s'", Arg.c_str());
		}

#if	TIMING
	if (opt_threads > 1)
		Die("--threads > 1 && TIMING");
#endif

	if (opt_compilerinfo)
		{
		CompilerInfo();
		exit(0);
		}
	SetLogFileName(opt_log);
	}

static unsigned GetStructPack()
	{
	struct
		{
		char a;
		char b;
		} x;

	return (unsigned) (&x.b - &x.a);
	}

static void CompilerInfo()
	{
	printf("%u bits\n", BITS);

#ifdef __GNUC__
	printf("__GNUC__\n");
#endif

#ifdef __APPLE__
	printf("__APPLE__\n");
#endif

#ifdef _MSC_VER
	printf("_MSC_VER %d\n", _MSC_VER);
#endif

#define x(t)	printf("sizeof(" #t ") = %d\n", (int) sizeof(t));
	x(int)
	x(long)
	x(float)
	x(double)
	x(void *)
	x(off_t)
	x(size_t)
#undef x

	printf("pack(%u)\n", GetStructPack());

#ifdef _FILE_OFFSET_BITS
    printf("_FILE_OFFSET_BITS = %d\n", _FILE_OFFSET_BITS);
#else
    printf("_FILE_OFFSET_BITS not defined\n");
#endif

	exit(0);
	}

bool EndsWith(const string &s, const string &t)
	{
	size_t ns = s.size();
	size_t nt = t.size();
	if (nt > ns)
		return false;

	for (uint i = 0; i < nt; ++i)
		{
		char tc = t[i];
		char sc = s[ns - nt + i];
		if (tc != sc)
			return false;
		}
	return true;
	}

bool StartsWith(const char *S, const char *T)
	{
	for (;;)
		{
		char t = *T++;
		if (t == 0)
			return true;
		char s = *S++;
		if (s != t)
			return false;
		}
	}

bool StartsWith(const string &S, const char *T)
	{
	return StartsWith(S.c_str(), T);
	}

bool StartsWith(const string &s, const string &t)
	{
	return StartsWith(s.c_str(), t.c_str());
	}

void SplitWhite(const string &Str, vector<string> &Fields)
	{
	Fields.clear();
	const unsigned Length = (unsigned) Str.size();
	string s;
	string Field;
	for (unsigned i = 0; i < Length; ++i)
		{
		char c = Str[i];
		if (isspace(c))
			{
			if (!Field.empty())
				Fields.push_back(Field);
			Field.clear();
			}
		else
			Field.push_back(c);
		}
	if (!Field.empty())
		Fields.push_back(Field);
	}

void Split(const string &Str, vector<string> &Fields, char Sep)
	{
	Fields.clear();
	const unsigned Length = (unsigned) Str.size();
	string s;
	for (unsigned i = 0; i < Length; ++i)
		{
		char c = Str[i];
		if ((Sep == 0 && isspace(c)) || c == Sep)
			{
			if (!s.empty() || Sep != 0)
				Fields.push_back(s);
			s.clear();
			}
		else
			s.push_back(c);
		}
	if (!s.empty())
		Fields.push_back(s);
	}

const char *g_GitVer = 
#include "gitver.txt"
		;

void Version(FILE *f)
	{
	if (f == 0)
		return;
	const char *Flags = ""
#if	DEBUG
	"D"
#endif
#if TIMING
	"T"
#endif
	;
	fprintf(f, "\n");

	fprintf(f, "reseek v%s.%s%s [%s]\n", MY_VERSION, GetPlatform(), Flags, g_GitVer);
	}

void PrintHelp()
	{
	extern const char *usage_txt[];
	extern int g_n_usage_txt;
	for (int i = 0; i < g_n_usage_txt; ++i)
		puts(usage_txt[i]);
	puts("");
	}

void Help()
	{
	PrintProgramInfo(stdout);
	PrintCopyright(stdout);
	PrintHelp();
	exit(0);
	}

void PrintProgramInfo(FILE *f)
	{
	if (f == 0)
		return;
	Version(f);
	}

void PrintCopyright(FILE *f)
	{
	if (f == 0)
		return;

	fprintf(f, "(C) Copyright 2024 Robert C. Edgar\n");
	fprintf(f, "\n");
	}

void PrintCmdLine(FILE *f)
	{
	if (f == 0)
		return;
	for (unsigned i = 0; i < SIZE(g_Argv); ++i)
		fprintf(f, "%s ", g_Argv[i].c_str());
	fprintf(f, "\n");
	}

void GetCmdLine(string &s)
	{
	s.clear();
	for (unsigned i = 0; i < SIZE(g_Argv); ++i)
		{
		if (i > 0)
			s += " ";
		s += g_Argv[i];
		}
	}

char *mystrsave(const char *s)
	{
	unsigned n = unsigned(strlen(s));
	char *t = myalloc(char, n+1);
	memcpy(t, s, n+1);
	return t;
	}

unsigned myipow(unsigned x, unsigned y)
	{
	unsigned result = 1;
	for (unsigned k = 0; k < y; ++k)
		{
		if (result > UINT_MAX/x)
			Die("myipow(%u, %u), overflow", x, y);
		result *= x;
		}
	return result;
	}

void LogInt(unsigned i, unsigned w)
    {
    if (w == UINT_MAX)
        {
        if (i < 9999)
            Log("%u", i);
        else
            Log("%u (%s)", i, IntToStr(i));
        }
    else
        {
        if (i < 9999)
            Log("%*u", w, i);
        else
            Log("%*u (%s)", w, i, IntToStr(i));
        }
    }

void Logu(unsigned u, unsigned w, unsigned prefixspaces)
	{
	for (unsigned i = 0; i < prefixspaces; ++i)
		Log(" ");
	if (u == UINT_MAX)
		Log("%*.*s", w, w, "*");
	else
		Log("%*u", w, u);
	}

void Logf(float x, unsigned w, unsigned prefixspaces)
	{
	for (unsigned i = 0; i < prefixspaces; ++i)
		Log(" ");
	if (x == FLT_MAX)
		Log("%*.*s", w, w, "*");
	else
		Log("%*.2f", w, x);
	}

static uint32 g_SLCG_state = 1;

// Numerical values used by Microsoft C, according to wikipedia:
// http://en.wikipedia.org/wiki/Linear_congruential_generator
static uint32 g_SLCG_a = 214013;
static uint32 g_SLCG_c = 2531011;

// Simple Linear Congruential Generator
// Bad properties; used just to initialize the better generator.
static uint32 SLCG_rand()
	{
	g_SLCG_state = g_SLCG_state*g_SLCG_a + g_SLCG_c;
	return g_SLCG_state;
	}

static void SLCG_srand(uint32 Seed)
	{
	g_SLCG_state = Seed;
	for (int i = 0; i < 10; ++i)
		SLCG_rand();
	}

/***
A multiply-with-carry random number generator, see:
http://en.wikipedia.org/wiki/Multiply-with-carry

The particular multipliers used here were found on
the web where they are attributed to George Marsaglia.
***/

static bool g_InitRandDone = false;
static uint32 g_X[5];

uint32 RandInt32()
	{
	InitRand();

	uint64 Sum = 2111111111*(uint64) g_X[3] + 1492*(uint64) g_X[2] +
	  1776*(uint64) g_X[1] + 5115*(uint64) g_X[0] + g_X[4];
	g_X[3] = g_X[2];
	g_X[2] = g_X[1];
	g_X[1] = g_X[0];
	g_X[4] = (uint32) (Sum >> 32);
	g_X[0] = (uint32) Sum;
	return g_X[0];
	}

unsigned randu32()
	{
	return (unsigned) RandInt32();
	}

void InitRand()
	{
	if (g_InitRandDone)
		return;
// Do this first to avoid recursion
	g_InitRandDone = true;

	unsigned Seed = (optset_randseed ? opt_randseed : (unsigned) (time(0)*getpid()));
	Log("RandSeed=%u\n", Seed);
	SLCG_srand(Seed);

	for (unsigned i = 0; i < 5; i++)
		g_X[i] = SLCG_rand();

	for (unsigned i = 0; i < 100; i++)
		RandInt32();
	}

unsigned GetCPUCoreCount()
	{
#ifdef _MSC_VER
	SYSTEM_INFO SI;
	GetSystemInfo(&SI);
	unsigned n = SI.dwNumberOfProcessors;
	if (n == 0 || n > 64)
		return 1;
	return n;
#else
	long n = sysconf(_SC_NPROCESSORS_ONLN);
	if (n <= 0)
		return 1;
	return (unsigned) n;
#endif
	}

// MUST COME AT END BECAUSE OF #undefs
#undef myalloc
#undef myfree

#if	RCE_MALLOC
#undef mymalloc
#undef myfree
#undef myfree2

static unsigned g_NewCalls;
static unsigned g_FreeCalls;
static double g_InitialMemUseBytes;
static double g_TotalAllocBytes;
static double g_TotalFreeBytes;
static double g_NetBytes;
static double g_MaxNetBytes;

void LogAllocStats()
	{
	Log("\n");
	Log("       Allocs  %u\n", g_NewCalls);
	Log("        Frees  %u\n", g_FreeCalls);
	Log("Initial alloc  %s\n", MemBytesToStr(g_InitialMemUseBytes));
	Log("  Total alloc  %s\n", MemBytesToStr(g_TotalAllocBytes));
	Log("   Total free  %s\n", MemBytesToStr(g_TotalFreeBytes));
	Log("    Net bytes  %s\n", MemBytesToStr(g_NetBytes));
	Log("Max net bytes  %s\n", MemBytesToStr(g_MaxNetBytes));
	Log("   Peak total  %s\n", MemBytesToStr(g_MaxNetBytes + g_InitialMemUseBytes));
	}

void *mymalloc(unsigned bytes, const char *FileName, int Line)
	{
	void *rce_malloc(unsigned bytes, const char *FileName, int Line);
	return rce_malloc(bytes, FileName, Line);
	}

void myfree(void *p, const char *FileName, int Line)
	{
	void rce_free(void *p, const char *FileName, int Line);
	rce_free(p, FileName, Line);
	}

void myfree2(void *p, unsigned bytes, const char *FileName, int Line)
	{
	void rce_free(void *p, const char *FileName, int Line);
	rce_free(p, FileName, Line);
	}

#else // RCE_MALLOC

void *mymalloc(size_t bytes)
	{
	void *p = malloc(bytes);
	if (0 == p)
		{
		double b = GetMemUseBytes();
		fprintf(stderr, "\nOut of memory mymalloc(%u), curr %.3g bytes",
		  (unsigned) bytes, b);
#if DEBUG && defined(_MSC_VER)
		asserta(_heapchk() == _HEAPOK);
#endif
		Die("Out of memory, mymalloc(%u), curr %.3g bytes\n",
		  (unsigned) bytes, b);
		}
	return p;
	}

void myfree(void *p)
	{
	if (p == 0)
		return;
	free(p);
	}

#endif // RCE_MALLOC

void Ps(string &Str, const char *Format, ...)
	{
	va_list ArgList;
	va_start(ArgList, Format);
	myvstrprintf(Str, Format, ArgList);
	va_end(ArgList);
	}

void Pf(FILE *f, const char *Format, ...)
	{
	if (f == 0)
		return;
	va_list ArgList;
	va_start(ArgList, Format);
	string Tmp;
	myvstrprintf(Tmp, Format, ArgList);
	va_end(ArgList);
	fprintf(f, "%s", Tmp.c_str());
	}

void Psa(string &Str, const char *Format, ...)
	{
	va_list ArgList;
	va_start(ArgList, Format);
	string Tmp;
	myvstrprintf(Tmp, Format, ArgList);
	va_end(ArgList);
	Str += Tmp;
	}

unsigned GetRequestedThreadCount()
	{
	static unsigned N = 1;
	static bool Done = false;
	if (Done)
		return N;
	unsigned CoreCount = GetCPUCoreCount();
	bool MsgDone = false;
	if (optset_threads)
		N = opt_threads;
	else
		N = CoreCount;
	if (N == 0)
		N = 1;
	Done = true;
	if (!MsgDone)
		{
		Progress("Starting %u threads (%u CPU cores)\n", N, CoreCount);
		MsgDone = true;
		}
	return N;
	}

void ToLower(string &Str)
	{
	unsigned n = SIZE(Str);
	for (uint i = 0; i < n; ++i)
		Str[i] = tolower(Str[i]);
	}

void ToUpper(string &Str)
	{
	unsigned n = SIZE(Str);
	for (uint i = 0; i < n; ++i)
		Str[i] = toupper(Str[i]);
	}

void StripWhiteSpace(string &Str)
	{
	unsigned n = SIZE(Str);
	unsigned FirstNonWhite = UINT_MAX;
	unsigned LastNonWhite = UINT_MAX;
	for (unsigned i = 0; i < n; ++i)
		{
		char c = Str[i];
		if (!isspace(c))
			{
			if (FirstNonWhite == UINT_MAX)
				FirstNonWhite = i;
			LastNonWhite = i;
			}
		}

	if (FirstNonWhite == UINT_MAX)
		return;

	string t;
	for (unsigned i = FirstNonWhite; i <= LastNonWhite; ++i)
		{
		char c = Str[i];
		t += c;
		}
	Str = t;
	}

const char *GetPctStr(double x, double y, string &s)
	{
	double Pct = GetPct(x, y);
	if (Pct < 0.1)
		Ps(s, "%.3g%%", Pct);
	else if (Pct < 1)
		Ps(s, "%.2f%%", Pct);
	else
		Ps(s, "%.1f%%", Pct);
	return s.c_str();
	}

// Fisher-Yates shuffle:
// To shuffle an array a of n elements (indices 0 .. n-1):
//  for i from n - 1 downto 1 do
//       j := random integer with 0 <= j <= i
//       exchange a[j] and a[i]
void Shuffle(vector<unsigned> &v)
	{
	const unsigned N = SIZE(v);
	for (unsigned i = N - 1; i >= 1; --i)
		{
		unsigned j = randu32()%(i + 1);
		
		unsigned vi = v[i];
		unsigned vj = v[j];

		v[i] = vj;
		v[j] = vi;
		}
	}

void ReadLinesFromFile(const string &FileName, vector<string> &Lines)
	{
	if (EndsWith(FileName, ".gz"))
		{
		void ReadLinesFromGzipFile(const string &FileName, vector<string> &Lines);
		ReadLinesFromGzipFile(FileName, Lines);
		return;
		}
	Lines.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	while (ReadLineStdioFile(f, Line))
		Lines.push_back(Line);
	CloseStdioFile(f);
	}


void Dirize(string &Dir)
	{
	if (!EndsWith(Dir, "/") && !EndsWith(Dir, "\\"))
		Dir += "/";
	}
