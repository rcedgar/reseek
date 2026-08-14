#include "myutils.h"

#include <thread>
#include <atomic>

void LogCoords16(const char *mem, uint chainLength);
float *GetCoordsFromMem(const char *mem, uint chainLength, uint entryLength);

//////////////////////////////////////////////////////////////////
// "C:\src\notebooks\2026-07-02_foldseek_db_reverse_engineer.docx"
//////////////////////////////////////////////////////////////////

/***
-rwxrwxrwx 1 bob bob  35789 Feb  2 18:23 hiqual				# aa sequences, ASCII
-rwxrwxrwx 1 bob bob      4 Feb  2 18:23 hiqual.dbtype		# 0x00000000
-rwxrwxrwx 1 bob bob    492 Feb  2 18:23 hiqual.index		# index to aa sequences (tsv)
-rwxrwxrwx 1 bob bob    409 Feb  2 18:23 hiqual.lookup		# e.g. 23      6f5p    23
-rwxrwxrwx 1 bob bob    302 Feb  2 18:23 hiqual.source		# e.g. 23      6f5p
-rwxrwxrwx 1 bob bob 234772 Feb  2 18:23 hiqual_ca			# C-alpha coords
-rwxrwxrwx 1 bob bob      4 Feb  2 18:23 hiqual_ca.dbtype	# 0x00000065
-rwxrwxrwx 1 bob bob    552 Feb  2 18:23 hiqual_ca.index	# index to C-alphas (tsv)
-rwxrwxrwx 1 bob bob    356 Feb  2 18:23 hiqual_h			# labels ("headers"), ASCII
-rwxrwxrwx 1 bob bob      4 Feb  2 18:23 hiqual_h.dbtype	# 0x0000000c
-rwxrwxrwx 1 bob bob    329 Feb  2 18:23 hiqual_h.index		# index to labels (tsv)
-rwxrwxrwx 1 bob bob  35789 Feb  2 18:23 hiqual_ss			# 3Di sequences, ASCII
-rwxrwxrwx 1 bob bob      4 Feb  2 18:23 hiqual_ss.dbtype	# 0x00000000
-rwxrwxrwx 1 bob bob    492 Feb  2 18:23 hiqual_ss.index	# index to 3Di
***/

static void CheckDBType(const string &Prefix,
						const string &Suffix,
						uint32_t ExpectedType)
	{
	const string FN(Prefix + Suffix);
	FILE *f = OpenStdioFile(FN);
	uint32_t Size = GetStdioFileSize32(f);
	if (Size != 4)
		{
		Warning("Size is %u bytes, expected 4: %s",
				Size, FN.c_str());
		CloseStdioFile(f);
		return;
		}
	uint32_t Type;
	ReadStdioFile(f, &Type, 4);
	if (Type != ExpectedType)
		Warning("Type is 0x%X, expected 0x%X", Type, ExpectedType);
	CloseStdioFile(f);
	}

static void ReadIndex(const string &FN,
					  vector<uint64> &Offsets,
					  vector<uint64> &Lengths)
	{
	Offsets.clear();
	Lengths.clear();
	FILE *f = OpenStdioFile(FN);
	ProgressFileInit(f, "Reading %s", FN.c_str());
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep();
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 3)
			Die("Expected 3 fields, got '%s' in %s",
				Line.c_str(), FN.c_str());
		Offsets.push_back(StrToUint64(Fields[1]));
		Lengths.push_back(StrToUint64(Fields[2]));
		}
	ProgressFileDone();
	CloseStdioFile(f);
	}

static void ReadRecord(FILE *f, uint64 Off, uint64 Len, string &s)
	{
	if (Len == 0)
		{
		s.clear();
		return;
		}
	asserta(Len <= UINT32_MAX);
	s.resize(size_t(Len));
	ReadStdioFile64(f, Off, &s[0], Len);
	while (!s.empty())
		{
		char c = s.back();
		if (c == 0 || c == '\n' || c == '\r')
			s.pop_back();
		else
			break;
		}
	}

static void WriteCal(FILE *f, const string &Label, const string &AA,
	const float *Coords, uint L)
	{
	if (f == 0)
		return;
	string buf;
	buf.reserve(Label.size() + 2 + size_t(L)*40);
	buf.push_back('>');
	buf += Label;
	buf.push_back('\n');
	char line[64];
	for (uint i = 0; i < L; ++i)
		{
		int n = snprintf(line, sizeof(line), "%c\t%.1f\t%.1f\t%.1f\n",
			AA[i], Coords[i], Coords[L + i], Coords[2*L + i]);
		asserta(n > 0 && n < (int) sizeof(line));
		buf.append(line, size_t(n));
		}
	WriteStdioFile64(f, buf.data(), buf.size());
	}

static void ConcatShardFiles(const string &dest, const vector<string> &shards)
	{
	FILE *fout = CreateStdioFile(dest);
	const uint BUF = 1u << 20;
	byte *buf = myalloc(byte, BUF);
	for (uint i = 0; i < SIZE(shards); ++i)
		{
		FILE *fin = OpenStdioFile(shards[i]);
		uint64 sz = GetStdioFileSize64(fin);
		uint64 pos = 0;
		while (pos < sz)
			{
			uint64 n = sz - pos;
			if (n > BUF)
				n = BUF;
			ReadStdioFile64(fin, pos, buf, n);
			WriteStdioFile64(fout, buf, n);
			pos += n;
			}
		CloseStdioFile(fin);
		DeleteStdioFile(shards[i]);
		}
	CloseStdioFile(fout);
	myfree(buf);
	}

struct fsdb_ctx
	{
	uint N = 0;
	bool want_fasta = false;
	bool want_3di = false;
	bool want_cal = false;
	const vector<uint64> *off_h = 0;
	const vector<uint64> *len_h = 0;
	const vector<uint64> *off_aa = 0;
	const vector<uint64> *len_aa = 0;
	const vector<uint64> *off_ss = 0;
	const vector<uint64> *len_ss = 0;
	const vector<uint64> *off_ca = 0;
	const vector<uint64> *len_ca = 0;
	FILE *fH = 0;
	FILE *fAA = 0;
	FILE *fSS = 0;
	FILE *fCA = 0;
	uint64 ca_size = 0;
	vector<FILE *> fasta_fs;
	vector<FILE *> di_fs;
	vector<FILE *> cal_fs;
	atomic<uint> done;
	};

static void ConvertRange(fsdb_ctx *ptrC, uint tid, uint begin, uint end)
	{
	fsdb_ctx &C = *ptrC;
	FILE *faa = C.want_fasta ? C.fasta_fs[tid] : 0;
	FILE *f3Di = C.want_3di ? C.di_fs[tid] : 0;
	FILE *fcal = C.want_cal ? C.cal_fs[tid] : 0;

	string Label;
	string AA;
	string SS;
	time_t LastTime = 0;
	for (uint SeqIdx = begin; SeqIdx < end; ++SeqIdx)
		{
		ReadRecord(C.fH, (*C.off_h)[SeqIdx], (*C.len_h)[SeqIdx], Label);

		if (C.want_fasta || C.want_cal)
			ReadRecord(C.fAA, (*C.off_aa)[SeqIdx], (*C.len_aa)[SeqIdx], AA);

		if (C.want_3di)
			{
			ReadRecord(C.fSS, (*C.off_ss)[SeqIdx], (*C.len_ss)[SeqIdx], SS);
			if ((C.want_fasta || C.want_cal) && SIZE(AA) != SIZE(SS))
				Die("aa/3Di sequence mismatch %u, %u >%s",
					SIZE(AA), SIZE(SS), Label.c_str());
			}

		if (faa != 0)
			SeqToFasta(faa, Label.c_str(), AA.c_str(), SIZE(AA));
		if (f3Di != 0)
			SeqToFasta(f3Di, Label.c_str(), SS.c_str(), SIZE(SS));

		if (fcal != 0)
			{
			const uint L = SIZE(AA);
			const uint64 CoordsOffset = (*C.off_ca)[SeqIdx];
			const uint64 CoordsLength = (*C.len_ca)[SeqIdx];
			if (CoordsOffset + CoordsLength > C.ca_size)
				Die("CA entry %u extends past end of file", SeqIdx);
			if (CoordsLength == 0 || CoordsLength > UINT_MAX)
				Die("Invalid CA entry length %llu seq %u",
					(unsigned long long) CoordsLength, SeqIdx);
			byte *Entry = myalloc(byte, (uint) CoordsLength);
			ReadStdioFile64(C.fCA, CoordsOffset, Entry, CoordsLength);
			float *Coords = GetCoordsFromMem((const char *) Entry, L,
				(uint) CoordsLength);
			WriteCal(fcal, Label, AA, Coords, L);
			if ((void *) Coords != (void *) Entry)
				myfree((void *) Coords);
			myfree(Entry);
			}

		C.done.fetch_add(1);
		if (tid == 0)
			{
			time_t Now = time(0);
			if (Now != LastTime)
				{
				uint d = C.done.load();
				Progress("%u / %u converted (%.1f%%)\r",
					d, C.N, GetPct(d, C.N));
				LastTime = Now;
				}
			}
		}
	}

void cmd_convert_foldseekdb()
	{
	const bool want_fasta = optset_fasta;
	const bool want_3di = optset_3di;
	const bool want_cal = optset_cal;
	if (!want_fasta && !want_3di && !want_cal)
		Die("Must set -fasta, -3di and/or -cal");

	const string &Prefix = g_Arg1;

	CheckDBType(Prefix, "_h.dbtype", 0xC);
	if (want_fasta || want_cal)
		CheckDBType(Prefix, ".dbtype", 0x0);
	if (want_3di)
		CheckDBType(Prefix, "_ss.dbtype", 0x0);
	if (want_cal)
		CheckDBType(Prefix, "_ca.dbtype", 0x65);

	vector<uint64> off_h, len_h, off_aa, len_aa, off_ss, len_ss, off_ca, len_ca;
	ReadIndex(Prefix + "_h.index", off_h, len_h);
	const uint N = SIZE(off_h);
	if (N == 0)
		Die("No entries in %s_h.index", Prefix.c_str());
	ProgressLog("%u entries in '%s_h.index'\n", N, Prefix.c_str());

	if (want_fasta || want_cal)
		{
		ReadIndex(Prefix + ".index", off_aa, len_aa);
		if (SIZE(off_aa) != N)
			Die("%u headers, %u aa index entries", N, SIZE(off_aa));
		}
	if (want_3di)
		{
		ReadIndex(Prefix + "_ss.index", off_ss, len_ss);
		if (SIZE(off_ss) != N)
			Die("%u headers, %u 3Di index entries", N, SIZE(off_ss));
		}
	if (want_cal)
		{
		ReadIndex(Prefix + "_ca.index", off_ca, len_ca);
		if (SIZE(off_ca) != N)
			Die("%u headers, %u CA index entries", N, SIZE(off_ca));
		}

	fsdb_ctx C;
	C.N = N;
	C.want_fasta = want_fasta;
	C.want_3di = want_3di;
	C.want_cal = want_cal;
	C.off_h = &off_h;
	C.len_h = &len_h;
	C.off_aa = &off_aa;
	C.len_aa = &len_aa;
	C.off_ss = &off_ss;
	C.len_ss = &len_ss;
	C.off_ca = &off_ca;
	C.len_ca = &len_ca;
	C.fH = OpenStdioFile(Prefix + "_h");
	if (want_fasta || want_cal)
		C.fAA = OpenStdioFile(Prefix);
	if (want_3di)
		C.fSS = OpenStdioFile(Prefix + "_ss");
	if (want_cal)
		{
		C.fCA = OpenStdioFile(Prefix + "_ca");
		C.ca_size = GetStdioFileSize64(C.fCA);
		asserta(N > 0);
		if (off_ca[N-1] + len_ca[N-1] != C.ca_size)
			Die("_ca file size %llu, LastOffset + LastLength = %llu",
				(unsigned long long) C.ca_size,
				(unsigned long long) (off_ca[N-1] + len_ca[N-1]));
		}
	C.done = 0;

	uint ThreadCount = GetRequestedThreadCount();
	if (ThreadCount < 1)
		ThreadCount = 1;
	if (ThreadCount > N)
		ThreadCount = N;

	vector<string> fasta_tmps;
	vector<string> di_tmps;
	vector<string> cal_tmps;
	C.fasta_fs.resize(ThreadCount);
	C.di_fs.resize(ThreadCount);
	C.cal_fs.resize(ThreadCount);

	const bool shard = (ThreadCount > 1);
	for (uint t = 0; t < ThreadCount; ++t)
		{
		if (want_fasta)
			{
			if (shard)
				{
				string fn;
				Ps(fn, "%s.tmp.%u", opt(fasta), t);
				fasta_tmps.push_back(fn);
				C.fasta_fs[t] = CreateStdioFile(fn);
				}
			else
				C.fasta_fs[t] = CreateStdioFile(opt(fasta));
			}
		if (want_3di)
			{
			if (shard)
				{
				string fn;
				Ps(fn, "%s.tmp.%u", opt(3di), t);
				di_tmps.push_back(fn);
				C.di_fs[t] = CreateStdioFile(fn);
				}
			else
				C.di_fs[t] = CreateStdioFile(opt(3di));
			}
		if (want_cal)
			{
			if (shard)
				{
				string fn;
				Ps(fn, "%s.tmp.%u", opt(cal), t);
				cal_tmps.push_back(fn);
				C.cal_fs[t] = CreateStdioFile(fn);
				}
			else
				C.cal_fs[t] = CreateStdioFile(opt(cal));
			}
		}

	vector<thread *> ts;
	const uint per = (N + ThreadCount - 1)/ThreadCount;
	for (uint t = 0; t < ThreadCount; ++t)
		{
		uint begin = t*per;
		uint end = begin + per;
		if (end > N)
			end = N;
		if (begin >= end)
			continue;
		ts.push_back(new thread(ConvertRange, &C, t, begin, end));
		}
	for (uint i = 0; i < SIZE(ts); ++i)
		{
		ts[i]->join();
		delete ts[i];
		}

	for (uint t = 0; t < ThreadCount; ++t)
		{
		CloseStdioFile(C.fasta_fs[t]);
		CloseStdioFile(C.di_fs[t]);
		CloseStdioFile(C.cal_fs[t]);
		}
	CloseStdioFile(C.fH);
	CloseStdioFile(C.fAA);
	CloseStdioFile(C.fSS);
	CloseStdioFile(C.fCA);

	if (shard)
		{
		if (want_fasta)
			{
			Progress("finalizing fasta... ");
			ConcatShardFiles(opt(fasta), fasta_tmps);
			Progress("done\n");
			}
		if (want_3di)
			{
			Progress("finalizing 3di... ");
			ConcatShardFiles(opt(3di), di_tmps);
			Progress("done\n");
			}
		if (want_cal)
			{
			Progress("finalizing cal... ");
			ConcatShardFiles(opt(cal), cal_tmps);
			Progress("done\n");
			}
		}

	ProgressLog("%u converted\n", N);
	}
