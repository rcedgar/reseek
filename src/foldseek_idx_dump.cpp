#include "myutils.h"

// Dump Foldseek/MMseqs createindex .idx k-mer inverted index to TSV.
//
// Input: Foldseek createindex output <db>.idx (+ .index / .dbtype).
//   Arg may be <db> or <db>.idx, e.g. 100_ss or 100_ss.idx
//
// Output line (non-empty k-mers only):
//   <kmer_letters>\t<kmer_int>\t<n_postings>\t<db_idx>\t<pos>\t...
//
// META fields are written as leading '#...' comment lines.
//
// Record keys from MMseqs PrefilteringIndexReader:
//   VERSION=0 META=1 ENTRIES=9 ENTRIESOFFSETS=10 ENTRIESNUM=12 SEQCOUNT=13
//   per-split keys are offset by 1000*split.

#pragma pack(push, 1)
struct FoldseekIndexEntryLocal
	{
	uint32_t seqId;
	uint16_t position_j;
	};
#pragma pack(pop)

static const uint KEY_VERSION = 0;
static const uint KEY_META = 1;
static const uint KEY_ENTRIES = 9;
static const uint KEY_ENTRIESOFFSETS = 10;
static const uint KEY_ENTRIESNUM = 12;
static const uint KEY_SEQCOUNT = 13;

// MMseqs AA / 3Di num2aa without trailing X (alphabetSize-1 for indexing).
static const char *NUM2AA20 = "ARNDCQEGHILKMFPSTWYV";

static void ReadKeyedIndex(const string &FN,
	map<uint, uint64> &KeyToOffset,
	map<uint, uint64> &KeyToLength)
	{
	KeyToOffset.clear();
	KeyToLength.clear();
	FILE *f = OpenStdioFile(FN);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		uint Key = StrToUint(Fields[0]);
		uint64 Offset = StrToUint64(Fields[1]);
		uint64 Length = StrToUint64(Fields[2]);
		asserta(Length > 0);
		asserta(KeyToOffset.find(Key) == KeyToOffset.end());
		KeyToOffset[Key] = Offset;
		KeyToLength[Key] = Length;
		}
	CloseStdioFile(f);
	asserta(SIZE(KeyToOffset) > 0);
	}

static byte *ReadKeyBytes(FILE *fData,
	const map<uint, uint64> &KeyToOffset,
	const map<uint, uint64> &KeyToLength,
	uint Key, uint64 &OutLength)
	{
	map<uint, uint64>::const_iterator itOff = KeyToOffset.find(Key);
	map<uint, uint64>::const_iterator itLen = KeyToLength.find(Key);
	asserta(itOff != KeyToOffset.end());
	asserta(itLen != KeyToLength.end());
	uint64 Offset = itOff->second;
	uint64 Length = itLen->second;
	asserta(Length > 0);
	asserta(Length < (uint64(1) << 40));
	byte *Buf = myalloc64(byte, Length);
	ReadStdioFile64(fData, Offset, Buf, Length);
	OutLength = Length;
	return Buf;
	}

static uint64 IPow(uint64 Base, uint Exp)
	{
	uint64 R = 1;
	for (uint i = 0; i < Exp; ++i)
		{
		asserta(R <= UINT64_MAX / Base);
		R *= Base;
		}
	return R;
	}

static void KmerIntToLetters(uint64 Kmer, int KmerSize, uint AlphabetSize,
	string &Out)
	{
	asserta(KmerSize > 0);
	asserta(KmerSize <= 32);
	asserta(AlphabetSize >= 2);
	Out.clear();
	Out.resize((size_t) KmerSize);
	uint64 Idx = Kmer;
	// Same order as Indexer::index2int (least-significant letter first in powers[]).
	vector<uint64> Powers((size_t) KmerSize);
	uint64 Pow = 1;
	for (int i = 0; i < KmerSize; ++i)
		{
		Powers[(size_t) i] = Pow;
		asserta(Pow <= UINT64_MAX / AlphabetSize);
		Pow *= AlphabetSize;
		}
	for (int i = KmerSize - 1; i >= 0; --i)
		{
		uint64 Digit = Idx / Powers[(size_t) i];
		Idx = Idx - Digit * Powers[(size_t) i];
		asserta(Digit < AlphabetSize);
		asserta(AlphabetSize == 20);
		asserta(Digit < 20);
		Out[(size_t) i] = NUM2AA20[Digit];
		}
	asserta(Idx == 0);
	}

void cmd_foldseek_idx_dump()
	{
	const string &Arg = g_Arg1;
	asserta(Arg.size() > 0);

	// Foldseek createindex writes <db>.idx alongside <db>.
	// Accept either ".../100_ss" or ".../100_ss.idx".
	string IdxPrefix = Arg;
	const size_t n = IdxPrefix.size();
	if (!(n >= 4 && IdxPrefix.substr(n - 4) == ".idx"))
		IdxPrefix += ".idx";

	const string DataFN = IdxPrefix;
	const string IndexFN = IdxPrefix + ".index";
	const string DbtypeFN = IdxPrefix + ".dbtype";

	FILE *fout = CreateStdioFile(opt(output));

	FILE *fDbtype = OpenStdioFile(DbtypeFN);
	asserta(GetStdioFileSize64(fDbtype) == 4);
	uint32_t Dbtype = UINT32_MAX;
	ReadStdioFile(fDbtype, &Dbtype, 4);
	CloseStdioFile(fDbtype);
	// Parameters::DBTYPE_INDEX_DB = 9
	asserta(Dbtype == 9);

	map<uint, uint64> KeyToOffset;
	map<uint, uint64> KeyToLength;
	ReadKeyedIndex(IndexFN, KeyToOffset, KeyToLength);

	FILE *fData = OpenStdioFile(DataFN);
	const uint64 DataSize = GetStdioFileSize64(fData);
	asserta(DataSize > 0);

	uint64 MetaLen = 0;
	byte *MetaBytes = ReadKeyBytes(fData, KeyToOffset, KeyToLength, KEY_META, MetaLen);
	asserta(MetaLen >= 12 * sizeof(int));
	const int *Meta = (const int *) MetaBytes;
	const int MaxSeqLength = Meta[0];
	const int KmerSize = Meta[1];
	const int CompBiasCorr = Meta[2];
	const int AlphabetSizeMeta = Meta[3];
	const int Masked = Meta[4];
	const int Spaced = Meta[5];
	const int KmerThr = Meta[6];
	const int SeqType = Meta[7];
	const int SrcSeqType = Meta[8];
	const int Headers1 = Meta[9];
	const int Headers2 = Meta[10];
	int Splits = Meta[11];
	if (Splits == 0)
		Splits = 1;
	asserta(KmerSize >= 1 && KmerSize <= 32);
	asserta(AlphabetSizeMeta >= 2 && AlphabetSizeMeta <= 32);
	asserta(Splits >= 1 && Splits <= 1024);
	asserta(MaxSeqLength > 0);

	// AA / nucleotide indexes drop X (or N) from the k-mer alphabet.
	// Foldseek 3Di uses DBTYPE_AMINO_ACIDS (0).
	const int DBTYPE_AMINO_ACIDS = 0;
	const int DBTYPE_NUCLEOTIDES = 1;
	int AdjustAlphabetSize = AlphabetSizeMeta;
	if (SeqType == DBTYPE_AMINO_ACIDS || SeqType == DBTYPE_NUCLEOTIDES)
		AdjustAlphabetSize = AlphabetSizeMeta - 1;
	asserta(AdjustAlphabetSize >= 2);
	asserta(AdjustAlphabetSize == 20);

	const uint64 TableSize = IPow((uint64) AdjustAlphabetSize, (uint) KmerSize);
	asserta(TableSize > 0);

	ProgressLog("# foldseek_idx_dump\n");
	ProgressLog("# prefix\t%s\n", IdxPrefix.c_str());
	ProgressLog("# dbtype\t%u\n", Dbtype);
	ProgressLog("# MaxSeqLength\t%d\n", MaxSeqLength);
	ProgressLog("# KmerSize\t%d\n", KmerSize);
	ProgressLog("# CompBiasCorr\t%d\n", CompBiasCorr);
	ProgressLog("# AlphabetSize\t%d\n", AlphabetSizeMeta);
	ProgressLog("# AdjustAlphabetSize\t%d\n", AdjustAlphabetSize);
	ProgressLog("# TableSize\t%llu\n", (unsigned long long) TableSize);
	ProgressLog("# Masked\t%d\n", Masked);
	ProgressLog("# Spaced\t%d\n", Spaced);
	ProgressLog("# KmerThr\t%d\n", KmerThr);
	ProgressLog("# SequenceType\t%d\n", SeqType);
	ProgressLog("# SourceSeqType\t%d\n", SrcSeqType);
	ProgressLog("# Headers1\t%d\n", Headers1);
	ProgressLog("# Headers2\t%d\n", Headers2);
	ProgressLog("# Splits\t%d\n", Splits);
	ProgressLog("# num2aa\t%s\n", NUM2AA20);
	ProgressLog("# fields\tkmer_letters\tkmer_int\tn_postings\t(db_idx\tpos)*\n");

	if (KeyToOffset.find(KEY_VERSION) != KeyToOffset.end())
		{
		uint64 VerLen = 0;
		byte *VerBytes = ReadKeyBytes(fData, KeyToOffset, KeyToLength, KEY_VERSION, VerLen);
		asserta(VerLen >= 1);
		string Version;
		for (uint64 i = 0; i < VerLen; ++i)
			{
			char c = (char) VerBytes[i];
			if (c == 0)
				break;
			Version.push_back(c);
			}
		ProgressLog("# IndexVersion\t%s\n", Version.c_str());
		myfree(VerBytes);
		}

	ProgressLog("META k=%d alphabet=%d (adj %d) tableSize=%llu splits=%d\n",
		KmerSize, AlphabetSizeMeta, AdjustAlphabetSize,
		(unsigned long long) TableSize, Splits);
	if (fout == 0) return;

	uint64 TotalNonEmpty = 0;
	uint64 TotalPostingsWritten = 0;

	for (int Split = 0; Split < Splits; ++Split)
		{
		const uint SplitOffset = (uint) Split * 1000;
		const uint KeyEntries = SplitOffset + KEY_ENTRIES;
		const uint KeyOffsets = SplitOffset + KEY_ENTRIESOFFSETS;
		const uint KeyEntriesNum = SplitOffset + KEY_ENTRIESNUM;
		const uint KeySeqCount = SplitOffset + KEY_SEQCOUNT;

		uint64 EntriesNumLen = 0;
		byte *EntriesNumBytes = ReadKeyBytes(fData, KeyToOffset, KeyToLength,
			KeyEntriesNum, EntriesNumLen);
		asserta(EntriesNumLen >= sizeof(uint64_t));
		const uint64 EntriesNum = *(const uint64_t *) EntriesNumBytes;
		myfree(EntriesNumBytes);

		uint64 SeqCountLen = 0;
		byte *SeqCountBytes = ReadKeyBytes(fData, KeyToOffset, KeyToLength,
			KeySeqCount, SeqCountLen);
		asserta(SeqCountLen >= sizeof(uint64_t));
		const uint64 SeqCount = *(const uint64_t *) SeqCountBytes;
		myfree(SeqCountBytes);

		Log("# split\t%d\tentriesNum\t%llu\tseqCount\t%llu\n",
			Split,
			(unsigned long long) EntriesNum,
			(unsigned long long) SeqCount);

		uint64 OffsetsBytesLen = 0;
		byte *OffsetsBytes = ReadKeyBytes(fData, KeyToOffset, KeyToLength,
			KeyOffsets, OffsetsBytesLen);
		const uint64 NeedOffsetsBytes = (TableSize + 1) * sizeof(uint64_t);
		asserta(OffsetsBytesLen >= NeedOffsetsBytes);
		const uint64_t *Offsets = (const uint64_t *) OffsetsBytes;
		asserta(Offsets[0] == 0);
		asserta(Offsets[TableSize] == EntriesNum);

		uint64 EntriesBytesLen = 0;
		byte *EntriesBytes = ReadKeyBytes(fData, KeyToOffset, KeyToLength,
			KeyEntries, EntriesBytesLen);
		const uint64 NeedEntriesBytes = EntriesNum * sizeof(FoldseekIndexEntryLocal);
		asserta(sizeof(FoldseekIndexEntryLocal) == 6);
		asserta(EntriesBytesLen >= NeedEntriesBytes);
		const FoldseekIndexEntryLocal *Entries =
			(const FoldseekIndexEntryLocal *) EntriesBytes;

		ProgressLog("split %d/%d entriesNum=%llu seqCount=%llu\n",
			Split + 1, Splits,
			(unsigned long long) EntriesNum,
			(unsigned long long) SeqCount);

		string KmerLetters;
		for (uint64 Kmer = 0; Kmer < TableSize; ++Kmer)
			{
			const uint64 Start = Offsets[Kmer];
			const uint64 End = Offsets[Kmer + 1];
			asserta(End >= Start);
			asserta(End <= EntriesNum);
			const uint64 N = End - Start;
			if (N == 0)
				continue;

			KmerIntToLetters(Kmer, KmerSize, (uint) AdjustAlphabetSize, KmerLetters);
			fprintf(fout, "%s\t%llu\t%llu",
				KmerLetters.c_str(),
				(unsigned long long) Kmer,
				(unsigned long long) N);
			for (uint64 i = 0; i < N; ++i)
				{
				const FoldseekIndexEntryLocal &E = Entries[Start + i];
				fprintf(fout, "\t%u\t%u", E.seqId, (uint) E.position_j);
				}
			fprintf(fout, "\n");
			++TotalNonEmpty;
			TotalPostingsWritten += N;
			}

		myfree(OffsetsBytes);
		myfree(EntriesBytes);
		}

	myfree(MetaBytes);
	CloseStdioFile(fData);
	CloseStdioFile(fout);

	ProgressLog("%llu non-empty k-mers, %llu postings written\n",
		(unsigned long long) TotalNonEmpty,
		(unsigned long long) TotalPostingsWritten);
	}
