#include "myutils.h"

void LogCoords16(const char *mem, uint chainLength);
float *GetCoordsFromMem(const char *mem, uint chainLength, uint entryLength);

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

static void ReadNulTerminatedSeqs(const string &FN,
								  vector<string> &Seqs)
	{
	Seqs.clear();
	FILE *f = OpenStdioFile(FN);
	if (GetStdioFileSize64(f) == 0)
		Die("Empty file: %s", FN.c_str());
	string Seq;
	for (;;)
		{
		int c = fgetc(f);
		if (c == 0)
			{
			Seqs.push_back(Seq);
			Seq.clear();
			}
		else if (c == '\n' || c == '\r')
			continue;
		else
			Seq += c;
		if (c < 0)
			break;
		}
	CloseStdioFile(f);
	}

static void VerifyLookup(const string &Prefix,
						 const string &Suffix,
						 const vector<string> &Labels)
	{
	bool IsLookup = false;
	if (Suffix == ".lookup")
		IsLookup = true;
	else if (Suffix == ".source")
		IsLookup = false;
	else
		asserta(false);
	const string FN(Prefix + Suffix);
	string Line;
	vector<string> Fields;
	FILE *f = OpenStdioFile(FN);
	uint Idx = 0;
	uint ErrCount = 0;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		uint n = SIZE(Fields);
		if (IsLookup)
			{
			if (n != 3)
				Die("Expected 3 fields, got %u %s:%u",
					n, FN.c_str(), Idx + 1);
			}
		else
			{
			if (n != 2)
				Die("Expected 2 fields, got %u %s:%u",
					n, FN.c_str(), Idx + 1);
			}
		if (Idx >= SIZE(Labels))
			Die("More labels than expected in %s", FN.c_str());
		const string &Label = Fields[1];

		vector<string> Fields2;
		const string &Label_Idx = Labels[Idx];
		Split(Label_Idx, Fields2, 0);
		const string &ExpectedLabel = Fields2[0];

		if (Label != ExpectedLabel)
			{
			Log("Idx %u Label='%s' ExpectedLabel='%s' Label_Idx='%s'\n",
				Idx, Label.c_str(), ExpectedLabel.c_str(), Label_Idx.c_str());
			++ErrCount;
			}
		++Idx;
		}
	if (ErrCount > 0)
		Warning("%u errs checking lookup file '%s', see log", ErrCount, FN.c_str());
	CloseStdioFile(f);
	}

static void ReadIndex(const string &FN,
					  vector<uint> &Offsets,
					  vector<uint> &Lengths)
	{
	Offsets.clear();
	Lengths.clear();
	FILE *f = OpenStdioFile(FN);
	string Line;
	vector<string> Fields;
	uint ExpectedIdx = 0;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 3)
			Die("Expected 3 fields, got '%s' in %s",
				Line.c_str(), FN.c_str());
		uint Idx = StrToUint(Fields[0]);
		uint Offset = StrToUint(Fields[1]);
		uint Length = StrToUint(Fields[2]);
		if (Idx != ExpectedIdx)
			Die("Expected idx %u, got '%s' in %s",
				ExpectedIdx, Line.c_str(), FN.c_str());

		Offsets.push_back(Offset);
		Lengths.push_back(Length);

		++ExpectedIdx;
		}
	}

void cmd_convert_foldseekdb()
	{
	const string &Prefix = g_Arg1;

	CheckDBType(Prefix, ".dbtype", 0x0);
	CheckDBType(Prefix, "_h.dbtype", 0xC);
	CheckDBType(Prefix, "_ss.dbtype", 0x0);
	CheckDBType(Prefix, "_ca.dbtype", 0x65);

	vector<string> Labels;
	vector<string> SeqsAA;
	vector<string> Seqs3Di;
	Progress("Read labels\n");
	string FN = Prefix + "_h";
	ReadNulTerminatedSeqs(FN, Labels);
	const uint SeqCount = SIZE(Labels);
	Log("%u labels in '%s'\n", SeqCount, FN.c_str());

	Progress("Read aa seqs\n");
	FN = Prefix;
	ReadNulTerminatedSeqs(Prefix, SeqsAA);
	const uint SeqCountAA = SIZE(SeqsAA);
	Log("%u labels, %u aa seqs in '%s'\n", SeqCount, SeqCountAA, FN.c_str());
	if (SeqCountAA != SeqCount) Die("%u labels, %u aa seqs", SeqCount, SeqCountAA);

	Progress("Read 3Di seqs\n");
	FN = Prefix + "_ss";
	ReadNulTerminatedSeqs(FN, Seqs3Di);
	const uint SeqCount3Di = SIZE(Seqs3Di);
	Log("%u labels, %u 3Di seqs in '%s'\n", SeqCount, SeqCount3Di, FN.c_str());
	if (SeqCountAA != SeqCount) Die("%u labels, %u 3Di seqs", SeqCount, SeqCount3Di);

	Progress("Check lookup\n");
	VerifyLookup(Prefix, ".lookup", Labels);

	Progress("Check source\n");
	VerifyLookup(Prefix, ".source", Labels);

	Progress("Read CA coords\n");
	const string CAFN = Prefix + "_ca";
	uint64 Size64;
	byte *Data = ReadAllStdioFile64(CAFN, Size64);
	uint32_t Size32 = uint32_t(Size64);
	if (uint64_t(Size32) != Size64)
		Die("_ca file too big");

	vector<uint> Offsets_coords;
	vector<uint> Lengths_coords;
	Progress("Read CA index\n");
	ReadIndex(Prefix + "_ca.index", Offsets_coords, Lengths_coords);
	const uint SeqCount2 = SIZE(Offsets_coords);
	asserta(SIZE(Lengths_coords) == SeqCount2);
	if (SeqCount2 != SeqCount)
		Die("%u seqs in FASTA, %u seqs in ca.index",
			SeqCount, SeqCount2);

	asserta(SeqCount > 0);
	uint LastOffset = Offsets_coords[SeqCount-1];
	uint LastLength = Lengths_coords[SeqCount-1];
	if (LastOffset + LastLength != Size32)
		Die("_ca file size %u, LastOffset + LastLength = %u",
			Size32, LastOffset + LastLength);

	FILE *faa = CreateStdioFile(opt(fasta));
	FILE *f3Di = CreateStdioFile(opt(3di));
	FILE *fcal = CreateStdioFile(opt(cal));
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, SeqCount, "Write output files");
		const string &Label = Labels[SeqIdx];
		const string &SeqAA = SeqsAA[SeqIdx];
		const string &Seq3Di = Seqs3Di[SeqIdx];

		const uint Laa = SIZE(SeqAA);
		const uint L3Di = SIZE(Seq3Di);
		if (Laa != L3Di)
			Die("aa/3Di sequence mismatch %u, %u >%s",
				Laa, L3Di, Label.c_str());

		SeqToFasta(faa, Label, SeqAA);
		SeqToFasta(f3Di, Label, Seq3Di);

		if (fcal != 0)
			{
			uint CoordsOffset = Offsets_coords[SeqIdx];
			uint CoordsLength = Lengths_coords[SeqIdx];
			const char *mem = (const char *) (Data + CoordsOffset);
			float *Coords = GetCoordsFromMem(mem, Laa, CoordsLength);
		////////////////////////////////////////////
		//	char *mem2 = CoordsToMem(Coords, Laa);
		//	if (mem2 == 0)
		//		Warning("Overflow");
		//	else
		//		{
		//		float *Coords2 = GetCoordsFromMem(mem, Laa, CoordsLength);
		//		for (uint i = 0; i < Laa; ++i)
		//			{
		//			double dx = fabs(Coords[i] - Coords2[i]);
		//			double dy = fabs(Coords[Laa+i] - Coords2[Laa+i]);
		//			double dz = fabs(Coords[2*Laa+i] - Coords2[2*Laa+i]);
		//			asserta(dx < 0.2);
		//			asserta(dy < 0.2);
		//			asserta(dz < 0.2);
		//			}
		//		}
		////////////////////////////////////////////
			fprintf(fcal, ">%s\n", Label.c_str());
			for (uint i = 0; i < Laa; ++i)
				fprintf(fcal, "%c	%.1f\t%.1f\t%.1f\n",
						SeqAA[i], Coords[i], Coords[Laa+i], Coords[2*Laa+i]);
			if ((void *) Coords != (void *) mem)
				myfree((void *) Coords);
			}
		}
	CloseStdioFile(faa);
	CloseStdioFile(f3Di);
	CloseStdioFile(fcal);
	}
