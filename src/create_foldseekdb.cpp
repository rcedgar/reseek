#include "myutils.h"
#include "chainreader2.h"
#include "seqdb.h"

char *CoordsToMem(const float *coords, int L, uint &membytes);
float *GetCoordsFromMem(const char *mem, uint chainLength, uint entryLength);

static void Write_dbtype(const string &Prefix,
						 const string &Suffix,
						 uint32_t Value)
	{
	FILE *f = CreateStdioFile(Prefix + Suffix + ".dbtype");
	WriteStdioFile(f, &Value, 4);
	CloseStdioFile(f);
	}

void cmd_create_foldseekdb()
	{
	uint MinChainLength = 1;
	if (optset_minchainlength)
		MinChainLength = opt_minchainlength;

	asserta(optset_output);
	asserta(optset_3di);

	const string Prefix = string(opt_output);

	SeqDB SS;
	SS.FromFasta(opt_3di);
	const uint SeqCount = SS.GetSeqCount();
	for (uint i = 0; i < SeqCount; ++i)
		TruncateAtFirstWhiteSpace(SS.m_Labels[i]);
	SS.SetLabelToIndex();

	Write_dbtype(Prefix, "", 0);
	Write_dbtype(Prefix, "_h", 0xc);
	Write_dbtype(Prefix, "_ca", 0x65);
	Write_dbtype(Prefix, "_ss", 0);

	FILE *fSeqs = CreateStdioFile(Prefix);
	FILE *fLabels = CreateStdioFile(Prefix + "_h");
	FILE *fSource = CreateStdioFile(Prefix + ".source");
	FILE *fCA = CreateStdioFile(Prefix + "_ca");
	FILE *fSeqs3Di = CreateStdioFile(Prefix + "_ss");

	FILE *fSeqsLookup = CreateStdioFile(Prefix + ".lookup");
	FILE *fSeqsIndex = CreateStdioFile(Prefix + ".index");
	FILE *fSeqs3DiIndex = CreateStdioFile(Prefix + "_ss.index");
	FILE *fLabelsIndex = CreateStdioFile(Prefix + "_h.index");
	FILE *fCAIndex = CreateStdioFile(Prefix + "_ca.index");

	ChainReader2 CR;
	CR.Open(g_Arg1);

	uint32_t Idx = 0;
	uint64_t SeqOffset = 0;
	uint64_t LabelOffset = 0;
	uint64_t CAOffset = 0;
	const char nl_null[2] = { '\n', 0 };
	for (;;)
		{
		PDBChain *Chain = CR.GetNext();
		if (Chain == 0)
			break;

		string Label = Chain->m_Label.c_str();
		TruncateAtFirstWhiteSpace(Label);

		string Seq3Di;
		bool Found = SS.GetSeqByLabel(Label, Seq3Di, false);
		if (!Found)
			Die("Missing 3Di sequence >%s", Label.c_str());
		const uint SeqLen3Di = SIZE(Seq3Di);

		const string &Seq = Chain->m_Seq.c_str();

		const uint SeqLen = Chain->GetSeqLength();
		const uint LabelLen = SIZE(Label);

		if (SeqLen != SeqLen3Di)
			Die("Seqence length mismatch, aa=%u 3Di=%u >%s",
				SeqLen, SeqLen3Di, Label.c_str());

		WriteStdioFile(fLabels, Label.c_str(), LabelLen);
		WriteStdioFile(fLabels, nl_null, 2);

		WriteStdioFile(fSeqs, Seq.c_str(), SeqLen);
		WriteStdioFile(fSeqs, nl_null, 2);

		WriteStdioFile(fSeqs3Di, Seq3Di.c_str(), SeqLen);
		WriteStdioFile(fSeqs3Di, nl_null, 2);

		fprintf(fSeqsLookup, "%u\t%s\t%u\n",
				Idx, Label.c_str(), Idx);

		fprintf(fSource, "%u\t%s\n",
				Idx, Label.c_str());

		fprintf(fSeqsIndex, "%u\t%llu\t%u\n",
				Idx, SeqOffset, SeqLen+2);

		fprintf(fSeqs3DiIndex, "%u\t%llu\t%u\n",
				Idx, SeqOffset, SeqLen+2);

		fprintf(fLabelsIndex, "%u\t%llu\t%u\n",
				Idx, LabelOffset, LabelLen+2);

		SeqOffset += SeqLen + 2;
		LabelOffset += LabelLen + 2;

		float *Coords = myalloc(float, 3*SeqLen);

		for (uint i = 0; i < SeqLen; ++i)
			{
			float X = Chain->m_Xs[i];
			float Y = Chain->m_Ys[i];
			float Z = Chain->m_Zs[i];

			Coords[i] = X;
			Coords[SeqLen+i] = Y;
			Coords[2*SeqLen+i] = Z;
			}

		uint membytes;
		char *mem = CoordsToMem(Coords, SeqLen, membytes);
		if (mem == 0)
			{
			uint floatbytes = 3*SeqLen*sizeof(float);
			Log("Overflow %s\n", Label.c_str());
			fprintf(fCAIndex, "%u\t%llu\t%u\n",
					Idx, CAOffset, floatbytes+2);
			WriteStdioFile(fCA, Coords, floatbytes);
			WriteStdioFile(fCA, nl_null, 2);
			CAOffset += floatbytes + 2;
			}
		else
			{
			fprintf(fCAIndex, "%u\t%llu\t%u\n",
					Idx, CAOffset, membytes+2);
			WriteStdioFile(fCA, mem, membytes);
			WriteStdioFile(fCA, nl_null, 2);
			CAOffset += membytes + 2;
			}

		myfree(Coords);
		delete Chain;
		++Idx;
		}

	CloseStdioFile(fLabels);
	CloseStdioFile(fSeqs);
	CloseStdioFile(fSeqsLookup);
	CloseStdioFile(fSeqsIndex);
	CloseStdioFile(fSeqs3DiIndex);
	CloseStdioFile(fLabelsIndex);
	CloseStdioFile(fCAIndex);
	CloseStdioFile(fCA);
	CloseStdioFile(fSeqs3Di);
	CloseStdioFile(fSource);
	}
