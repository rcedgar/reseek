#if 0
#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"
#include "seqinfo.h"

void cmd_core_blocks()
	{
	const string &InputFN = opt(input);
	const string &OutputFN = opt(output);
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	const uint ChainCount = SIZE(Chains);
	map<string, uint> LabelToChainIdx;
	vector<string> Fields;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		string Label = Chains[ChainIdx]->m_Label;
		Split(Label, Fields, '/');
		LabelToChainIdx[Fields[0]] = ChainIdx;
		}

	FEATURE F = StrToFeature(opt(feature));
	DSS D;

/***
core_block      9       47      a.1.1
ARAVSIMKA       124     d1xg0c_
GAFQELLLS       125     d1or4a_
IAAFTFTRD       127     d1b8da_
IEALKYIKA       129     d3l0fa_
***/
	FILE *fIn = OpenStdioFile(InputFN);
	FILE *fOut = CreateStdioFile(OutputFN);
	string Line;
	set<string> NotFound;
	uint DiffSeq = 0;
	for (;;)
		{
		bool Ok = ReadLineStdioFile(fIn, Line);
		if (!Ok)
			break;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 4);
		asserta(Fields[0] == "core_block");
		uint Length = StrToUint(Fields[1]);
		uint NrSeqs = StrToUint(Fields[2]);
		const string &Name = Fields[3];
		vector<string> Seqs;
		vector<uint> PosVec;
		vector<string> Labels;
		for (uint i = 0; i < NrSeqs; ++i)
			{
			Ok = ReadLineStdioFile(fIn, Line);
			asserta(Ok);
			Split(Line, Fields, '\t');
			asserta(SIZE(Fields) == 3);
			const string &Seq = Fields[0];
			const uint Pos = StrToUint(Fields[1]);
			const string &Label = Fields[2];
			map<string, uint>::const_iterator iter = LabelToChainIdx.find(Label);
			if (iter == LabelToChainIdx.end())
				{
				NotFound.insert(Label);
				continue;
				}
			const PDBChain &Chain = *Chains[iter->second];
			const string &ChainSeq = Chain.m_Seq;
			D.Init(Chain);
			string ChainRowAA;
			for (uint k = Pos; k < Pos + Length; ++k)
				ChainRowAA += ChainSeq[k];
			if (ChainRowAA != Seq)
				{
				++DiffSeq;
				continue;
				}

			string ChainRow;
			for (uint k = Pos; k < Pos + Length; ++k)
				ChainRow += 'A' + D.GetFeature(F, k);

			Labels.push_back(Label);
			Seqs.push_back(ChainRow);
			PosVec.push_back(Pos);
			}
		const uint N = SIZE(Labels);
		asserta(SIZE(Seqs) == N);
		asserta(SIZE(PosVec) == N);
		if (N < 8)
			continue;
		fprintf(fOut, "core_block\t%u\t%u\t%s\n", Length, N, Name.c_str());
		for (uint i = 0; i < N; ++i)
			fprintf(fOut, "%s\t%u\t%s\n", Seqs[i].c_str(), PosVec[i], Labels[i].c_str());
		}
	CloseStdioFile(fIn);
	CloseStdioFile(fOut);
	ProgressLog("%u not found, %u diff seq\n", SIZE(NotFound), DiffSeq);
	}
#endif 