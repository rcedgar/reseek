#include "myutils.h"
#include "seqdb.h"
#include "chainreader2.h"
#include "alpha.h"
#include <set>

double GetDALIScore(const PDBChain& Q, const PDBChain& T,
	const vector<uint>& PosQs, const vector<uint>& PosTs);

static void GetAlignedPositions(const string& RowQ, const string& RowR,
	vector<uint>& PosQs, vector<uint>& PosRs)
{
	PosQs.clear();
	PosRs.clear();
	const uint ColCount = SIZE(RowQ);
	asserta(SIZE(RowR) == ColCount);
	uint PosQ = 0;
	uint PosR = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
	{
		char q = RowQ[Col];
		char r = RowR[Col];
		bool gq = isgap(q);
		bool gr = isgap(r);
		if (gq && gr)
			continue;
		else if (!gq && !gr)
		{
			if (isupper(q) && isupper(r))
			{
				PosQs.push_back(PosQ);
				PosRs.push_back(PosR);
			}
			else
				asserta(islower(q) && islower(r));
			++PosQ;
			++PosR;
		}
		else if (!gq && gr)
			++PosQ;
		else if (gq && !gr)
			++PosR;
		else
			asserta(false);
	}
}

static void GetUngappedSeq(const string& Row, string& Seq)
{
	for (uint i = 0; i < SIZE(Row); ++i)
	{
		char c = Row[i];
		if (!isgap(c))
			Seq += toupper(c);
	}
}

void cmd_daliscore_msa()
{
	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);
	FILE* fOut = CreateStdioFile(opt_output);
	ChainReader2 CR;
	CR.Open(opt_input);
	vector<PDBChain*> Chains;
	map<string, uint> LabelToChainIdx;
	for (;;)
	{
		PDBChain* ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		uint Idx = SIZE(Chains);
		Chains.push_back(ptrChain);
		LabelToChainIdx[ptrChain->m_Label] = Idx;
	}

	set<string> NotFound;
	const uint SeqCount = MSA.GetSeqCount();
	uint PairCount = 0;
	double Sum = 0;
	for (uint SeqIdx1 = 0; SeqIdx1 < SeqCount; ++SeqIdx1)
	{
		const string& Label1 = string(MSA.GetLabel(SeqIdx1));
		map<string, uint>::const_iterator iter1 = LabelToChainIdx.find(Label1);
		if (iter1 == LabelToChainIdx.end())
		{
			NotFound.insert(Label1);
			continue;
		}
		uint ChainIdx1 = iter1->second;
		PDBChain& Chain1 = *Chains[ChainIdx1];
		uint L1 = Chain1.GetSeqLength();
		const string& Row1 = MSA.GetSeq(SeqIdx1);
		string U1;
		GetUngappedSeq(Row1, U1);
		asserta(U1 == Chain1.m_Seq);
		for (uint SeqIdx2 = SeqIdx1; SeqIdx2 < SeqCount; ++SeqIdx2)
		{
			const string& Label2 = string(MSA.GetLabel(SeqIdx2));
			map<string, uint>::const_iterator iter2 = LabelToChainIdx.find(Label2);
			if (iter2 == LabelToChainIdx.end())
			{
				NotFound.insert(Label2);
				continue;
			}
			uint ChainIdx2 = iter2->second;
			PDBChain& Chain2 = *Chains[ChainIdx2];
			uint L2 = Chain2.GetSeqLength();
			const string& Row2 = MSA.GetSeq(SeqIdx2);
			string U2;
			GetUngappedSeq(Row2, U2);
			asserta(U2 == Chain2.m_Seq);

			string Path;
			vector<uint> Pos1s, Pos2s;
			GetAlignedPositions(Row1, Row2, Pos1s, Pos2s);
			double Z = 10*GetDALIScore(Chain1, Chain2, Pos1s, Pos2s);
			if (fOut != 0)
				fprintf(fOut, "z\t%.1f\t%s\t%s\n", Z, Label1.c_str(), Label2.c_str());
			if (SeqIdx1 != SeqIdx2)
			{
				++PairCount;
				Sum += Z;
			}
		}
	}
	double MeanZ = (PairCount == 0 ? 0 : Sum / PairCount);
	string Name;
	GetStemName(g_Arg1, Name);
	ProgressLog("Z  %.1f  %s\n", MeanZ, Name.c_str());
	if (fOut != 0)
		fprintf(fOut, "Z\t%.1f\t%s\n", MeanZ, Name.c_str());
	CloseStdioFile(fOut);
}
