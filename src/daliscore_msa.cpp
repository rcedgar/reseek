#include "myutils.h"
#include "seqdb.h"
#include "chainreader2.h"
#include "alpha.h"
#include <set>

double GetDALIScore(const PDBChain& Q, const PDBChain& T,
	const vector<uint>& PosQs, const vector<uint>& PosTs);

double GetDALIZFromScoreAndLengths(double DALIScore, uint QL, uint TL)
	{
	double n12 = sqrt(QL * TL);
	double x = min(n12, 400.0);
	double mean = 7.9494 + 0.70852 * x + 2.5895e-4 * x * x - 1.9156e-6 * x * x * x;
	if (n12 > 400)
		mean += n12 - 400.0;
	double sigma = 0.5 * mean;
	double z = (DALIScore - mean) / max(1.0, sigma);
	return z;
	}

void GetAlignedPositions(const string& RowQ, const string& RowR,
	vector<uint>& PosQs, vector<uint>& PosRs, vector<bool>* ptrCore)
	{
	PosQs.clear();
	PosRs.clear();
	const uint ColCount = SIZE(RowQ);
	asserta(SIZE(RowR) == ColCount);
	uint PosQ = 0;
	uint PosR = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		if (ptrCore != 0 && !(*ptrCore)[Col])
			continue;
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

void GetUngappedSeq(const string& Row, string& Seq)
	{
	for (uint i = 0; i < SIZE(Row); ++i)
		{
		char c = Row[i];
		if (!isgap(c))
			Seq += toupper(c);
		}
	}

static void GetCore(const SeqDB& MSA, vector<bool>& IsCore)
	{
	IsCore.clear();
	const uint SeqCount = MSA.GetSeqCount();
	const uint ColCount = MSA.GetColCount();
	const uint MaxGaps = SeqCount / 10 + 1;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint n = MSA.GetGapCount(Col);
		uint nu = MSA.GetUpperCount(Col);
		uint nl = MSA.GetLowerCount(Col);
		if (nu != 0 && nl != 0)
			Die("Mixed case");
		IsCore.push_back(n <= MaxGaps && nl == 0);
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

	const uint SeqCount = MSA.GetSeqCount();
	vector<bool> IsCore;
	GetCore(MSA, IsCore);

	set<string> NotFound;
	uint PairCount = 0;
	double SumZ = 0;
	double SumZCore = 0;
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

			for (int iCore = 0; iCore <= 1; ++iCore)
				{
				if (fOut != 0 && iCore == 0)
					{
					if (SeqIdx1 == SeqIdx2)
						fprintf(fOut, "label_i=%s\t(self-align)",
							Label1.c_str());
					else
						fprintf(fOut, "label_i=%s\tlabel_j=%s",
							Label1.c_str(), Label2.c_str());
					}
				string Path;
				vector<uint> Pos1s, Pos2s;
				if (iCore == 0)
					{
					GetAlignedPositions(Row1, Row2, Pos1s, Pos2s, 0);
					double Score = GetDALIScore(Chain1, Chain2, Pos1s, Pos2s);
					double Z = GetDALIZFromScoreAndLengths(Score, L1, L2);
					SumZ += Z;
					if (fOut != 0)
						fprintf(fOut, "\tz-all=%.1f", Z);
					}
				else if (iCore == 1)
					{
					GetAlignedPositions(Row1, Row2, Pos1s, Pos2s, &IsCore);
					double Score = GetDALIScore(Chain1, Chain2, Pos1s, Pos2s);
					double Z = GetDALIZFromScoreAndLengths(Score, L1, L2);
					SumZCore += Z;
					if (fOut != 0)
						fprintf(fOut, "\tz-core=%.1f", Z);
					}
				else
					asserta(false);
				if (SeqIdx1 != SeqIdx2)
					{
					++PairCount;
					}
				}
			if (fOut != 0)
				fprintf(fOut, "\n");
			}
		}
	double MeanZ = (PairCount == 0 ? 0 : SumZ / PairCount);
	double MeanZCore = (PairCount == 0 ? 0 : SumZCore / PairCount);
	string Name;
	GetStemName(g_Arg1, Name);
	ProgressLog("Z  all %.1f, core %.1f %s\n", MeanZ, MeanZCore, Name.c_str());
	if (fOut != 0)
		{
		fprintf(fOut, "Z-all=%.1f\tZ-core=%.1f\tfile=%s\n",
			MeanZ, MeanZCore, Name.c_str());
		}
	CloseStdioFile(fOut);
	}
