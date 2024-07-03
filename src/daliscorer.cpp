#include "myutils.h"
#include "seqdb.h"
#include "chainreader2.h"
#include "alpha.h"
#include <set>
#include "daliscorer.h"

static void GetUngappedSeq(const string& Row, string& Seq)
	{
	for (uint i = 0; i < SIZE(Row); ++i)
		{
		char c = Row[i];
		if (!isgap(c))
			Seq += toupper(c);
		}
	}

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
				if (ptrCore == 0 || (*ptrCore)[Col])
					{
					PosQs.push_back(PosQ);
					PosRs.push_back(PosR);
					}
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

void DALIScorer::LoadChains(const string &FN)
	{
	ChainReader2 CR;
	CR.Open(FN);

	for (;;)
		{
		PDBChain* ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		uint Idx = SIZE(m_Chains);
		m_Chains.push_back(ptrChain);
		m_SeqToChainIdx[ptrChain->m_Seq] = Idx;
		}
	}

void DALIScorer::SetCore()
	{
	const SeqDB &MSA = *m_MSA;
	m_ColIsCore.clear();
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
		m_ColIsCore.push_back(n <= MaxGaps && nl == 0);
		}
	}

void DALIScorer::SetSeqIdxToChainIdx()
	{
	m_SeqIdxToChainIdx.clear();
	m_NotFoundLabels.clear();
	const uint SeqCount = m_MSA->GetSeqCount();
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		const string &Label = string(m_MSA->GetLabel(SeqIdx));
		string Seq;
		m_MSA->GetSeq_StripGaps(SeqIdx, Seq, true);
		map<string, uint>::const_iterator iter = m_SeqToChainIdx.find(Seq);
		if (iter == m_SeqToChainIdx.end())
			{
			m_NotFoundLabels.insert(Label);
			m_SeqIdxToChainIdx.push_back(UINT_MAX);
			}
		else
			m_SeqIdxToChainIdx.push_back(iter->second);
		}
	}

void DALIScorer::Run(const string &Name, const SeqDB &MSA)
	{
	ClearMSA();

	m_Name = Name;
	m_MSA = &MSA;

	SetCore();
	SetSeqIdxToChainIdx();

	const uint SeqCount = MSA.GetSeqCount();
	double SumScore = 0;
	double SumScore_core = 0;
	double SumZ = 0;
	double SumZ_core = 0;
	uint PairCount = 0;
	for (uint SeqIdx1 = 0; SeqIdx1 < SeqCount; ++SeqIdx1)
		{
		uint ChainIdx1 = m_SeqIdxToChainIdx[SeqIdx1];
		if (ChainIdx1 == UINT_MAX)
			continue;
		asserta(ChainIdx1 < SIZE(m_Chains));
		const PDBChain &Chain1 = *m_Chains[ChainIdx1];
		const string &Row1 = MSA.GetSeq(SeqIdx1);
		const uint L1 = Chain1.GetSeqLength();
		const string &Label1 = Chain1.m_Label;
		string U1;
		GetUngappedSeq(Row1, U1);
		asserta(U1 == Chain1.m_Seq);
		for (uint SeqIdx2 = SeqIdx1; SeqIdx2 < SeqCount; ++SeqIdx2)
			{
			uint ChainIdx2 = m_SeqIdxToChainIdx[SeqIdx2];
			if (ChainIdx2 == UINT_MAX)
				continue;
			asserta(ChainIdx2 < SIZE(m_Chains));
			const PDBChain &Chain2 = *m_Chains[ChainIdx2];
			const string &Row2 = MSA.GetSeq(SeqIdx2);
			const uint L2 = Chain2.GetSeqLength();
			const string &Label2 = Chain2.m_Label;
			string U2;
			GetUngappedSeq(Row2, U2);
			asserta(U2 == Chain2.m_Seq);

			vector<uint> Pos1s;
			vector<uint> Pos2s;
			GetAlignedPositions(Row1, Row2, Pos1s, Pos2s, 0);
			double Score = GetDALIScore(Chain1, Chain2, Pos1s, Pos2s);
			double Z = GetDALIZFromScoreAndLengths(Score, L1, L2);

			vector<uint> Pos1s_core;
			vector<uint> Pos2s_core;
			GetAlignedPositions(Row1, Row2, Pos1s_core, Pos2s_core, &m_ColIsCore);
			double Score_core = GetDALIScore(Chain1, Chain2, Pos1s_core, Pos2s_core);
			double Z_core = GetDALIZFromScoreAndLengths(Score_core, L1, L2);

			m_MSALabels1.push_back(Label1);
			m_MSALabels2.push_back(Label2);

			m_Scores.push_back(Score);
			m_Scores_core.push_back(Score_core);

			m_Zs.push_back(Z);
			m_Zs_core.push_back(Z_core);

			if (SeqIdx1 != SeqIdx2)
				{
				++PairCount;
				SumScore += Score;
				SumScore_core += Score_core;
				SumZ += Z;
				SumZ_core += Z_core;
				}
			}
		}
	if (PairCount == 0)
		{
		m_MeanScore = DBL_MAX;
		m_MeanScore_core = DBL_MAX;
		m_MeanZ = DBL_MAX;
		m_MeanZ_core = DBL_MAX;
		}
	else
		{
		m_MeanScore = SumScore/PairCount;
		m_MeanScore_core = SumScore_core/PairCount;
		m_MeanZ = SumZ/PairCount;
		m_MeanZ_core = SumZ_core/PairCount;;
		}
	}

void DALIScorer::ToFev(FILE *f) const
	{
	if (f == 0)
		return;
	const uint PairCount = SIZE(m_MSALabels1);
	uint NotFound = SIZE(m_NotFoundLabels);
	fprintf(f, "set=%s\tpairs=%u\tnotfound=%u",
	  m_Name.c_str(), PairCount, NotFound);
	for (set<string>::const_iterator iter = m_NotFoundLabels.begin();
	  iter != m_NotFoundLabels.end(); ++iter)
		fprintf(f, "\t%s", iter->c_str());
	fprintf(f, "\n");
	for (uint i = 0; i < PairCount; ++i)
		{
		const string &Label1 = m_MSALabels1[i];
		const string &Label2 = m_MSALabels2[i];
		fprintf(f, "set=%s", m_Name.c_str());
		fprintf(f, "\tpair=%u", i);
		fprintf(f, "\tlabel1=%s", Label1.c_str());
		fprintf(f, "\tlabel1=%s", Label2.c_str());
		fprintf(f, "\tself=%s", (Label1 == Label2 ? "yes" : "no"));
		fprintf(f, "\tscore=%.1f", m_Scores[i]);
		fprintf(f, "\tscore_core=%.1f", m_Scores_core[i]);
		fprintf(f, "\tZ=%.2f", m_Zs[i]);
		fprintf(f, "\tZ_core=%.2f", m_Zs_core[i]);
		fprintf(f, "\n");
		}
	fprintf(f, "set=%s", m_Name.c_str());
	fprintf(f, "\tavg_score=%.1f", m_MeanScore);
	fprintf(f, "\tavg_score_core=%.1f", m_MeanScore_core);
	fprintf(f, "\tavg_Z_core=%.2f", m_MeanZ_core);
	fprintf(f, "\tavg_Z=%.2f", m_MeanZ);
	fprintf(f, "\n");
	}
