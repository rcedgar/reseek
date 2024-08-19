#include "myutils.h"
#include "seqdb.h"
#include "chainreader2.h"
#include "alpha.h"
#include "daliscorer.h"
#include "sort.h"
#include <set>

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
  vector<uint>& PosQs, vector<uint>& PosRs, const vector<bool>* ptrCore)
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
				{
				if (!opt_mixedcasecolok)
					asserta(islower(q) && islower(r));
				}
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

	time_t t0 = time(0);
	Progress("Reading chains\r");
	for (;;)
		{
		time_t t = time(0);
		if (t > t0)
			{
			Progress("Reading chains (%u)\r", SIZE(m_Chains));
			t0 = t;
			}
		PDBChain* ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		uint Idx = SIZE(m_Chains);
		m_Chains.push_back(ptrChain);
		m_SeqToChainIdx[ptrChain->m_Seq] = Idx;
		}
	Progress("Reading chains done (%u)\n", SIZE(m_Chains));
	{//@@
	for (map<string, uint>::const_iterator iter = m_SeqToChainIdx.begin();
		iter != m_SeqToChainIdx.end(); ++iter)
			{
			Log("[%u] %s\n", iter->second, iter->first);
			}

	}//@@
	}

void DALIScorer::SetCore()
	{
	m_CoreColCount = 0;
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
		bool IsCore =  n <= MaxGaps && nl == 0;
		if (IsCore)
			++m_CoreColCount;
		m_ColIsCore.push_back(IsCore);
		}
	}

bool DALIScorer::SetSeqIdxToChainIdx(bool MissingSeqOk)
	{
	m_SeqIdxToChainIdx.clear();
	m_NotFoundLabels.clear();
	const uint SeqCount = m_MSA->GetSeqCount();
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		const string &Label = string(m_MSA->GetLabel(SeqIdx));
		string Seq;
		m_MSA->GetSeq_StripGaps(SeqIdx, Seq, true);
		Log(">%s\n", Label.c_str());//@@
		Log("%s\n", Seq.c_str());//@@
		map<string, uint>::const_iterator iter = m_SeqToChainIdx.find(Seq);
		if (iter == m_SeqToChainIdx.end())
			{
			if (!MissingSeqOk)
				{
				Warning("Sequence not matched >%s", Label.c_str());
				return false;
				}
			m_NotFoundLabels.insert(Label);
			m_SeqIdxToChainIdx.push_back(UINT_MAX);
			}
		else
			{
			uint ChainIdx = iter->second;
			m_SeqIdxToChainIdx.push_back(ChainIdx);
			}
		}
	return true;
	}

bool DALIScorer::SetMSA(const string &Name, const SeqDB &MSA,
  bool DoCore, bool MissingSeqOk)
	{
	ClearMSA();
	m_Name = Name;
	m_MSA = &MSA;
	bool Ok = SetSeqIdxToChainIdx(MissingSeqOk);
	if (!Ok)
		{
		ClearMSA();
		return false;
		}

	m_DoCore = DoCore;
	if (m_DoCore)
		SetCore();
	else
		m_ColIsCore.clear();
	SetColToPosVec(DoCore);
	SetDistMxs();
	return true;
	}

bool DALIScorer::GetDALIRowPair(uint SeqIdx1, uint SeqIdx2,
  double &Score, double &Z) const
	{
	const SeqDB &MSA = *m_MSA;
	Score = DBL_MAX;
	Z = DBL_MAX;
	uint ChainIdx1 = m_SeqIdxToChainIdx[SeqIdx1];
	if (ChainIdx1 == UINT_MAX)
		return false;
	uint ChainIdx2 = m_SeqIdxToChainIdx[SeqIdx2];
	if (ChainIdx2 == UINT_MAX)
		return false;

	asserta(ChainIdx1 < SIZE(m_Chains));
	const PDBChain &Chain1 = *m_Chains[ChainIdx1];
	const string &Row1 = MSA.GetSeq(SeqIdx1);
	const uint L1 = Chain1.GetSeqLength();
	const string &Label1 = Chain1.m_Label;
	string U1;
	GetUngappedSeq(Row1, U1);
	asserta(U1 == Chain1.m_Seq);
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
	if (m_DoCore)
		GetAlignedPositions(Row1, Row2, Pos1s, Pos2s, &m_ColIsCore);
	else
		GetAlignedPositions(Row1, Row2, Pos1s, Pos2s, 0);
	Score = GetDALIScore_ChainPair(ChainIdx1, ChainIdx2, Pos1s, Pos2s);
	Z = GetDALIZFromScoreAndLengths(Score, L1, L2);
	return true;
	}

double DALIScorer::GetZ() const
	{
	const SeqDB &MSA = *m_MSA;
	const uint SeqCount = MSA.GetSeqCount();
	double SumZ = 0;
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
		for (uint SeqIdx2 = SeqIdx1 + 1; SeqIdx2 < SeqCount; ++SeqIdx2)
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
			if (m_DoCore)
				GetAlignedPositions(Row1, Row2, Pos1s, Pos2s, &m_ColIsCore);
			else
				GetAlignedPositions(Row1, Row2, Pos1s, Pos2s, 0);
			double Score = 
			  GetDALIScore_ChainPair(ChainIdx1, ChainIdx2, Pos1s, Pos2s);
			double Z = GetDALIZFromScoreAndLengths(Score, L1, L2);

			if (SeqIdx1 != SeqIdx2)
				{
				++PairCount;
				SumZ += Z;
				}
			}
		}
	if (PairCount == 0)
		return 0;
	return SumZ/PairCount;
	}

double DALIScorer::GetZ_Rows() const
	{
	const SeqDB &MSA = *m_MSA;
	const uint SeqCount = MSA.GetSeqCount();
	double SumZ = 0;
	uint PairCount = 0;
	for (uint SeqIdx1 = 0; SeqIdx1 < SeqCount; ++SeqIdx1)
		{
		for (uint SeqIdx2 = SeqIdx1 + 1; SeqIdx2 < SeqCount; ++SeqIdx2)
			{
			double Score, Z;
			bool Ok = GetDALIRowPair(SeqIdx1, SeqIdx2, Score, Z);
			if (!Ok)
				continue;
			SumZ += Z;
			++PairCount;
			}
		}
	if (PairCount == 0)
		return 0;
	return SumZ/PairCount;
	}

void DALIScorer::SetColToPosVec(bool Core)
	{
	const uint SeqCount = m_MSA->GetSeqCount();
	m_ColToPosVec.resize(SeqCount);
#pragma omp for
	for (int SeqIdx = 0; SeqIdx < int(SeqCount); ++SeqIdx)
		GetColToPos(uint(SeqIdx), m_ColToPosVec[SeqIdx], Core);
	}

void DALIScorer::GetColToPos(uint SeqIdx, vector<uint> &ColToPos, bool Core)
	{
	ColToPos.clear();
	const string &Row = m_MSA->GetSeq(SeqIdx);
	const uint ColCount = m_MSA->GetColCount();
	asserta(SIZE(Row) == ColCount);
	if (Core)
		asserta(SIZE(m_ColIsCore) == ColCount);
	uint Pos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Row[Col];

		if (islower(c) || isgap(c) || (Core && !m_ColIsCore[Col]))
			ColToPos.push_back(UINT_MAX);
		else
			ColToPos.push_back(Pos);

		if (!isgap(c))
			++Pos;
		}
	}

double DALIScorer::GetDALIScorePosPair_ById(
  uint ChainIdX, uint PosX1, uint PosX2,
  uint ChainIdY, uint PosY1, uint PosY2) const
	{
	const PDBChain &ChainX = *m_Chains[ChainIdX];
	const PDBChain &ChainY = *m_Chains[ChainIdY];
	//double dij_X = ChainX.GetDist(PosX1, PosX2);
	//double dij_Y = ChainY.GetDist(PosY1, PosY2);
	double dij_X = GetDist(ChainIdX, PosX1, PosX2);
	double dij_Y = GetDist(ChainIdY, PosY1, PosY2);
	double x = DALI_dpscorefun(dij_X, dij_Y);
	return x;
	}

//double DALIScorer::GetDALIScorePosPair(
//	const PDBChain &Q, uint PosQi, uint PosQj,
//	const PDBChain &T, uint PosTi, uint PosTj) const
//	{
//	double dij_Q = Q.GetDist(PosQi, PosQj);
//	double dij_T = T.GetDist(PosTi, PosTj);
//	double x = DALI_dpscorefun(dij_Q, dij_T);
//	return x;
//	}

double DALIScorer::GetDALIScoreColPair(uint Col1, uint Col2) const
	{
	const uint SeqCount = GetSeqCount();
	double Sum = 0;
	for (uint SeqX = 0; SeqX < SeqCount; ++SeqX)
		{
		uint ChainIdX = m_SeqIdxToChainIdx[SeqX];
		if (ChainIdX == UINT_MAX)
			continue;
		asserta(ChainIdX < SIZE(m_Chains));
		const PDBChain &ChainX = *m_Chains[ChainIdX];

		const vector<uint> &ColToPosX = m_ColToPosVec[SeqX];
		uint PosX1 = ColToPosX[Col1];
		uint PosX2 = ColToPosX[Col2];
		if (PosX1 == UINT_MAX || PosX2 == UINT_MAX)
			continue;

		for (uint SeqY = SeqX+1; SeqY < SeqCount; ++SeqY)
			{
			uint ChainIdY = m_SeqIdxToChainIdx[SeqY];
			if (ChainIdY == UINT_MAX)
				continue;
			asserta(ChainIdY < SIZE(m_Chains));
			const PDBChain &ChainY = *m_Chains[ChainIdY];

			const vector<uint> &ColToPosY = m_ColToPosVec[SeqY];
			uint PosY1 = ColToPosY[Col1];
			uint PosY2 = ColToPosY[Col2];
			if (PosY1 == UINT_MAX || PosY2 == UINT_MAX)
				continue;

			//Sum += GetDALIScorePosPair(
			//  ChainX, PosX1, PosX2, ChainY, PosY1, PosY2);
			Sum += GetDALIScorePosPair_ById(
			  ChainIdX, PosX1, PosX2, ChainIdY, PosY1, PosY2);
			}
		}
	return Sum;
	}

double DALIScorer::GetDiagScore() const
	{
	const uint SeqCount = GetSeqCount();
	double Sum = 0;
	for (uint i = 0; i < SeqCount; ++i)
		for (uint j = i+1; j < SeqCount; ++j)
			Sum += GetDiagScoreSeqPair(i, j);
	return Sum;
	}

double DALIScorer::GetDiagScoreSeqPair(uint SeqIdx1, uint SeqIdx2) const
	{
	if (m_SeqIdxToChainIdx[SeqIdx1] == UINT_MAX ||
	  m_SeqIdxToChainIdx[SeqIdx2] == UINT_MAX)
		return 0;
	const uint ColCount = GetColCount();
	const vector<uint> &ColToPosIdx1 = m_ColToPosVec[SeqIdx1];
	const vector<uint> &ColToPosIdx2 = m_ColToPosVec[SeqIdx2];
	asserta(SIZE(ColToPosIdx1) == ColCount);
	asserta(SIZE(ColToPosIdx2) == ColCount);
	uint Lali = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		if (ColToPosIdx1[Col] != UINT_MAX && ColToPosIdx2[Col] != UINT_MAX)
			++Lali;
	double DiagScore = Lali*g_DALI_Theta;
	return DiagScore;
	}

double DALIScorer::GetSumScore_Rows() const
	{
	const uint SeqCount = GetSeqCount();
	double SumScore = 0;
	uint PairCount = 0;
	for (uint SeqIdx1 = 0; SeqIdx1 < SeqCount; ++SeqIdx1)
		{
		for (uint SeqIdx2 = SeqIdx1 + 1; SeqIdx2 < SeqCount; ++SeqIdx2)
			{
			double Score, Z;
			bool Ok = GetDALIRowPair(SeqIdx1, SeqIdx2, Score, Z);
			if (Ok)
				SumScore += Score;
			}
		}
	return SumScore;
	}

double DALIScorer::GetSumScore_Cols() const
	{
	const uint SeqCount = GetSeqCount();
	const uint ColCount = GetColCount();
	double SumScore = 0;
	double SumScore_core = 0;
	double SumZ = 0;
	double SumZ_core = 0;
	uint PairCount = 0;
	double SumColScores = 0;
	for (uint Col1 = 0; Col1 < ColCount; ++Col1)
		{
		double SumCol1 = 0;
		for (uint Col2 = 0; Col2 < ColCount; ++Col2)
			{
			if (Col2 == Col1)
				continue;
			SumCol1 += GetDALIScoreColPair(Col1, Col2);
			}
		SumColScores += SumCol1;
		}
	double DiagScore = GetDiagScore();
	return SumColScores + DiagScore;
	}


const vector<vector<double> > &DALIScorer::GetDistMx(uint ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_DistMxVec));
	return m_DistMxVec[ChainIdx];
	}

double DALIScorer::GetDist(uint ChainId, uint Pos1, uint Pos2) const
	{
	assert(ChainId < SIZE(m_DistMxVec));
	const vector<vector<double> > &DistMx = m_DistMxVec[ChainId];
	const uint L = SIZE(DistMx);
	assert(Pos1 < L);
	const vector<double> &Row = DistMx[Pos1];
	assert(Pos2 < SIZE(Row));
	double d = Row[Pos2];
	return d;
	}

void DALIScorer::SetDistMx(uint ChainId)
	{
	asserta(ChainId < SIZE(m_DistMxVec));
	asserta(ChainId < SIZE(m_Chains));
	const PDBChain &Chain = *m_Chains[ChainId];
	const uint L = Chain.GetSeqLength();
	vector<vector<double> > &DistMx = m_DistMxVec[ChainId];
	DistMx.resize(L);
	for (uint i = 0; i < L; ++i)
		{
		DistMx[i].resize(L);
		DistMx[i][i] = 0;
		}

	for (uint i = 0; i < L; ++i)
		{
		for (uint j = 0; j < i; ++j)
			{
			double d = Chain.GetDist(i, j);
			DistMx[i][j] = d;
			DistMx[j][i] = d;
			}
		}
	}

void DALIScorer::SetDistMxs()
	{
	const uint ChainCount = SIZE(m_Chains);
	m_DistMxVec.clear();
	m_DistMxVec.resize(ChainCount);

	const uint N = SIZE(m_SeqIdxToChainIdx);
#pragma omp for
	for (int ii = 0; ii < int(N); ++ii)
		{
		uint i = uint(ii);
		uint ChainIdx = m_SeqIdxToChainIdx[i];
		if (ChainIdx == UINT_MAX)
			continue;
		SetDistMx(ChainIdx);
		}
	}
