#include "myutils.h"
#include "dssaligner.h"
#include "xdpmem.h"
#include <map>
#include <set>
#include "timing.h"

struct XDATA
	{
	const vector<vector<byte> > *ptrProfile1 = 0;
	const vector<vector<byte> > *ptrProfile2 = 0;
	DSSAligner *ptrDSS = 0;
	};

static float SubFn(void *UserData, uint Pos1, uint Pos2)
	{
	XDATA *XD = (XDATA *) UserData;
	DSSAligner &D = *XD->ptrDSS;
	float Score = D.GetScorePosPair(*XD->ptrProfile1, *XD->ptrProfile2, Pos1, Pos2);
	return Score;
	}

void MakeKmerToCoords(const vector<uint> &Kmers,
  map<uint, vector<uint> > &KmerToCoords)
	{
	KmerToCoords.clear();
	const uint n = SIZE(Kmers);
	for (uint Pos = 0; Pos < n; ++Pos)
		{
		uint Kmer = Kmers[Pos];
		if (KmerToCoords.find(Kmer) != KmerToCoords.end())
			KmerToCoords[Kmer].push_back(Pos);
		else
			{
			vector<uint> v;
			v.push_back(Pos);
			KmerToCoords[Kmer] = v;
			}
		}
	}

uint DSSAligner::GetU(const vector<uint> &Kmers1, const vector<uint> &Kmers2) const
	{
	set<uint> Set1;
	for (uint i = 0; i < SIZE(Kmers1); ++i)
		{
		uint Kmer = Kmers1[i];
		if (Kmer != UINT_MAX)
			Set1.insert(Kmer);
		}
	uint U = 0;
	for (uint i = 0; i < SIZE(Kmers2); ++i)
		{
		uint Kmer = Kmers2[i];
		if (Kmer != UINT_MAX && Set1.find(Kmer) != Set1.end())
			++U;
		}
	return U;
	}

bool DSSAligner::IsDiagPair(uint PosQi, uint PosQj, uint PosRi, uint PosRj) const
	{
	int spacing = abs(int(PosQi) - int(PosQj));
	if (spacing < m_Params->m_MinDiagSpacing || 
	  spacing > m_Params->m_MaxDiagSpacing)
		return false;
	int di = int(PosQi) - int(PosRi);
	int dj = int(PosQj) - int(PosRj);
	if (di == dj)
		return true;
	return false;
	}

bool DSSAligner::GetDiagPairs(
  const vector<uint> &SeedPosQs, const vector<uint> &SeedPosRs,
  vector<uint> &DiagPosQs, vector<uint> &DiagPosRs,
  vector<uint> &Lengths) const
	{
	DiagPosQs.clear();
	DiagPosRs.clear();
	Lengths.clear();
	const uint n = SIZE(SeedPosQs);
	asserta(SIZE(SeedPosRs) == n);
	for (uint i = 0; i < n; ++i)
		{
		uint PosQi = SeedPosQs[i];
		uint PosRi = SeedPosRs[i];
		for (uint j = i+1; j < n; ++j)
			{
			uint PosQj = SeedPosQs[j];
			uint PosRj = SeedPosRs[j];
			if (IsDiagPair(PosQi, PosQj, PosRi, PosRj))
				{
				DiagPosQs.push_back(min(PosQi, PosQj));
				DiagPosRs.push_back(min(PosRi, PosRj));
				uint Length = max(PosQi, PosQj) - min(PosQi, PosQj);
				Lengths.push_back(Length);
				}
			}
		}
	return SIZE(Lengths) > 0;
	}

void DSSAligner::GetSeeds(
	const vector<vector<byte> > &ProfileQ,
	const vector<vector<byte> > &ProfileR,
	const vector<uint> &KmersQ,
	const vector<uint> &KmersR,
	vector<uint> &SeedPosQs,
	vector<uint> &SeedPosRs) const
	{
	const uint QL = SIZE(ProfileQ[0]);
	const uint RL = SIZE(ProfileR[0]);
	asserta(SIZE(KmersQ) < QL);
	asserta(SIZE(KmersR) < RL);
	map<uint, vector<uint> > KmerToCoordsQ;
	map<uint, vector<uint> > KmerToCoordsR;
	MakeKmerToCoords(KmersQ, KmerToCoordsQ);
	MakeKmerToCoords(KmersR, KmerToCoordsR);
	const string &PatternStr = m_Params->m_PatternStr;
	const uint PatternLength = SIZE(PatternStr);
	const uint MAXNQNR = m_Params->m_MAXNQNR;
	const float MinKmerScore = m_Params->m_MinKmerScore;
	set<pair<uint, uint> > PairSet;
	for (map<uint, vector<uint> >::const_iterator iterq = KmerToCoordsQ.begin();
	  iterq != KmerToCoordsQ.end(); ++iterq)
		{
		uint Kmer = iterq->first;
		if (Kmer == UINT_MAX)
			continue;
		map<uint, vector<uint> >::const_iterator iterr = KmerToCoordsR.find(Kmer);
		if (iterr == KmerToCoordsR.end())
			continue;
		asserta(iterr->first == Kmer);
		const vector<uint> &CoordsQ = iterq->second;
		const vector<uint> &CoordsR = iterr->second;
		const uint nq = SIZE(CoordsQ);
		const uint nr = SIZE(CoordsR);
		if (nq*nr > MAXNQNR)
			continue;
		for (uint iq = 0; iq < nq; ++iq)
			{
			uint PosQ = CoordsQ[iq];
			for (uint ir = 0; ir < nr; ++ir)
				{
				uint PosR = CoordsR[ir];
				float Score = GetScoreSegPair(
				  ProfileQ, ProfileR, PosQ, PosR, PatternLength);
				if (Score >= MinKmerScore)
					{
					const pair<uint, uint> Pair(PosQ, PosR);
					asserta(PairSet.find(Pair) == PairSet.end());
					asserta(PosQ < QL);
					asserta(PosR < RL);
					SeedPosQs.push_back(PosQ);
					SeedPosRs.push_back(PosR);
					PairSet.insert(Pair);
					}
				}
			}
		}
	}

void DSSAligner::GetDiagSeedPairs(const vector<vector<byte> > &ProfileQ,
  const vector<vector<byte> > &ProfileR,
  const vector<uint> &KmersQ, const vector<uint> &KmersR,
  vector<uint> &DiagPosQs, vector<uint> &DiagPosRs, vector<uint> &DiagLengths) const
	{
	StartTimer(GetDiagSeedPairs);
	DiagPosQs.clear();
	DiagPosRs.clear();
	DiagLengths.clear();
	vector<uint> SeedPosQs;
	vector<uint> SeedPosRs;
	GetSeeds(ProfileQ, ProfileR, KmersQ, KmersR, SeedPosQs, SeedPosRs);
	GetDiagPairs(SeedPosQs, SeedPosRs, DiagPosQs, DiagPosRs, DiagLengths);
	EndTimer(GetDiagSeedPairs);
	}

//bool DSSAligner::MergeDiag(
//  uint PosQ1, uint PosR1, uint Len1,
//  uint &PosQ2, uint &PosR2, uint &Len2) const
//	{
//	if (PosQ1 == PosQ2)
//		return false;
//	int d1 = int(PosQ1) - int(PosR1);
//	int d2 = int(PosQ2) - int(PosR2);
//	if (d1 != d2)
//		return false;
//	if (PosQ1 < PosQ2)
//		{
//		Len2 += PosQ2 - PosQ1;
//		PosQ2 = PosQ1;
//		}
//	else
//		Len2 += PosQ1 - PosQ2;
//	return true;
//	}
//
//void DSSAligner::MergeDiagPairs(
//  const vector<uint> &ArgInDiagPosQs, const vector<uint> &ArgInDiagPosRs, const vector<uint> &ArgInLengths,
//  vector<uint> &OutDiagPosQs, vector<uint> &OutDiagPosRs, vector<uint> &OutLengths) const
//	{
//	OutDiagPosQs.clear();
//	OutDiagPosRs.clear();
//	OutLengths.clear();
//
//	vector<uint> InDiagPosQs = ArgInDiagPosQs;
//	vector<uint> InDiagPosRs = ArgInDiagPosRs;
//	vector<uint> InLengths = ArgInLengths;
//
//	const uint N = SIZE(InDiagPosQs);
//	asserta(SIZE(InDiagPosRs) == N);
//	asserta(SIZE(InLengths) == N);
//	for (uint i = 0; i < N; ++i)
//		{
//		uint InPosQ = InDiagPosQs[i];
//		uint InPosR = InDiagPosRs[i];
//		uint InLength = InLengths[i];
//		bool ThisMerged = false;
//		for (uint j = 0; j < SIZE(OutDiagPosQs); ++j)
//			{
//			bool Merged = MergeDiag(InPosQ, InPosR, InLength,
//				OutDiagPosQs[j], OutDiagPosRs[j], OutLengths[j]);
//			if (Merged)
//				{
//				ThisMerged = true;
//				break;
//				}
//			}
//		if (!ThisMerged)
//			{
//			OutDiagPosQs.push_back(InPosQ);
//			OutDiagPosRs.push_back(InPosR);
//			OutLengths.push_back(InLength);
//			}
//		}
//	}

float DSSAligner::AlignX(
	const PDBChain &ChainQ, const PDBChain &ChainR,
	const vector<uint> &KmersQ, const vector<uint> &KmersR,
	const vector<vector<byte> > &ProfileQ, 
	const vector<vector<byte> > &ProfileR)
	{
	m_BestHSPLo1 = UINT_MAX;
	m_BestHSPLo2 = UINT_MAX;
	m_BestHSPHi1 = UINT_MAX;
	m_BestHSPHi2 = UINT_MAX;
	m_BestHSPScore = FLT_MAX;
	m_BestHSPPath.clear();
	m_HSPScores.clear();

	const float X = m_Params->m_X;
	const float Open = m_Params->m_GapOpen;
	const float Ext = m_Params->m_GapExt;

	vector<uint> DiagPosQs;
	vector<uint> DiagPosRs;
	vector<uint> DiagLengths;
	GetDiagSeedPairs(ProfileQ, ProfileR, KmersQ, KmersR,
		DiagPosQs, DiagPosRs, DiagLengths);
	const uint DiagCount = SIZE(DiagLengths);
	asserta(SIZE(DiagPosQs) == DiagCount);
	asserta(SIZE(DiagPosRs) == DiagCount);
	if (DiagCount == 0)
		return -1;

	uint LQ = ChainQ.GetSeqLength();
	uint LR = ChainR.GetSeqLength();
	vector<uint> HSPLoQs;
	vector<uint> HSPLoRs;
	vector<uint> HSPHiQs;
	vector<uint> HSPHiRs;
	vector<string> HSPPaths;
	const float MinHSPScore = 8;
	float TopScore = 0;
	for (uint DiagIndex = 0; DiagIndex < DiagCount; ++DiagIndex)
		{
		uint LoQ = DiagPosQs[DiagIndex];
		uint LoR = DiagPosRs[DiagIndex];
		const uint HSPCount = SIZE(m_HSPScores);
		bool Skip = false;
		for (uint HSPIndex = 0; HSPIndex < HSPCount; ++HSPIndex)
			{
			const uint HSPLoQ = HSPLoQs[HSPIndex];
			const uint HSPHiQ = HSPHiQs[HSPIndex];
			if (LoQ >= HSPLoQ && LoQ <= HSPHiQ)
				{
				Skip = true;
				break;
				}

			const uint HSPLoR = HSPLoRs[HSPIndex];
			const uint HSPHiR = HSPHiRs[HSPIndex];
			if (LoR >= HSPLoR && LoR <= HSPHiR)
				{
				Skip = true;
				break;
				}
			}
		if (Skip)
			continue;

		string FwdPath;
		string BwdPath;
		XDATA XD;
		XD.ptrProfile1 = &ProfileQ;
		XD.ptrProfile2 = &ProfileR;
		XD.ptrDSS = this;
		float ScoreFwd = 0;
		float ScoreBwd = 0;
		uint FwdLoQ = LoQ + 2;
		uint FwdLoR = LoR + 2;
		uint BwdHiQ = LoQ + 1;
		uint BwdHiR = LoR + 1;
		ScoreFwd = XDropFwd(m_Mem, X, Open, Ext, SubFn, &XD,
		  FwdLoQ, LQ, FwdLoR, LR, FwdPath);
		uint ColCount = SIZE(FwdPath);
		ScoreBwd = XDropBwd(m_Mem, X, Open, Ext, SubFn, &XD,
		  BwdHiQ, LQ, BwdHiR, LR, BwdPath);
		float HSPScore = ScoreFwd + ScoreBwd;
		if (HSPScore < MinHSPScore)
			continue;

		uint HSPLoQ, HSPLoR, HSPHiQ, HSPHiR;
		string HSPPath;
		MergeFwdBwd(LQ, LR,
		  FwdLoQ, FwdLoR, FwdPath,
		  BwdHiQ, BwdHiR, BwdPath,
		  HSPLoQ, HSPLoR, HSPHiQ, HSPHiR, HSPPath);
		if (HSPScore > TopScore)
			{
			TopScore = HSPScore;
			m_BestHSPLo1 = HSPLoQ;
			m_BestHSPLo2 = HSPLoR;
			m_BestHSPHi1 = HSPHiQ;
			m_BestHSPHi2 = HSPHiR;
			m_BestHSPScore = HSPScore;
			m_BestHSPPath = HSPPath;
			}
		HSPLoQs.push_back(HSPLoQ);
		HSPLoRs.push_back(HSPLoR);
		HSPHiQs.push_back(HSPHiQ);
		HSPHiRs.push_back(HSPHiR);
		m_HSPScores.push_back(HSPScore);
		HSPPaths.push_back(HSPPath);
		}
	return TopScore;
	}