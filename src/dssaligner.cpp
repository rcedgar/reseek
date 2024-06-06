#include "myutils.h"
#include "mx.h"
#include "dssaligner.h"
#include "pdbchain.h"
#include "alpha.h"
#include "tma.h"
#include "xdpmem.h"
#include "timing.h"
#include "combomx.h"
#include "cigar.h"
#include <thread>
#include <set>
#include <mutex>

mutex DSSAligner::m_OutputLock;
mutex DSSAligner::m_StatsLock;

uint SWFastPinopGapless(const int8_t * const *AP, uint LA,
  const int8_t *B, uint LB);
void PrettyAln(FILE *f, const PDBChain &A, const PDBChain &B,
  uint LoA, uint LoB, const string &Path, float Evalue);
void LogAln(const char *A, const char *B, const char *Path, unsigned ColCount);
float SWFast(XDPMem &Mem, const Mx<float> &SMx, uint LA, uint LB,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path);
float SWFastGapless(XDPMem &Mem, const Mx<float> &SMx, uint LA, uint LB,
  uint &Besti, uint &Bestj);
float SWFastGaplessProf(XDPMem &Mem, const float * const *ProfA, uint LA,
  const byte *B, uint LB, uint &Besti, uint &Bestj);
double GetDALIScore_Path(const PDBChain &Q, const PDBChain &T,
  const string &Path, uint LoQ, uint LoT);
float SWFastGapless(XDPMem &Mem, const Mx<float> &SMx, uint LA, uint LB,
  uint &Besti, uint &Bestj);
int SWFastGapless_Int(XDPMem &Mem, const Mx<int8_t> &SMx, uint LA, uint LB,
  uint &Besti, uint &Bestj);
void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);
float SWFastGaplessProfb(float *DProw_, const float * const *ProfA, uint LA, const byte *B, uint LB);

uint DSSAligner::m_AlnCount;
uint DSSAligner::m_SWCount;
uint DSSAligner::m_ComboFilterCount;
uint DSSAligner::m_UFilterCount;
uint DSSAligner::m_ParasailSaturateCount;

uint GetU(const vector<uint> &KmersQ, const vector<uint> &KmersR)
	{
	set<uint> SetQ;
	set<uint> SetR;
	for (uint i = 0; i < SIZE(KmersQ); ++i)
		SetQ.insert(KmersQ[i]);
	for (uint i = 0; i < SIZE(KmersR); ++i)
		SetR.insert(KmersR[i]);
	uint U = 0;
	for (set<uint>::const_iterator iter = SetQ.begin();
	  iter != SetQ.end(); ++iter)
		if (SetR.find(*iter) != SetR.end())
			++U;
	return U;
	}

uint GetUBits(const vector<uint> &KmerBitsQ, const vector<uint> &KmerBitsR)
	{
	const uint DictSize = 36*36;
	const uint DictSizeWords = 1 + (DictSize - 1)/32;
	asserta(SIZE(KmerBitsQ) == DictSizeWords);
	asserta(SIZE(KmerBitsR) == DictSizeWords);
	uint n = 0;
	for (uint i = 0; i < DictSizeWords; ++i)
		{
		uint WordQ = KmerBitsQ[i];
		uint WordR = KmerBitsR[i];
		uint And = (WordQ & WordR);
		uint Bit = 1;
		for (uint j = 0; j < 32; ++j)
			{
			if (And & Bit)
				++n;
			Bit <<= 1;
			}
		}
	return n;
	}

uint GetMatchColCount(const string &Path)
	{
	uint M = 0;
	for (uint i = 0; i < SIZE(Path); ++i)
		if (Path[i] == 'M')
			++M;
	return M;
	}

float DSSAligner::GetEvaluePath(
  const PDBChain &ChainA, const PDBChain &ChainB,
  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB,
  uint LoA, uint LoB, const string &Path) const
	{
	const uint LA = ChainA.GetSeqLength();
	float DPScore = GetDPScorePath(ProfileA, ProfileB, LoA, LoB, Path);
	float AlnDALIScore = 0;
	if (m_Params->m_DALIw != 0)
		AlnDALIScore = (float)
		  GetDALIScore_Path(ChainA, ChainB, Path, LoA, LoB);
	float ScoreDiff = DPScore;
	uint M = GetMatchColCount(Path);
	float TestStat = ScoreDiff + M*m_Params->m_FwdMatchScore +
	  m_Params->m_DALIw*AlnDALIScore/(m_Params->m_Lambda + LA);
	float E = m_Params->GetEvalue(TestStat);
	return E; 
	}

int DSSAligner::GetComboDPScorePathInt(const vector<byte> &ComboLettersA,
  const vector<byte> &ComboLettersB, uint LoA, uint LoB,
  const string &Path) const
	{
	uint Sum = 0;
	uint PosA = LoA;
	uint PosB = LoB;
	const int Open = -m_Params->m_ParaComboGapOpen;
	const int Ext = -m_Params->m_ParaComboGapExt;
	const float FwdMatchScore = m_Params->m_FwdMatchScore;
	const uint ColCount = SIZE(Path);
	extern int8_t IntScoreMx_Combo[36][36];
 
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			{
			asserta(PosA < SIZE(ComboLettersA));
			asserta(PosB < SIZE(ComboLettersB));
			byte a = ComboLettersA[PosA];
			byte b = ComboLettersB[PosB];
			Sum += IntScoreMx_Combo[a][b];
			++PosA;
			++PosB;
			break;
			}

		case 'D':
			if (Col != 0 && Path[Col-1] == 'D')
				Sum += Ext;
			else
				Sum += Open;
			++PosA;
			break;

		case 'I':
			if (Col != 0 && Path[Col-1] == 'I')
				Sum += Ext;
			else
				Sum += Open;
			++PosB;
			break;

		default:
			asserta(false);
			}
		}
	return Sum;
	}

float DSSAligner::GetComboDPScorePath(const vector<byte> &LettersA,
	const vector<byte> &LettersB, uint LoA, uint LoB,
	float GapOpen, float GapExt, const string &Path) const
	{
	asserta(GapOpen <= 0);
	asserta(GapExt <= 0);
	float Sum = 0;
	uint PosA = LoA;
	uint PosB = LoB;
	const uint ColCount = SIZE(Path);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			{
			asserta(PosA < SIZE(LettersA));
			asserta(PosB < SIZE(LettersB));
			uint LetterA = LettersA[PosA];
			uint LetterB = LettersB[PosB];
			asserta(LetterA < 36);
			asserta(LetterB < 36);
			extern float ScoreMx_Combo[36][36];
			Sum += ScoreMx_Combo[LetterA][LetterB];
			++PosA;
			++PosB;
			break;
			}

		case 'D':
			if (Col != 0 && Path[Col-1] == 'D')
				Sum += GapExt;
			else
				Sum += GapOpen;
			++PosA;
			break;

		case 'I':
			if (Col != 0 && Path[Col-1] == 'I')
				Sum += GapExt;
			else
				Sum += GapOpen;
			++PosB;
			break;

		default:
			asserta(false);
			}
		}
	return Sum;
	}

//float DSSAligner::GetDPScorePathVecs(
//	const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB,
//	uint LoA, uint LoB, const string &Path,
//	vector<uint> &PosAs, vector<uint> &PosBs, vector<float> &ColScores) const
//	{
//	const uint ColCount = SIZE(Path);
//	PosAs.resize(ColCount);
//	PosBs.resize(ColCount);
//	ColScores.resize(ColCount);
//	float Sum = 0;
//	uint PosA = LoA;
//	uint PosB = LoB;
//	const float Open = m_Params->m_GapOpen;
//	const float Ext = m_Params->m_GapExt;
//	for (uint Col = 0; Col < ColCount; ++Col)
//		{
//		char c = Path[Col];
//		switch (c)
//			{
//		case 'M':
//			{
//			float ColScore = GetScorePosPair(ProfileA, ProfileB, PosA, PosB);
//			Sum += ColScore;
//			ColScores[Col] = ColScore;
//			PosAs[Col] = PosA;
//			PosBs[Col] = PosB;
//			++PosA;
//			++PosB;
//			break;
//			}
//
//		case 'D':
//			{
//			ColScores[Col] = Open/2;
//			PosAs[Col] = PosA;
//			PosBs[Col] = UINT_MAX;
//			if (Col != 0 && Path[Col-1] == 'D')
//				Sum += Ext;
//			else
//				Sum += Open;
//			++PosA;
//			break;
//			}
//
//		case 'I':
//			{
//			ColScores[Col] = Open/2;
//			PosAs[Col] = UINT_MAX;
//			PosBs[Col] = PosB;
//			if (Col != 0 && Path[Col-1] == 'I')
//				Sum += Ext;
//			else
//				Sum += Open;
//			++PosB;
//			break;
//			}
//
//		default:
//			asserta(false);
//			}
//		}
//	return Sum;
//	}

// GetDPScorePath calculates AlnScore which is optimized by SWFast.
// The test statistic is 
//	TS = AlnScore + ColCount*g_FwdMatchScore + g_DALIw*FwdDaliScore;
float DSSAligner::GetDPScorePath(const vector<vector<byte> > &ProfileA,
  const vector<vector<byte> > &ProfileB, uint LoA, uint LoB,
  const string &Path) const
	{
	float Sum = 0;
	uint PosA = LoA;
	uint PosB = LoB;
	const float Open = m_Params->m_GapOpen;
	const float Ext = m_Params->m_GapExt;
	const uint ColCount = SIZE(Path);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			Sum += GetScorePosPair(ProfileA, ProfileB, PosA, PosB);
			++PosA;
			++PosB;
			break;

		case 'D':
			if (Col != 0 && Path[Col-1] == 'D')
				Sum += Ext;
			else
				Sum += Open;
			++PosA;
			break;

		case 'I':
			if (Col != 0 && Path[Col-1] == 'I')
				Sum += Ext;
			else
				Sum += Open;
			++PosB;
			break;

		default:
			asserta(false);
			}
		}
	return Sum;
	}

float DSSAligner::GetScorePosPair(const vector<vector<byte> > &ProfileA,
  const vector<vector<byte> > &ProfileB, uint PosA, uint PosB) const
	{
	const DSSParams &Params = *m_Params;
	const uint FeatureCount = Params.GetFeatureCount();
	asserta(SIZE(ProfileA) == FeatureCount);
	asserta(SIZE(ProfileB) == FeatureCount);
	float MatchScore = 0;
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		asserta(PosA < SIZE(ProfileA[FeatureIdx]));
		asserta(PosB < SIZE(ProfileB[FeatureIdx]));
		//float w = m_Params->m_Weights[FeatureIdx];
		FEATURE F = m_Params->m_Features[FeatureIdx];
		uint AlphaSize = g_AlphaSizes2[F];
		//float **ScoreMx = g_ScoreMxs2[F];
		float **ScoreMx = m_Params->m_ScoreMxs[F];
		const vector<byte> &ProfRowA = ProfileA[FeatureIdx];
		const vector<byte> &ProfRowB = ProfileB[FeatureIdx];
		byte ia = ProfRowA[PosA];
		assert(ia < AlphaSize);
		const float *ScoreMxRow = ScoreMx[ia];
		byte ib = ProfRowB[PosB];
		assert(ib < AlphaSize);
		MatchScore += ScoreMxRow[ib];
		}
	return MatchScore;
	}

float DSSAligner::GetScoreSegPair(const vector<vector<byte> > &ProfileA,
  const vector<vector<byte> > &ProfileB, uint PosA, uint PosB, uint n) const
	{
	float Score = 0;
	for (uint i = 0; i < n; ++i)
		Score += GetScorePosPair(ProfileA, ProfileB, PosA+i, PosB+i);
	return Score;
	}

void DSSAligner::SetSMx_YesRev()
	{
	const DSSParams &Params = *m_Params;
	const vector<vector<byte> > &ProfileA = *m_ProfileA;
	const vector<vector<byte> > &ProfileB = *m_ProfileB;
	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	Mx<float> &SMx = m_SMx;
	SMx.Alloc("SMx", LA, LB);
	Mx<float> &RevSMx = m_RevSMx;
	RevSMx.Alloc("RevSMx", LA, LB);
	StartTimer(SetSMx_YesRev);
	float **Sim = SMx.GetData();
	float **RevSim = RevSMx.GetData();

	const uint FeatureCount = Params.GetFeatureCount();
	asserta(SIZE(ProfileA) == FeatureCount);
	asserta(SIZE(ProfileB) == FeatureCount);

// Special case first feature because = not += and
	FEATURE F0 = m_Params->m_Features[0];
	uint AlphaSize0 = g_AlphaSizes2[F0];
	float **ScoreMx0 = m_Params->m_ScoreMxs[F0];
	const vector<byte> &ProfRowA = ProfileA[0];
	const vector<byte> &ProfRowB = ProfileB[0];
	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		byte ia = ProfRowA[PosA];
		float *SimRow = Sim[PosA];
		float *RevSimRow = RevSim[PosA];
		const float *ScoreMxRow = ScoreMx0[ia];

		for (uint PosB = 0; PosB < LB; ++PosB)
			{
			byte ib = ProfRowB[PosB];
			assert(ia < AlphaSize0 && ib < AlphaSize0);
			float MatchScore = ScoreMxRow[ib];
			SimRow[PosB] = MatchScore;
			RevSimRow[LB-1-PosB] = MatchScore;
			}
		}

	for (uint FeatureIdx = 1; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = m_Params->m_Features[FeatureIdx];
		uint AlphaSize = g_AlphaSizes2[F];
		float **ScoreMx = m_Params->m_ScoreMxs[F];
		const vector<byte> &ProfRowA = ProfileA[FeatureIdx];
		const vector<byte> &ProfRowB = ProfileB[FeatureIdx];
		for (uint PosA = 0; PosA < LA; ++PosA)
			{
			byte ia = ProfRowA[PosA];
			float *SimRow = Sim[PosA];
			float *RevSimRow = RevSim[PosA];
			const float *ScoreMxRow = ScoreMx[ia];

			for (uint PosB = 0; PosB < LB; ++PosB)
				{
				byte ib = ProfRowB[PosB];
				float MatchScore = 0;
				assert(ia < AlphaSize && ib < AlphaSize);
				MatchScore = ScoreMxRow[ib];
				SimRow[PosB] += MatchScore;
				RevSimRow[LB-1-PosB] += MatchScore;
				}
			}
		}
		EndTimer(SetSMx_YesRev);
#if DEBUG
	{
	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		for (uint PosB = 0; PosB < LB; ++PosB)
			{
			float MatchScore = Sim[PosA][PosB];
			float MatchScore2 = GetScorePosPair(ProfileA, ProfileB, PosA, PosB);
			asserta(feq(MatchScore2, MatchScore));
			}
		}
	}
#endif
	}

void DSSAligner::SetSMx_Combo_Int()
	{
	const vector<byte> &ComboLettersA = *m_ComboLettersA;
	const vector<byte> &ComboLettersB = *m_ComboLettersB;
	const uint LA = SIZE(ComboLettersA);
	const uint LB = SIZE(ComboLettersB);
	const uint AS = 36; // g_AlphaSizes2[FEATURE_Combo]; @@TODO
	asserta(AS == 36);

	Mx<int8_t> &SMx = m_SMx_Int;
	SMx.Alloc("SMx_Int", LA, LB);
	Mx<int8_t> &RevSMx = m_RevSMx_Int;
	RevSMx.Alloc("RevSMx_Int", LA, LB);
	StartTimer(SetSMx_Combo_Int);
	int8_t **Sim = SMx.GetData();
	int8_t **RevSim = RevSMx.GetData();

	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		byte LetterA = ComboLettersA[PosA];
		asserta(LetterA < AS);
		const int8_t *ComboMxRow = IntScoreMx_Combo[LetterA];
		int8_t *SimRow = Sim[PosA];
		int8_t *RevSimRow = RevSim[PosA];
		for (uint PosB = 0; PosB < LB; ++PosB)
			{
			byte LetterB = ComboLettersB[PosB];
			asserta(LetterB < AS);
			int8_t Score = ComboMxRow[LetterB];
			SimRow[PosB] = Score;
			RevSimRow[LB-1-PosB] = Score;
			}
		}
	EndTimer(SetSMx_Combo_Int);
	}

void DSSAligner::SetComboQP()
	{
	StartTimer(SetComboQP);
	const vector<byte> &ComboLettersA = *m_ComboLettersA;
	uint LA = SIZE(ComboLettersA);
	uint n = SIZE(m_ProfCombo);
	if (LA > n)
		{
		m_ProfCombo.resize(LA);
		m_ProfComboRev.resize(LA);
		}
	const uint AS = 36; // g_AlphaSizes2[FEATURE_Combo]; @@TODO
	asserta(AS == 36);
	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		byte LetterA = ComboLettersA[PosA];
		asserta(LetterA < AS);
		const float *ComboMxRow = ScoreMx_Combo[LetterA];
		m_ProfCombo[PosA] = ComboMxRow;
		m_ProfComboRev[LA-PosA-1] = ComboMxRow;
		}
	EndTimer(SetComboQP);
	}

void DSSAligner::SetComboQPi()
	{
	const vector<byte> &ComboLettersA = *m_ComboLettersA;
	uint LA = SIZE(ComboLettersA);
	uint n = SIZE(m_ProfComboi);
	if (LA > n)
		{
		m_ProfComboi.resize(LA);
		m_ProfComboRevi.resize(LA);
		}
	const uint AS = 36; // g_AlphaSizes2[FEATURE_Combo]; @@TODO
	asserta(AS == 36);
	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		byte LetterA = ComboLettersA[PosA];
		asserta(LetterA < AS);
		const int8_t *ComboMxRow = IntScoreMx_Combo[LetterA];
		m_ProfComboi[PosA] = ComboMxRow;
		m_ProfComboRevi[LA-PosA-1] = ComboMxRow;
		}
	}

void DSSAligner::SetSMx_Combo()
	{
	const vector<byte> &ComboLettersA = *m_ComboLettersA;
	const vector<byte> &ComboLettersB = *m_ComboLettersB;
	const uint LA = SIZE(ComboLettersA);
	const uint LB = SIZE(ComboLettersB);
	const uint AS = 36; // g_AlphaSizes2[FEATURE_Combo]; @@TODO
	asserta(AS == 36);

	Mx<float> &SMx = m_SMx;
	SMx.Alloc("SMx", LA, LB);
	Mx<float> &RevSMx = m_RevSMx;
	RevSMx.Alloc("RevSMx", LA, LB);
	StartTimer(SetSMx_Combo);
	float **Sim = SMx.GetData();
	float **RevSim = RevSMx.GetData();

	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		byte LetterA = ComboLettersA[PosA];
		asserta(LetterA < AS);
		const float *ComboMxRow = ScoreMx_Combo[LetterA];
		float *SimRow = Sim[PosA];
		float *RevSimRow = RevSim[PosA];
		for (uint PosB = 0; PosB < LB; ++PosB)
			{
			byte LetterB = ComboLettersB[PosB];
			asserta(LetterB < AS);
			float Score = ComboMxRow[LetterB];
			SimRow[PosB] = Score;
			RevSimRow[LB-1-PosB] = Score;
			}
		}
	EndTimer(SetSMx_Combo);
	}

void DSSAligner::SetSMx_NoRev()
	{
	const DSSParams &Params = *m_Params;
	const vector<vector<byte> > &ProfileA = *m_ProfileA;
	const vector<vector<byte> > &ProfileB = *m_ProfileB;
	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	Mx<float> &SMx = m_SMx;
	SMx.Alloc("SMx", LA, LB);
	StartTimer(SetSMx_NoRev);
	float **Sim = SMx.GetData();

	const uint FeatureCount = Params.GetFeatureCount();
	asserta(SIZE(ProfileA) == FeatureCount);
	asserta(SIZE(ProfileB) == FeatureCount);

// Special case first feature because = not += and
	FEATURE F0 = m_Params->m_Features[0];
	uint AlphaSize0 = g_AlphaSizes2[F0];
	float **ScoreMx0 = m_Params->m_ScoreMxs[F0];
	const vector<byte> &ProfRowA = ProfileA[0];
	const vector<byte> &ProfRowB = ProfileB[0];
	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		byte ia = ProfRowA[PosA];
		float *SimRow = Sim[PosA];
		assert(ia < AlphaSize0);
		const float *ScoreMxRow = ScoreMx0[ia];

		for (uint PosB = 0; PosB < LB; ++PosB)
			{
			byte ib = ProfRowB[PosB];
			assert(ia < AlphaSize0 && ib < AlphaSize0);
			float MatchScore = ScoreMxRow[ib];
			SimRow[PosB] = MatchScore;
			}
		}

	for (uint FeatureIdx = 1; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = m_Params->m_Features[FeatureIdx];
		uint AlphaSize = g_AlphaSizes2[F];
		float **ScoreMx = m_Params->m_ScoreMxs[F];
		const vector<byte> &ProfRowA = ProfileA[FeatureIdx];
		const vector<byte> &ProfRowB = ProfileB[FeatureIdx];
		for (uint PosA = 0; PosA < LA; ++PosA)
			{
			byte ia = ProfRowA[PosA];
			assert(ia < AlphaSize);
			const float *ScoreMxRow = ScoreMx[ia];
			float *SimRow = Sim[PosA];

			for (uint PosB = 0; PosB < LB; ++PosB)
				{
				byte ib = ProfRowB[PosB];
				assert(ib < AlphaSize);
				float MatchScore = ScoreMxRow[ib];
				SimRow[PosB] += MatchScore;
				}
			}
		}
	EndTimer(SetSMx_NoRev);
// GetScorePosPair bug?
#if 0 // DEBUG
	{
	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		for (uint PosB = 0; PosB < LB; ++PosB)
			{
			float MatchScore = Sim[PosA][PosB];
			float MatchScore2 = GetScorePosPair(ProfileA, ProfileB, PosA, PosB);
			asserta(feq(MatchScore2, MatchScore));
			}
		}
	}
#endif
	}

float DSSAligner::GetComboScore()
	{
	float ComboScore = AlignComboQP(*m_ComboLettersA, *m_ComboLettersB);
	return ComboScore;
	}

bool DSSAligner::ComboFilter()
	{
	float MCS = m_Params->m_Omega;
	if (MCS <= 0)
		return true;
	float ComboScore = GetComboScore(); // AlignComboQP(*m_ComboLettersA, *m_ComboLettersB);
	if (ComboScore < MCS)
		return false;
	return true;
	}

bool DSSAligner::UFilter()
	{
	uint MinU = m_Params->m_MinU;
	if (MinU == 0)
		return true;
	uint U = GetUBits(*m_ComboKmerBitsA, *m_ComboKmerBitsB);
	if (U < MinU)
		return false;
	return true;
	}

void DSSAligner::Align_ComboFilter(
  const PDBChain &ChainA, const PDBChain &ChainB,
  const vector<byte> &ComboLettersA, const vector<byte> &ComboLettersB,
  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB)
	{
	SetQuery(ChainA, &ProfileA, 0, &ComboLettersA);
	SetTarget(ChainB, &ProfileB, 0, &ComboLettersB);

	m_EvalueAB = FLT_MAX;
	m_EvalueBA = FLT_MAX;
	m_PathAB.clear();

	m_StatsLock.lock();
	++m_AlnCount;
	m_StatsLock.unlock();

	bool ComboFilterOk = ComboFilter();
	if (!ComboFilterOk)
		{
		m_StatsLock.lock();
		++m_ComboFilterCount;
		m_StatsLock.unlock();
		return;
		}
	Align_NoAccel();
	}

void DSSAligner::SetQuery(
	const PDBChain &Chain,
	const vector<vector<byte> > *ptrProfile,
	const vector<uint> *ptrComboKmerBits,
	const vector<byte> *ptrComboLetters)
	{
	m_ChainA = &Chain;
	m_ProfileA = ptrProfile;
	m_ComboKmerBitsA = ptrComboKmerBits;
	m_ComboLettersA = ptrComboLetters;
	if (m_Params->m_UsePara)
		SetComboQP_Para();
	else
		SetComboQP();
	}

void DSSAligner::SetTarget(
	const PDBChain &Chain,
	const vector<vector<byte> > *ptrProfile,
	const vector<uint> *ptrComboKmerBits,
	const vector<byte> *ptrComboLetters)
	{
	m_ChainB = &Chain;
	m_ProfileB = ptrProfile;
	m_ComboKmerBitsB = ptrComboKmerBits;
	m_ComboLettersB = ptrComboLetters;
	}

void DSSAligner::AlignComboOnly()
	{
	m_EvalueAB = FLT_MAX;
	m_EvalueBA = FLT_MAX;
	m_PathAB.clear();

	m_StatsLock.lock();
	++m_AlnCount;
	m_StatsLock.unlock();

	bool UFilterOk = UFilter();
	if (!UFilterOk)
		{
		m_StatsLock.lock();
		++m_UFilterCount;
		m_StatsLock.unlock();
		return;
		}
	bool ComboFilterOk = ComboFilter();
	if (!ComboFilterOk)
		{
		m_StatsLock.lock();
		++m_ComboFilterCount;
		m_StatsLock.unlock();
		return;
		}
	AlignComboPath();
	}

void DSSAligner::AlignQueryTarget()
	{
	m_EvalueAB = FLT_MAX;
	m_EvalueBA = FLT_MAX;
	m_PathAB.clear();

	m_StatsLock.lock();
	++m_AlnCount;
	m_StatsLock.unlock();

	bool UFilterOk = UFilter();
	if (!UFilterOk)
		{
		m_StatsLock.lock();
		++m_UFilterCount;
		m_StatsLock.unlock();
		return;
		}
	bool ComboFilterOk = ComboFilter();
	if (!ComboFilterOk)
		{
		m_StatsLock.lock();
		++m_ComboFilterCount;
		m_StatsLock.unlock();
		return;
		}
	Align_NoAccel();
	}

void DSSAligner::CalcEvalue()
	{
// MinFwdScore small speedup by avoiding call to GetDALIScore_Path()
	float AlnDALIScore = 0;
	if (m_AlnFwdScore >= m_Params->m_MinFwdScore)
		{
		StartTimer(DALIScore);
		if (m_Params->m_DALIw != 0)
			AlnDALIScore = (float)
			  GetDALIScore_Path(*m_ChainA, *m_ChainB, m_PathAB, m_LoA, m_LoB);
		EndTimer(DALIScore);
		}

	const float DALIw = m_Params->m_DALIw;
	const float FwdMatchScore = m_Params->m_FwdMatchScore;
	uint M = GetMatchColCount(m_PathAB);
	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();
	uint Lambda = m_Params->m_Lambda;
	float StatTop = m_AlnFwdScore + M*FwdMatchScore + DALIw*AlnDALIScore;

	if (StatTop < 0)
		StatTop = 0;
	m_TestStatisticAB = StatTop/(LA + Lambda);
	m_TestStatisticBA = StatTop/(LB + Lambda);

	m_EvalueAB = m_Params->GetEvalue(m_TestStatisticAB);
	m_EvalueBA = m_Params->GetEvalue(m_TestStatisticBA);
	}

void DSSAligner::Align_NoAccel()
	{
	m_PathAB.clear();
	m_LoA = UINT_MAX;
	m_LoB = UINT_MAX;

	SetSMx_NoRev();

	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	StartTimer(SWFwd);
	uint Leni, Lenj;
	m_AlnFwdScore = SWFast(m_Mem, m_SMx, LA, LB,
	  m_Params->m_GapOpen, m_Params->m_GapExt,
	  m_LoA, m_LoB, Leni, Lenj, m_PathAB);
	EndTimer(SWFwd);

	CalcEvalue();
	}

void DSSAligner::AlignComboPath()
	{
	m_PathAB.clear();
	m_LoA = UINT_MAX;
	m_LoB = UINT_MAX;
	m_AlnFwdScore = 0;
	float ComboFwdScore = AlignComboQP_Para_Path(m_LoA, m_LoB, m_PathAB);
	if (ComboFwdScore < 0)
		{
		m_StatsLock.lock();
		++m_ParasailSaturateCount;
		m_StatsLock.unlock();
		SetSMx_NoRev();
		Mx<float> &SMx = m_SMx;
		Mx<float> &RevSMx = m_RevSMx;
		ComboFwdScore = 
		  AlignCombo(*m_ComboLettersA, *m_ComboLettersB, m_LoA, m_LoB, m_PathAB);
		m_AlnFwdScore =
		  GetDPScorePath(*m_ProfileA, *m_ProfileB, m_LoA, m_LoB, m_PathAB);
		}
	CalcEvalue();
	}

void DSSAligner::ToAln(FILE *f, float MaxEvalue)
	{
	if (f == 0)
		return;
	if (m_EvalueAB > MaxEvalue)
		return;
	PrettyAln(f, *m_ChainA, *m_ChainB, m_LoA, m_LoB, m_PathAB, m_EvalueAB);
	}

void DSSAligner::ToAlnBA(FILE *f, float MaxEvalue)
	{
	if (f == 0)
		return;
	if (m_EvalueBA > MaxEvalue)
		return;
	string PathBA;
	for (uint i = 0; i < SIZE(m_PathAB); ++i)
		{
		char c = m_PathAB[i];
		if (c == 'D')
			c = 'I';
		else if (c == 'I')
			c = 'D';
		PathBA.push_back(c);
		}
	PrettyAln(f, *m_ChainB, *m_ChainA, m_LoA, m_LoB, PathBA, m_EvalueBA);
	}

void DSSAligner::ToFasta2(FILE *f, float MaxEvalue)
	{
	if (f == 0)
		return;
	if (m_EvalueAB > MaxEvalue)
		return;
	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	const char *LabelA = m_ChainA->m_Label.c_str();
	const char *LabelB = m_ChainB->m_Label.c_str();

	const string &SeqA = m_ChainA->m_Seq;
	const string &SeqB = m_ChainB->m_Seq;

	uint MaxLo = max(m_LoA, m_LoB);

	string RowA;
	string RowB;

	for (uint i = 0; i < MaxLo - m_LoA; ++i)
		RowA.push_back('.');
	for (uint i = 0; i < m_LoA; i++)
		RowA.push_back(SeqA[i]);

	for (uint i = 0; i < MaxLo - m_LoB; ++i)
		RowB.push_back('.');
	for (uint i = 0; i < m_LoB; i++)
		RowB.push_back(SeqB[i]);

	const uint ColCount = SIZE(m_PathAB);
	uint PosA = m_LoA;
	uint PosB = m_LoB;
	uint IdCount = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = m_PathAB[Col];
		switch (c)
			{
		case 'M':
			{
			asserta(PosA < SIZE(SeqA));
			asserta(PosB < SIZE(SeqB));
			char a = SeqA[PosA++];
			char b = SeqB[PosB++];
			if (a == b)
				++IdCount;
			RowA += a;
			RowB += b;
			break;
			}

		case 'D':
			{
			asserta(PosA < SIZE(SeqA));
			RowA += SeqA[PosA++];
			RowB += '-';
			break;
			}

		case 'I':
			{
			asserta(PosB < SIZE(SeqB));
			RowA += '-';
			RowB += SeqB[PosB++];
			break;
			}

		default:
			asserta(false);
			}
		}

	string LabelAx = LabelA;
	Psa(LabelAx, " E=%.3g Id=%.1f%%",
	  m_EvalueAB, GetPct(IdCount, ColCount));
	LabelAx += " (";
	LabelAx += LabelB;
	LabelAx += ")";
	m_OutputLock.lock();
	SeqToFasta(f, LabelAx.c_str(), RowA);
	SeqToFasta(f, LabelB, RowB);
	m_OutputLock.unlock();
	}

void DSSAligner::ToTsv(FILE *f, float MaxEvalue)
	{
	if (f == 0)
		return;
	if (m_EvalueAB > MaxEvalue)
		return;
	string CIGAR;
	PathToCIGAR(m_PathAB.c_str(), CIGAR);
	uint M, D, I;
	GetPathCounts(m_PathAB, M, D, I);
	m_OutputLock.lock();
	fprintf(f, "%.3g", m_EvalueAB);
	fprintf(f, "\t%s", m_ChainA->m_Label.c_str());
	fprintf(f, "\t%s", m_ChainB->m_Label.c_str());
	fprintf(f, "\t%u", m_LoA + 1);
	fprintf(f, "\t%u", m_LoA + M + D);
	fprintf(f, "\t%u", m_LoB + 1);
	fprintf(f, "\t%u", m_LoB + M + I);
	fprintf(f, "\t%s\n", CIGAR.c_str());
	m_OutputLock.unlock();
	}

void DSSAligner::ToTsvBA(FILE *f, float MaxEvalue)
	{
	if (f == 0)
		return;
	if (m_EvalueBA > MaxEvalue)
		return;
	string CIGAR;
	PathToCIGAR(m_PathAB.c_str(), CIGAR, true);
	uint M, D, I;
	GetPathCounts(m_PathAB, M, D, I);
	m_OutputLock.lock();
	fprintf(f, "%.3g", m_EvalueBA);
	fprintf(f, "\t%s", m_ChainB->m_Label.c_str());
	fprintf(f, "\t%s", m_ChainA->m_Label.c_str());
	fprintf(f, "\t%u", m_LoB + 1);
	fprintf(f, "\t%u", m_LoB + M + I);
	fprintf(f, "\t%u", m_LoA + 1);
	fprintf(f, "\t%u", m_LoA + M + D);
	fprintf(f, "\t%s\n", CIGAR.c_str());
	m_OutputLock.unlock();
	}

float DSSAligner::AlignCombo(
  const vector<byte> &LettersA, const vector<byte> &LettersB,
  uint &LoA, uint &LoB, string &Path)
	{
	m_ComboLettersA = &LettersA;
	m_ComboLettersB = &LettersB;
	SetSMx_Combo();
	uint LA = SIZE(LettersA);
	uint LB = SIZE(LettersB);
//float SWFast(XDPMem &Mem, const Mx<float> &SMx, uint LA, uint LB,
//  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
//  string &Path)
	float GapOpen = -(float) m_Params->m_ParaComboGapOpen;
	float GapExt = -(float) m_Params->m_ParaComboGapExt;
	uint Loi, Loj, Leni, Lenj;

	StartTimer(SWFast);
	float FwdScore = SWFast(m_Mem, m_SMx, LA, LB, GapOpen, GapExt,
	  Loi, Loj, Leni, Lenj, Path);
	EndTimer(SWFast);

	LoA = Loi;
	LoB = Loj;
	return FwdScore;
	}

void DSSAligner::AllocDProw(uint LB)
	{
	uint Size = 2*LB + 4;
	if (Size <= m_DProwSize)
		return;
	if (m_DProw != 0)
		myfree(m_DProw);
	m_DProwSize = 2*(LB + 100);
	m_DProw = myalloc(float, m_DProwSize);
	}

float DSSAligner::AlignComboQP(const vector<byte> &LettersA,
  const vector<byte> &LettersB)
	{
	m_ComboLettersA = &LettersA;
	m_ComboLettersB = &LettersB;
	if (m_Params->m_UsePara)
		{
		//SetComboQP_Para();
		float ScorePara = AlignComboQP_Para();
		return ScorePara;
		}

	uint LA = SIZE(LettersA);
	uint LB = SIZE(LettersB);
	//SetComboQP();
	AllocDProw(LB);
	StartTimer(SWFastGaplessProfb);
	float Scorefb = SWFastGaplessProfb(m_DProw, m_ProfCombo.data(), LA, LettersB.data(), LB);
	EndTimer(SWFastGaplessProfb);
	return Scorefb;
	}

float DSSAligner::AlignCombo_Int(const vector<byte> &LettersA,
  const vector<byte> &LettersB)
	{
	m_ComboLettersA = &LettersA;
	m_ComboLettersB = &LettersB;
	StartTimer(SetComboQPi);
	SetComboQPi();
	EndTimer(SetComboQPi);
	uint LA = SIZE(LettersA);
	uint LB = SIZE(LettersB);
	StartTimer(SWFastGapless_Profi);
	float FwdScore = (float) SWFastPinopGapless(m_ProfComboi.data(), LA,
	  (const int8_t *) LettersB.data(), LB);
	float RevScore = (float) SWFastPinopGapless(m_ProfComboRevi.data(), LA,
		(const int8_t *) LettersB.data(), LB);
	EndTimer(SWFastGapless_Profi);
	return (FwdScore - RevScore)/2;
	}

void DSSAligner::Stats()
	{
	ProgressLog("Alns %s, Ufil %.1f%%, CFil %.1f%%\n",
	  MemBytesToStr(m_AlnCount),
	  GetPct(m_UFilterCount, m_AlnCount),
	  GetPct(m_ComboFilterCount, m_AlnCount - m_UFilterCount));
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
