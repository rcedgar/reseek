#include "myutils.h"
#include "mx.h"
#include "dssaligner.h"
#include "pdbchain.h"
#include "alpha.h"
#include "xdpmem.h"
#include "timing.h"
#include "mumx.h"
#include "cigar.h"
#include <thread>
#include <set>
#include <mutex>

mutex DSSAligner::m_OutputLock;
mutex DSSAligner::m_StatsLock;

uint SWFastPinopGapless(const int8_t * const *AP, uint LA,
  const int8_t *B, uint LB);
void LogAln(const char *A, const char *B, const char *Path, unsigned ColCount);
float SWFast(XDPMem &Mem, const Mx<float> &SMx, uint LA, uint LB,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path);
float SWFastGapless(XDPMem &Mem, const Mx<float> &SMx, uint LA, uint LB,
  uint &Besti, uint &Bestj);
//float SWFastGaplessProf(XDPMem &Mem, const float * const *ProfA, uint LA,
//  const byte *B, uint LB, uint &Besti, uint &Bestj);
float SWFastGapless(XDPMem &Mem, const Mx<float> &SMx, uint LA, uint LB,
  uint &Besti, uint &Bestj);
int SWFastGapless_Int(XDPMem &Mem, const Mx<int8_t> &SMx, uint LA, uint LB,
  uint &Besti, uint &Bestj);
void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);
float SWFastGaplessProfb(float *DProw_, const float * const *ProfA,
	uint LA, const byte *B, uint LB);

uint DSSAligner::m_AlnCount;
uint DSSAligner::m_SWCount;
uint DSSAligner::m_MuFilterCount;
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

void InvertPath(const string &Path, string &InvPath)
	{
	InvPath.clear();
	const uint n = SIZE(Path);
	InvPath.reserve(n);
	for (uint i = 0; i < n; ++i)
		{
		char c = Path[i];
		if (c == 'M')
			InvPath += c;
		else if (c == 'D')
			InvPath += 'I';
		else if (c == 'I')
			InvPath += 'D';
		}
	}

void GetPathCounts(const string &Path, uint &M, uint &D, uint &I)
	{
	M = 0;
	D = 0;
	I = 0;
	for (uint i = 0; i < SIZE(Path); ++i)
		{
		char c = Path[i];
		if (c == 'M')
			++M;
		else if (c == 'D')
			++D;
		else if (c == 'I')
			++I;
		}
	}

DSSAligner::DSSAligner()
	{
	if (optset_columns)
		{
		vector<string> Fields;
		Split(string(opt_columns), Fields, '+');
		const uint n = SIZE(Fields);
		if (n == 0)
			Die("Empty -columns option");
		for (uint i = 0; i < n; ++i)
			{
			if (Fields[i] == "std")
				{
				m_UFs.push_back(UF_query);
				m_UFs.push_back(UF_target);
				m_UFs.push_back(UF_qlo);
				m_UFs.push_back(UF_qhi);
				m_UFs.push_back(UF_ql);
				m_UFs.push_back(UF_tlo);
				m_UFs.push_back(UF_thi);
				m_UFs.push_back(UF_tl);
				m_UFs.push_back(UF_pctid);
				m_UFs.push_back(UF_evalue);
				}
			else
				{
				USERFIELD UF = StrToUF(Fields[i]);
				m_UFs.push_back(UF);
				}
			}
		}
	else
		{
		m_UFs.push_back(UF_aq);
		m_UFs.push_back(UF_query);
		m_UFs.push_back(UF_target);
		m_UFs.push_back(UF_evalue);
		}
	}

int DSSAligner::GetMuDPScorePathInt(const vector<byte> &MuLettersA,
  const vector<byte> &MuLettersB, uint LoA, uint LoB,
  const string &Path) const
	{
	uint Sum = 0;
	uint PosA = LoA;
	uint PosB = LoB;
	const int Open = -m_Params->m_ParaMuGapOpen;
	const int Ext = -m_Params->m_ParaMuGapExt;
	const float FwdMatchScore = m_Params->m_FwdMatchScore;
	const uint ColCount = SIZE(Path);
	extern int8_t IntScoreMx_Mu[36][36];
 
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			{
			asserta(PosA < SIZE(MuLettersA));
			asserta(PosB < SIZE(MuLettersB));
			byte a = MuLettersA[PosA];
			byte b = MuLettersB[PosB];
			Sum += IntScoreMx_Mu[a][b];
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

float DSSAligner::GetMuDPScorePath(const vector<byte> &LettersA,
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
			extern float ScoreMx_Mu[36][36];
			Sum += ScoreMx_Mu[LetterA][LetterB];
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

// GetDPScorePath calculates AlnScore which is optimized by SWFast.
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
			{
			float ColScore = GetScorePosPair(ProfileA, ProfileB, PosA, PosB);
			Sum += ColScore;
			if (opt_tracedpscorepath)
				Log("M  %5u  %5u  %10.3g  %10.3g\n",
				  PosA, PosB, ColScore, Sum);
			++PosA;
			++PosB;
			break;
			}

		case 'D':
			if (Col != 0 && Path[Col-1] == 'D')
				{
				Sum += Ext;
				if (opt_tracedpscorepath)
					Log("De %5u  %5u  %10.3g  %10.3g\n",
					  PosA, PosB, Ext, Sum);
				}
			else
				{
				Sum += Open;
				if (opt_tracedpscorepath)
					Log("Do %5u  %5u  %10.3g  %10.3g\n",
					  PosA, PosB, Open, Sum);
				}
			++PosA;
			break;

		case 'I':
			if (Col != 0 && Path[Col-1] == 'I')
				{
				Sum += Ext;
				if (opt_tracedpscorepath)
					Log("Ie %5u  %5u  %10.3g  %10.3g\n",
					  PosA, PosB, Ext, Sum);
				}
			else
				{
				Sum += Open;
				if (opt_tracedpscorepath)
					Log("Io %5u  %5u  %10.3g  %10.3g\n",
					  PosA, PosB, Open, Sum);
				}
			++PosB;
			break;

		default:
			asserta(false);
			}
		}
	if (opt_tracedpscorepath)
		Log("Total score %.3g\n", Sum);
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
	SMx.Alloc(LA, LB);
	Mx<float> &RevSMx = m_RevSMx;
	RevSMx.Alloc(LA, LB);
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

void DSSAligner::SetSMx_Mu_Int()
	{
	const vector<byte> &MuLettersA = *m_MuLettersA;
	const vector<byte> &MuLettersB = *m_MuLettersB;
	const uint LA = SIZE(MuLettersA);
	const uint LB = SIZE(MuLettersB);
	const uint AS = 36;
	asserta(AS == 36);

	m_SMx_Int.Clear();
	m_RevSMx_Int.Clear();
	Mx<int8_t> &SMx = m_SMx_Int;
	SMx.Alloc(LA, LB);
	Mx<int8_t> &RevSMx = m_RevSMx_Int;
	RevSMx.Alloc(LA, LB);
	StartTimer(SetSMx_Mu_Int);
	int8_t **Sim = SMx.GetData();
	int8_t **RevSim = RevSMx.GetData();

	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		byte LetterA = MuLettersA[PosA];
		asserta(LetterA < AS);
		const int8_t *MuMxRow = IntScoreMx_Mu[LetterA];
		int8_t *SimRow = Sim[PosA];
		int8_t *RevSimRow = RevSim[PosA];
		for (uint PosB = 0; PosB < LB; ++PosB)
			{
			byte LetterB = MuLettersB[PosB];
			asserta(LetterB < AS);
			int8_t Score = MuMxRow[LetterB];
			SimRow[PosB] = Score;
			RevSimRow[LB-1-PosB] = Score;
			}
		}
	EndTimer(SetSMx_Mu_Int);
	}

void DSSAligner::SetMuQP()
	{
	StartTimer(SetMuQP);
	const vector<byte> &MuLettersA = *m_MuLettersA;
	uint LA = SIZE(MuLettersA);
	uint n = SIZE(m_ProfMu);
	if (LA > n)
		{
		m_ProfMu.resize(LA);
		m_ProfMuRev.resize(LA);
		}
	const uint AS = 36;
	asserta(AS == 36);
	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		byte LetterA = MuLettersA[PosA];
		asserta(LetterA < AS);
		const float *MuMxRow = ScoreMx_Mu[LetterA];
		m_ProfMu[PosA] = MuMxRow;
		m_ProfMuRev[LA-PosA-1] = MuMxRow;
		}
	EndTimer(SetMuQP);
	}

void DSSAligner::SetMuQPi()
	{
	const vector<byte> &MuLettersA = *m_MuLettersA;
	uint LA = SIZE(MuLettersA);
	uint n = SIZE(m_ProfMui);
	if (LA > n)
		{
		m_ProfMui.resize(LA);
		m_ProfMuRevi.resize(LA);
		}
	const uint AS = 36;
	asserta(AS == 36);
	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		byte LetterA = MuLettersA[PosA];
		asserta(LetterA < AS);
		const int8_t *MuMxRow = IntScoreMx_Mu[LetterA];
		m_ProfMui[PosA] = MuMxRow;
		m_ProfMuRevi[LA-PosA-1] = MuMxRow;
		}
	}

void DSSAligner::SetSMx_Mu()
	{
	const vector<byte> &MuLettersA = *m_MuLettersA;
	const vector<byte> &MuLettersB = *m_MuLettersB;
	const uint LA = SIZE(MuLettersA);
	const uint LB = SIZE(MuLettersB);
	const uint AS = 36;
	asserta(AS == 36);

	m_SMx.Clear();
	m_RevSMx.Clear();
	Mx<float> &SMx = m_SMx;
	SMx.Alloc(LA, LB);
	Mx<float> &RevSMx = m_RevSMx;
	RevSMx.Alloc(LA, LB);
	StartTimer(SetSMx_Mu);
	float **Sim = SMx.GetData();
	float **RevSim = RevSMx.GetData();

	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		byte LetterA = MuLettersA[PosA];
		asserta(LetterA < AS);
		const float *MuMxRow = ScoreMx_Mu[LetterA];
		float *SimRow = Sim[PosA];
		float *RevSimRow = RevSim[PosA];
		for (uint PosB = 0; PosB < LB; ++PosB)
			{
			byte LetterB = MuLettersB[PosB];
			asserta(LetterB < AS);
			float Score = MuMxRow[LetterB];
			SimRow[PosB] = Score;
			RevSimRow[LB-1-PosB] = Score;
			}
		}
	EndTimer(SetSMx_Mu);
	}

void DSSAligner::SetSMx_Box(int Lo_i, int Hi_i, int Lo_j, int Hi_j)
	{
	const DSSParams &Params = *m_Params;
	const vector<vector<byte> > &ProfileA = *m_ProfileA;
	const vector<vector<byte> > &ProfileB = *m_ProfileB;
	const int SegLenA = Hi_i - Lo_i + 1;
	const int SegLenB = Hi_j - Lo_j + 1;
	asserta(SegLenA > 0 && SegLenB > 0);

// Memory blows up with grow-only strategy due to tail of long chains
	m_SMx.Clear();
	Mx<float> &SMx = m_SMx;
	SMx.Alloc((uint) SegLenA, (uint) SegLenB);
	StartTimer(SetSMx_Box);
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
	//for (uint PosA = 0; PosA < LA; ++PosA)
	for (int kA = 0; kA < SegLenA; ++kA)
		{
		uint PosA = uint(Lo_i + kA);
		byte ia = ProfRowA[PosA];
		assert(ia < AlphaSize0);
		const float *ScoreMxRow = ScoreMx0[ia];

		float *SimRow = Sim[kA];
		for (int kB = 0; kB < SegLenB; ++kB)
			{
			uint PosB = uint(Lo_j + kB);
			byte ib = ProfRowB[PosB];
			assert(ia < AlphaSize0 && ib < AlphaSize0);
			float MatchScore = ScoreMxRow[ib];
			SimRow[kB] = MatchScore;
			}
		}

	for (uint FeatureIdx = 1; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = m_Params->m_Features[FeatureIdx];
		uint AlphaSize = g_AlphaSizes2[F];
		float **ScoreMx = m_Params->m_ScoreMxs[F];
		const vector<byte> &ProfRowA = ProfileA[FeatureIdx];
		const vector<byte> &ProfRowB = ProfileB[FeatureIdx];
		for (int kA = 0; kA < SegLenA; ++kA)
			{
			uint PosA = uint(Lo_i + kA);
			byte ia = ProfRowA[PosA];
			assert(ia < AlphaSize);
			const float *ScoreMxRow = ScoreMx[ia];

			float *SimRow = Sim[kA];
			for (int kB = 0; kB < SegLenB; ++kB)
				{
				uint PosB = uint(Lo_j + kB);
				byte ib = ProfRowB[PosB];
				assert(ib < AlphaSize);
				float MatchScore = ScoreMxRow[ib];
				SimRow[kB] += MatchScore;
				}
			}
		}
	EndTimer(SetSMx_Box);
	}

void DSSAligner::SetSMx_NoRev()
	{
	const DSSParams &Params = *m_Params;
	const vector<vector<byte> > &ProfileA = *m_ProfileA;
	const vector<vector<byte> > &ProfileB = *m_ProfileB;
	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

// Memory blows up with grow-only strategy due to tail of long chains
	m_SMx.Clear();
	Mx<float> &SMx = m_SMx;
	SMx.Alloc(LA, LB);
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

float DSSAligner::GetMuScore()
	{
	float MuScore = AlignMuQP(*m_MuLettersA, *m_MuLettersB);
	return MuScore;
	}

bool DSSAligner::MuFilter()
	{
	if (m_MuLettersA == 0 || m_MuLettersB == 0)
		return true;
	float MCS = m_Params->m_Omega;
	if (MCS <= 0)
		return true;
	float MuScore = GetMuScore(); // AlignMuQP(*m_MuLettersA, *m_MuLettersB);
	if (MuScore < MCS)
		return false;
	return true;
	}

bool DSSAligner::UFilter()
	{
	uint MinU = uint(round(m_Params->m_MinU));
	if (MinU == 0)
		return true;
	if (m_MuKmerBitsA == 0 || m_MuKmerBitsB == 0)
		return true;
	uint U = GetUBits(*m_MuKmerBitsA, *m_MuKmerBitsB);
	if (U < MinU)
		return false;
	return true;
	}

void DSSAligner::Align_MuFilter(
  const PDBChain &ChainA, const PDBChain &ChainB,
  const vector<byte> &MuLettersA, const vector<uint> &MuKmersA,
  const vector<byte> &MuLettersB,const vector<uint> &MuKmersB,
  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB)
	{
	SetQuery(ChainA, &ProfileA, 0, &MuLettersA, &MuKmersA, FLT_MAX);
	SetTarget(ChainB, &ProfileB, 0, &MuLettersB, &MuKmersB, FLT_MAX);

	//m_EvalueA = FLT_MAX;
	//m_EvalueB = FLT_MAX;
	//m_Path.clear();
	ClearAlign();

	m_StatsLock.lock();
	++m_AlnCount;
	m_StatsLock.unlock();

	bool MuFilterOk = MuFilter();
	if (!MuFilterOk)
		{
		m_StatsLock.lock();
		++m_MuFilterCount;
		m_StatsLock.unlock();
		return;
		}
	Align_NoAccel();
	}

void DSSAligner::SetParams(const DSSParams &Params)
	{
	m_Params = &Params;
	m_MKF.m_Params = &Params;
	}

void DSSAligner::SetQuery(
	const PDBChain &Chain,
	const vector<vector<byte> > *ptrProfile,
	const vector<uint> *ptrMuKmerBits,
	const vector<byte> *ptrMuLetters,
	const vector<uint> *ptrMuKmers,
	float SelfRevScore)
	{
	m_ChainA = &Chain;
	m_ProfileA = ptrProfile;
	m_MuKmerBitsA = ptrMuKmerBits;
	m_MuLettersA = ptrMuLetters;
	m_MuKmersA = ptrMuKmers;
	m_SelfRevScoreA = SelfRevScore;
	if (ptrMuLetters != 0 && m_Params->m_Omega > 0)
		{
		if (m_Params->m_UsePara)
			SetMuQP_Para();
		else
			SetMuQP();
		}

	if (ptrMuKmers != 0)
		{
		asserta(ptrMuLetters != 0);
		m_MKF.MuKmerResetQ();
		m_MKF.MuKmerSetQ(Chain, ptrMuLetters, ptrMuKmers);
		}
	}

void DSSAligner::SetTarget(
	const PDBChain &Chain,
	const vector<vector<byte> > *ptrProfile,
	const vector<uint> *ptrMuKmerBits,
	const vector<byte> *ptrMuLetters,
	const vector<uint> *ptrMuKmers,
	float SelfRevScore)
	{
	m_ChainB = &Chain;
	m_ProfileB = ptrProfile;
	m_MuKmerBitsB = ptrMuKmerBits;
	m_MuKmersB = ptrMuKmers;
	m_MuLettersB = ptrMuLetters;
	m_SelfRevScoreB = SelfRevScore;
	}

void DSSAligner::AlignMuOnly()
	{
	//m_EvalueA = FLT_MAX;
	//m_EvalueB = FLT_MAX;
	//m_Path.clear();
	ClearAlign();

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
	bool MuFilterOk = MuFilter();
	if (!MuFilterOk)
		{
		m_StatsLock.lock();
		++m_MuFilterCount;
		m_StatsLock.unlock();
		return;
		}
	AlignMuPath();
	}

void DSSAligner::AlignQueryTarget()
	{
	//m_EvalueA = FLT_MAX;
	//m_EvalueB = FLT_MAX;
	//m_Path.clear();
	ClearAlign();

	m_StatsLock.lock();
	++m_AlnCount;
	m_StatsLock.unlock();

	if (m_Params->m_MinU > 0)
		{
		bool UFilterOk = UFilter();
		if (!UFilterOk)
			{
			m_StatsLock.lock();
			++m_UFilterCount;
			m_StatsLock.unlock();
			return;
			}
		}

	if (m_Params->m_Omega > 0)
		{
		bool MuFilterOk = MuFilter();
		if (!MuFilterOk)
			{
			m_StatsLock.lock();
			++m_MuFilterCount;
			m_StatsLock.unlock();
			return;
			}
		}

	Align_NoAccel();
	}

void DSSAligner::CalcEvalue_AAOnly()
	{
	static const float Log2 = logf(2);
	static const float GappedLambda = 0.267f;
	static const float LogGappedK = logf(0.0410f);

	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();
	const float DBSize = m_Params->m_DBSize;

	float Score = m_AlnFwdScore;
	float BitScore = (Score*GappedLambda - LogGappedK)/Log2;
	float NM_A = float(LA)*float(DBSize);
	float NM_B = float(LB)*float(DBSize);
	m_QualityA = 0;
	m_QualityB = 0;
	m_EvalueA = NM_A/powf(2, BitScore);
	m_EvalueB = NM_B/powf(2, BitScore);
	}

void DSSAligner::CalcEvalue()
	{
	if (m_Params->m_AAOnly)
		{
		CalcEvalue_AAOnly();
		return;
		}

// Threshold enables small speedup by avoiding LDDT and self-rev
	if (m_AlnFwdScore < m_Params->m_MinFwdScore)
		return;

	StartTimer(CalcEvalue)
	float LDDT = GetLDDT();
	float RevDPScore = 0;
	asserta(m_SelfRevScoreA != -FLT_MAX);
	asserta(m_SelfRevScoreB != -FLT_MAX);
	RevDPScore = 0;
	if (m_SelfRevScoreA != FLT_MAX && m_SelfRevScoreB != FLT_MAX)
		RevDPScore = (m_SelfRevScoreA + m_SelfRevScoreB)/2;
	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();
	float L = float(LA + LB)/2;

	const float dpw = 1.7f;
	const float lddtw = 0.13f;
	const float ladd = 250.0f;
	const float revtsw = 2.0f;

	m_NewTestStatisticA = lddtw*LDDT;
	m_NewTestStatisticA += (dpw*m_AlnFwdScore - revtsw*RevDPScore)/(L + ladd);

	m_NewTestStatisticB = m_NewTestStatisticA;

	uint M, D, I;
	GetPathCounts(m_Path, M, D, I);
	m_HiA = m_LoA + M + D - 1;
	m_HiB = m_LoB + M + I - 1;
	m_Ids = M;
	m_Gaps = D + I;

	const float a = 5.0f;
	const float b = -40.0f;
	float logE = a + b*m_NewTestStatisticA;
	float DBSize = m_Params->m_DBSize;
//a = 5
//b = -40
//def get_qual(ts):
//	LogE = a + b*ts
//	if LogE < -20:
//		return 1
//	c = LogE
//	# x = math.exp(0.25*c)
//	x = 10**(c/10)
//	Qual = 1/(1 + x)
//	return Qual

	float Qual = 0;
	if (logE < -20)
		Qual = 1;
	else
		{
		float x = powf(10, logE/10);
		Qual = 1/(1 + x/2);
		}

	float E_scop = expf(logE)/11211;
	float E = E_scop*DBSize;
	m_QualityA = Qual;
	m_QualityB = Qual;
	m_EvalueA = E;
	m_EvalueB = E;
	EndTimer(CalcEvalue)
	}

void DSSAligner::ClearAlign()
	{
	m_Path.clear();
	m_LoA = UINT_MAX;
	m_LoB = UINT_MAX;
	m_HiA = UINT_MAX;
	m_HiB = UINT_MAX;
	m_Ids = UINT_MAX;
	m_Gaps = UINT_MAX;
	m_EvalueA = FLT_MAX;
	m_EvalueB = FLT_MAX;
	m_TestStatisticA = -FLT_MAX;
	m_TestStatisticB = -FLT_MAX;
	m_NewTestStatisticA = -FLT_MAX;
	m_NewTestStatisticB = -FLT_MAX;
	m_AlnFwdScore = -FLT_MAX;
	}

void DSSAligner::Align_Box(int Lo_i, int Hi_i, int Lo_j, int Hi_j)
	{
	ClearAlign();

	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	Lo_i -= 64;
	Hi_i += 64;
	Lo_j -= 64;
	Hi_j += 64;

	if (Lo_i < 0) Lo_i = 0;
	if (Lo_j < 0) Lo_j = 0;

	if (Hi_i >= int(LA)) Hi_i = int(LA) - 1;
	if (Hi_j >= int(LB)) Hi_j = int(LB) - 1;

	//SetSMx_NoRev();
	SetSMx_Box(Lo_i, Hi_i, Lo_j, Hi_j);

	const uint SegLenA = uint(Hi_i - Lo_i + 1);
	const uint SegLenB = uint(Hi_j - Lo_j + 1);

	asserta(SegLenA <= LA && SegLenB <= LB);

	StartTimer(SW_Box);
	uint LoA, LoB, Leni, Lenj;
	m_AlnFwdScore = SWFast(m_Mem, m_SMx, SegLenA, SegLenB,
	  m_Params->m_GapOpen, m_Params->m_GapExt,
	  LoA, LoB, Leni, Lenj, m_Path);
	m_LoA = Lo_i + LoA;
	m_LoB = Lo_j + LoB;
	EndTimer(SW_Box);

	CalcEvalue();
	}

void DSSAligner::Align_NoAccel()
	{
	ClearAlign();
	SetSMx_NoRev();

	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	StartTimer(SWFwd);
	uint Leni, Lenj;
	m_AlnFwdScore = SWFast(m_Mem, m_SMx, LA, LB,
	  m_Params->m_GapOpen, m_Params->m_GapExt,
	  m_LoA, m_LoB, Leni, Lenj, m_Path);
	EndTimer(SWFwd);

	CalcEvalue();
	}

void DSSAligner::AlignMuPath()
	{
	//m_Path.clear();
	//m_LoA = UINT_MAX;
	//m_LoB = UINT_MAX;
	//m_AlnFwdScore = 0;
	ClearAlign();
	float MuFwdScore = AlignMuQP_Para_Path(m_LoA, m_LoB, m_Path);
	if (MuFwdScore < 0)
		{
		m_StatsLock.lock();
		++m_ParasailSaturateCount;
		m_StatsLock.unlock();
		SetSMx_NoRev();
		Mx<float> &SMx = m_SMx;
		Mx<float> &RevSMx = m_RevSMx;
		MuFwdScore = 
		  AlignMu(*m_MuLettersA, *m_MuLettersB, m_LoA, m_LoB, m_Path);
		m_AlnFwdScore =
		  GetDPScorePath(*m_ProfileA, *m_ProfileB, m_LoA, m_LoB, m_Path);
		}
	CalcEvalue();
	}

void DSSAligner::ToAln(FILE *f, bool Up) const
	{
	if (f == 0)
		return;
	if (Up)
		PrettyAln(f, *m_ChainA, *m_ChainB, *m_ProfileA, *m_ProfileB,
		  m_LoA, m_LoB, m_Path, m_QualityA, m_EvalueA);
	else
		{
		string Path;
		InvertPath(m_Path, Path);
		PrettyAln(f, *m_ChainB, *m_ChainA, *m_ProfileB, *m_ProfileA,
		  m_LoB, m_LoA, Path, m_QualityB, m_EvalueB);
		}
	}

void DSSAligner::ToFasta2(FILE *f, bool Global, bool Up) const
	{
	if (f == 0)
		return;

	string RowA, RowB;
	if (Up)
		{
		GetRow_A(RowA, Global);
		GetRow_B(RowB, Global);
		}
	else
		{
		GetRow_B(RowA, Global);
		GetRow_A(RowB, Global);
		}

	const string &LabelA = GetLabel(Up);
	const string &LabelB = GetLabel(!Up);
	float Evalue = GetEvalue(Up);
	float PctId = GetPctId();
	string LabelAx = LabelA;
	Psa(LabelAx, " E=%.3g Id=%.1f%%", Evalue, PctId);
	LabelAx += " (";
	LabelAx += LabelB;
	LabelAx += ")";

	m_OutputLock.lock();
	SeqToFasta(f, LabelAx.c_str(), RowA);
	SeqToFasta(f, LabelB, RowB);
	fputc('\n', f);
	m_OutputLock.unlock();
	}

void DSSAligner::ToTsv(FILE *f, bool Up)
	{
	if (f == 0)
		return;

	m_OutputLock.lock();
	const uint n = SIZE(m_UFs);
	asserta(n > 0);
	for (uint i = 0; i < n; ++i)
		{
		if (i > 0)
			fputc('\t', f);
		WriteUserField(f, m_UFs[i], Up);
		}
	fputc('\n', f);
	m_OutputLock.unlock();
	}

float DSSAligner::AlignMu(
  const vector<byte> &LettersA, const vector<byte> &LettersB,
  uint &LoA, uint &LoB, string &Path)
	{
	m_MuLettersA = &LettersA;
	m_MuLettersB = &LettersB;
	SetSMx_Mu();
	uint LA = SIZE(LettersA);
	uint LB = SIZE(LettersB);

	float GapOpen = -(float) m_Params->m_ParaMuGapOpen;
	float GapExt = -(float) m_Params->m_ParaMuGapExt;
	uint Loi, Loj, Leni, Lenj;

	StartTimer(SWFast);
	float FwdScore = SWFast(m_Mem, m_SMx, LA, LB, GapOpen, GapExt,
	  Loi, Loj, Leni, Lenj, Path);
	EndTimer(SWFast);
	m_Mem.Clear();
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

float DSSAligner::AlignMuQP(const vector<byte> &LettersA,
  const vector<byte> &LettersB)
	{
	m_MuLettersA = &LettersA;
	m_MuLettersB = &LettersB;
	if (m_Params->m_UsePara)
		{
		//SetMuQP_Para();
		float ScorePara = AlignMuQP_Para();
		return ScorePara;
		}

	uint LA = SIZE(LettersA);
	uint LB = SIZE(LettersB);
	//SetMuQP();
	AllocDProw(LB);
	StartTimer(SWFastGaplessProfb);
	float Scorefb = SWFastGaplessProfb(m_DProw, m_ProfMu.data(), LA, LettersB.data(), LB);
	EndTimer(SWFastGaplessProfb);
	return Scorefb;
	}

float DSSAligner::AlignMu_Int(const vector<byte> &LettersA,
  const vector<byte> &LettersB)
	{
	m_MuLettersA = &LettersA;
	m_MuLettersB = &LettersB;
	StartTimer(SetMuQPi);
	SetMuQPi();
	EndTimer(SetMuQPi);
	uint LA = SIZE(LettersA);
	uint LB = SIZE(LettersB);
	StartTimer(SWFastGapless_Profi);
	float FwdScore = (float) SWFastPinopGapless(m_ProfMui.data(), LA,
	  (const int8_t *) LettersB.data(), LB);
	float RevScore = (float) SWFastPinopGapless(m_ProfMuRevi.data(), LA,
		(const int8_t *) LettersB.data(), LB);
	EndTimer(SWFastGapless_Profi);
	return (FwdScore - RevScore)/2;
	}

void DSSAligner::Stats()
	{
	uint MuFilterInputCount = m_AlnCount - DSSAligner::m_UFilterCount;
	Log("DSSAligner::Stats() alns %s, ufil %.1f%%, mufil %.1f%% (sat %u)\n",
	  FloatToStr(m_AlnCount),
	  GetPct(m_UFilterCount, m_AlnCount),
	  GetPct(m_MuFilterCount, MuFilterInputCount),
	  m_ParasailSaturateCount);
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

void DSSAligner::GetRow(bool Up, bool Top, bool Global, string &Row) const
	{
	if (Up)
		{
		if (Top)
			GetRow_A(Row, Global);
		else
			GetRow_B(Row, Global);
		}
	else
		{
		if (Top)
			GetRow_B(Row, Global);
		else
			GetRow_A(Row, Global);
		}
	}

void DSSAligner::GetRow_A(string &Row, bool Global) const
	{
	Row.clear();
	const string &SeqA = m_ChainA->m_Seq;
	const string &SeqB = m_ChainB->m_Seq;
	const uint LA = SIZE(SeqA);
	const uint LB = SIZE(SeqB);
	const uint ColCount = SIZE(m_Path);
	if (Global)
		{
		for (uint i = m_LoA; i < m_LoB; ++i)
			Row += '.';
		for (uint i = 0; i < m_LoA; ++i)
			Row += tolower(SeqA[i]);
		}
	uint PosA = m_LoA;
	uint PosB = m_LoB;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = m_Path[Col];
		switch (c)
			{
		case 'M':
			{
			asserta(PosA <= m_HiA);
			asserta(PosB <= m_HiB);
			Row += SeqA[PosA++];
			++PosB;
			break;
			}

		case 'D':
			{
			asserta(PosA <= m_HiA);
			Row += SeqA[PosA++];
			break;
			}

		case 'I':
			Row += '-';
			++PosB;
			break;

		default:
			asserta(false);
			}
		}
	if (Global)
		{
		while (PosA < LA)
			{
			Row += tolower(SeqA[PosA++]);
			++PosB;
			}
		while (PosB++ < LB)
			Row += '.';
		}
	}

void DSSAligner::GetRow_B(string &Row, bool Global) const
	{
	Row.clear();
	const string &SeqA = m_ChainA->m_Seq;
	const string &SeqB = m_ChainB->m_Seq;
	const uint LA = SIZE(SeqA);
	const uint LB = SIZE(SeqB);
	const uint ColCount = SIZE(m_Path);
	if (Global)
		{
		for (uint i = m_LoB; i < m_LoA; ++i)
			Row += '.';
		for (uint i = 0; i < m_LoB; ++i)
			Row += tolower(SeqB[i]);
		}
	uint PosA = m_LoA;
	uint PosB = m_LoB;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = m_Path[Col];
		switch (c)
			{
		case 'M':
			{
			asserta(PosA <= m_HiA);
			asserta(PosB <= m_HiB);
			++PosA;
			Row += SeqB[PosB++];
			break;
			}

		case 'D':
			{
			++PosA;
			asserta(PosA <= m_HiA);
			Row += '-';
			break;
			}

		case 'I':
			{
			asserta(PosB <= m_HiB);
			Row += SeqB[PosB++];
			break;
			}

		default:
			asserta(false);
			}
		}
	if (Global)
		{
		while (PosB < LB)
			{
			Row += tolower(SeqB[PosB++]);
			++PosA;
			}
		while (PosA++ < LA)
			Row += '.';
		}
	}

void DSSAligner::GetPosABs(vector<uint> &PosAs,
  vector<uint> &PosBs) const
	{
	PosAs.clear();
	PosBs.clear();
	const uint Cols = SIZE(m_Path);
	uint PosA = m_LoA;
	uint PosB = m_LoB;
	for (uint Col = 0; Col < Cols; ++Col)
		{
		char c = m_Path[Col];
		switch (c)
			{
		case 'M':
			PosAs.push_back(PosA);
			PosBs.push_back(PosB);
			++PosA;
			++PosB;
			break;

		case 'D':
			++PosA;
			break;

		case 'I':
			++PosB;
			break;
			}
		}
	}

float DSSAligner::GetLDDT() const
	{
	double GetLDDT_mu(const PDBChain &Q, const PDBChain &T,
	  const vector<uint> &PosQs, const vector<uint> &PosTs,
	  bool DaliScorerCompatible);
	vector<uint> PosAs;
	vector<uint> PosBs;
	GetPosABs(PosAs, PosBs);
	double LDDT = 
	  GetLDDT_mu(*m_ChainA, *m_ChainB, PosAs, PosBs, false);
	return (float) LDDT;
	}

float DSSAligner::GetPctId() const
	{
	const string &Path = m_Path;
	const uint ColCount = SIZE(Path);
	uint PosA = m_LoA;
	uint PosB = m_LoB;
	const string &SeqA = m_ChainA->m_Seq;
	const string &SeqB = m_ChainB->m_Seq;
	const uint LA = SIZE(SeqA);
	const uint LB = SIZE(SeqB);
	uint N = 0;
	uint n = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			{
			asserta(PosA < LA);
			asserta(PosB < LB);
			char a = SeqA[PosA];
			char b = SeqB[PosB];
			++PosA;
			++PosB;
			++N;
			if (a == b)
				++n;
			break;
			}

		case 'D':
			++PosA;
			break;

		case 'I':
			++PosB;
			break;

		default:
			asserta(false);
			}
		}
	return N == 0 ? 0 : (n*100.0f)/N;
	}

float DSSAligner::GetKabsch(double t[3], double u[3][3], bool Up) const
	{
	double Kabsch(const PDBChain &ChainA, const PDBChain &ChainB,
		uint LoA, uint LoB, const string &Path,
		double t[3], double u[3][3]);

	if (Up)
		return (float) Kabsch(*m_ChainA, *m_ChainB, m_LoA, m_LoB, m_Path, t, u);
	else
		{
		string Path;
		InvertPath(m_Path, Path);
		return (float) Kabsch(*m_ChainB, *m_ChainA, m_LoB, m_LoA, Path, t, u);
		}
	}

void DSSAligner::AlignMKF()
	{
	ClearAlign();
	bool FoundHSP = m_MKF.MuKmerAln(*m_ChainA, *m_MuLettersB, *m_MuKmersB);

	if (FoundHSP && m_MKF.m_BestChainScore > 0)
		{
		SetSMx_Box(m_MKF.m_ChainLo_i, m_MKF.m_ChainHi_i, m_MKF.m_ChainLo_j, m_MKF.m_ChainHi_j);
		Align_Box(m_MKF.m_ChainLo_i, m_MKF.m_ChainHi_i, m_MKF.m_ChainLo_j, m_MKF.m_ChainHi_j);
		}
	}
