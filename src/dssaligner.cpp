#include "myutils.h"
#include "mx.h"
#include "dssaligner.h"
#include "pdbchain.h"
#include "alpha.h"
#include "xdpmem.h"
#include "timing.h"
#include "mumx.h"
#include "cigar.h"
#include "parasail.h"
#include "kabsch.h"
#include <thread>
#include <set>
#include <mutex>

mutex DSSAligner::m_OutputLock;

//uint SWFastPinopGapless(const int8_t * const *AP, uint LA,
//  const int8_t *B, uint LB);
float SWFast(XDPMem &Mem, const float * const *SMxData, uint LA, uint LB,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path);
void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);
float SWFastGaplessProfb(float *DProw_, const float * const *ProfA,
	uint LA, const byte *B, uint LB);

atomic<uint> DSSAligner::m_AlnCount;
atomic<uint> DSSAligner::m_XDropAlnCount;
atomic<uint> DSSAligner::m_XDropDiscardCount;
atomic<uint> DSSAligner::m_SWCount;
atomic<uint> DSSAligner::m_MuFilterInputCount;
atomic<uint> DSSAligner::m_MuFilterDiscardCount;
atomic<uint> DSSAligner::m_ParasailSaturateCount;

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

DSSAligner::~DSSAligner()
	{
	if (m_ProfPara != 0)
		parasail_profile_free((parasail_profile_t *) m_ProfPara);
	if (m_ProfParaRev != 0)
		parasail_profile_free((parasail_profile_t *) m_ProfParaRev);
	}

DSSAligner::DSSAligner()
	{
	if (optset_columns)
		{
		vector<string> Fields;
		Split(string(opt(columns)), Fields, '+');
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

float DSSAligner::GetDPScoreGivenPath(const vector<vector<byte> > &Profile1,
	const vector<vector<byte> > &Profile2, const string &Path) const
	{
	const uint ColCount = SIZE(Path);
	const uint L1 = SIZE(Profile1[0]);
	const uint L2 = SIZE(Profile2[0]);
	uint Pos1 = 0;
	uint Pos2 = 0;
	bool InGap = false;
	float GapOpen = m_Params->m_GapOpen;
	float GapExt = m_Params->m_GapExt;
	asserta(GapOpen < 0);
	asserta(GapExt < 0);
	float Score = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			Score += GetScorePosPair(Profile1, Profile2, Pos1++, Pos2++);
			InGap = false;
			break;

		case 'D':
			if (InGap)
				Score += GapExt;
			else
				{
				Score += GapOpen;
				InGap = true;
				}
			++Pos1;
			break;

		case 'I':
			if (InGap)
				Score += GapExt;
			else
				{
				Score += GapOpen;
				InGap = true;
				}
			++Pos2;
			break;

		default:
			asserta(false);
			}
		}
	asserta(Pos1 == L1);
	asserta(Pos2 == L2);
	return Score;
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
		uint AlphaSize = DSS::GetAlphaSize(F); // g_AlphaSizes2[F];
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

float DSSAligner::GetMegaHSPScore(uint Lo_i, uint Lo_j, uint Len)
	{
	StartTimer(GetMegaHSPScore);
	const DSSParams &Params = *m_Params;
	const vector<vector<byte> > &ProfileA = *m_ProfileA;
	const vector<vector<byte> > &ProfileB = *m_ProfileB;
	const uint FeatureCount = Params.GetFeatureCount();
	asserta(SIZE(ProfileA) == FeatureCount);
	asserta(SIZE(ProfileB) == FeatureCount);

// Special case first feature because = not += and
	FEATURE F0 = m_Params->m_Features[0];
	uint AlphaSize0 = DSS::GetAlphaSize(F0); // g_AlphaSizes2[F0];
	float **ScoreMx0 = m_Params->m_ScoreMxs[F0];
	const vector<byte> &ProfRowA = ProfileA[0];
	const vector<byte> &ProfRowB = ProfileB[0];
	float Total = 0;
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = m_Params->m_Features[FeatureIdx];
		uint AlphaSize = DSS::GetAlphaSize(F); // g_AlphaSizes2[F];
		float **ScoreMx = m_Params->m_ScoreMxs[F];
		const vector<byte> &ProfRowA = ProfileA[FeatureIdx];
		const vector<byte> &ProfRowB = ProfileB[FeatureIdx];
		for (uint k = 0; k < Len; ++k)
			{
			uint PosA = uint(Lo_i + k);
			byte ia = ProfRowA[PosA];
			assert(ia < AlphaSize);
			const float *ScoreMxRow = ScoreMx[ia];

			uint PosB = uint(Lo_j + k);
			byte ib = ProfRowB[PosB];
			assert(ib < AlphaSize);
			Total += ScoreMxRow[ib];
			}
		}
	EndTimer(GetMegaHSPScore);
	return Total;
	}

void DSSAligner::SetSMx_NoRev(const DSSParams &Params,
					  const vector<vector<byte> > &ProfileA,
					  const vector<vector<byte> > &ProfileB)
	{
	//const DSSParams &Params = *m_Params;
	//const vector<vector<byte> > &ProfileA = *m_ProfileA;
	//const vector<vector<byte> > &ProfileB = *m_ProfileB;
	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

// Memory blows up with grow-only strategy due to tail of long chains
	//m_SMx.Clear();
	//Mx<float> &SMx = m_SMx;
	//m_SMx.Alloc(LA, LB, __FILE__, __LINE__);
	AllocSMxData(LA, LB);
	StartTimer(SetSMx_NoRev);
	float **Sim = GetSMxData();

	const uint FeatureCount = Params.GetFeatureCount();
	asserta(SIZE(ProfileA) == FeatureCount);
	asserta(SIZE(ProfileB) == FeatureCount);

// Special case first feature because = not += and
	FEATURE F0 = m_Params->m_Features[0];
	uint AlphaSize0 = DSS::GetAlphaSize(F0); // g_AlphaSizes2[F0];
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
		uint AlphaSize = DSS::GetAlphaSize(F); // g_AlphaSizes2[F];
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

float DSSAligner::GetMuScore()
	{
	float MuScore = AlignMuQP(*m_MuLettersA, *m_MuLettersB);
	return MuScore;
	}

bool DSSAligner::MuDPFilter()
	{
	if (m_MuLettersA == 0 || m_MuLettersB == 0)
		return true;
	float MCS = m_Params->m_Omega;
	if (MCS <= 0)
		return true;
	float MuScore = GetMuScore();
	if (MuScore < MCS)
		return false;
	return true;
	}

void DSSAligner::SetParams(const DSSParams &Params)
	{
	m_Params = &Params;
	m_MKF.SetParams(Params);
	}

void DSSAligner::UnsetQuery()
	{
	m_MKF.ResetQ();
	m_ChainA = 0;
	m_ProfileA = 0;
	m_MuLettersA = 0;
	m_MuKmersA = 0;
	m_SelfRevScoreA = 0;
	}

void DSSAligner::SetQuery(
	const PDBChain &Chain,
	const vector<vector<byte> > *ptrProfile,
	const vector<byte> *ptrMuLetters,
	const vector<uint> *ptrMuKmers,
	float SelfRevScore)
	{
	if (ptrMuKmers != 0)
		{
		asserta(ptrMuLetters != 0);
		m_MKF.SetQ(Chain.m_Label, ptrMuLetters, ptrMuKmers);
		}

	m_ChainA = &Chain;
	m_ProfileA = ptrProfile;
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
	}

void DSSAligner::SetTarget(
	const PDBChain &Chain,
	const vector<vector<byte> > *ptrProfile,
	const vector<byte> *ptrMuLetters,
	const vector<uint> *ptrMuKmers,
	float SelfRevScore)
	{
	m_ChainB = &Chain;
	m_ProfileB = ptrProfile;
	m_MuKmersB = ptrMuKmers;
	m_MuLettersB = ptrMuLetters;
	m_SelfRevScoreB = SelfRevScore;
	}

bool DSSAligner::DoMKF() const
	{
	if (m_MuLettersA == 0 || m_MuLettersB == 0)
		return false;
	if (m_MuKmersA == 0 || m_MuKmersB == 0)
		return false;
	uint LA = m_ChainA->GetSeqLength();
	uint LB = m_ChainB->GetSeqLength();
	if (LA >= m_Params->m_MKFL)
		return true;
	if (LB >= m_Params->m_MKFL)
		return true;
	return false;
	}

void DSSAligner::AlignQueryTarget_Trace()
	{
	Log("\n");
	Log("______________________________________\n");
	Log("A>%s(%u)\n", m_ChainA->m_Label.c_str(), m_ChainA->GetSeqLength());
	Log("B>%s(%u)\n", m_ChainB->m_Label.c_str(), m_ChainB->GetSeqLength());

	ClearAlign();

	if (DoMKF())
		{
		Log("DoMKF()=true\n");
		AlignMKF();
		Log("m_MKF.m_BestChainScore=%d\n", m_MKF.m_BestChainScore);
		Log("m_XDropScore=%.1f\n", m_XDropScore);
		Log("AlnFwdScore=%.3g\n", m_AlnFwdScore);
		float E = m_EvalueA;
		if (E > 1e5)
			Log("EvalueA=%.3g\n", E);
		else
			Log("EvalueA=%.1f\n", E);
		Log("Path=(%u)%.10s...\n", SIZE(m_Path), m_Path.c_str());
		return;
		}

	if (m_Params->m_Omega > 0)
		{
		Log("Omega > 0\n");
		++m_MuFilterInputCount;
		bool MuFilterOk = MuDPFilter();
		Log("MuFilterOk=%c\n", tof(MuFilterOk));
		if (!MuFilterOk)
			{
			++m_MuFilterDiscardCount;
			return;
			}
		}

	SetSMx_NoRev(*m_Params, *m_ProfileA, *m_ProfileB);

	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	uint Leni, Lenj;
	m_AlnFwdScore = SWFast(m_Mem, GetSMxData(), LA, LB,
	  m_Params->m_GapOpen, m_Params->m_GapExt,
	  m_LoA, m_LoB, Leni, Lenj, m_Path);

	CalcEvalue();

	Log("AlnFwdScore=%.3g\n", m_AlnFwdScore);
	float E = m_EvalueA;
	if (E > 1e5)
		Log("EvalueA=%.3g\n", E);
	else
		Log("EvalueA=%.1f\n", E);
	Log("Path=(%u)%.10s...\n", SIZE(m_Path), m_Path.c_str());
	}

// RCE 'D'/'I' convention is reverse of CIGAR
static char Fix(char c)
	{
	if (c == 'M')
		return 'M';
	else if (c == 'D')
		return 'I';
	else if (c == 'I')
		return 'D';
	asserta(false);
	return 0;
	}

void DSSAligner::GetCIGAR(string &CIGAR) const
	{
	CIGAR.clear();
	const uint ColCount = SIZE(m_Path);
	if (ColCount == 0)
		return;

	char LastC = m_Path[0];
	uint n = 1;

	if (m_LoA > 0)
		Psa(CIGAR, "%uS", m_LoA);
	if (m_LoB > 0)
		Psa(CIGAR, "%uT", m_LoB);

	for (uint i = 1; ; ++i)
		{
		char c = m_Path[i];
		if (c == 0)
			break;

		if (c == LastC)
			{
			++n;
			continue;
			}
		else
			{
			assert(n > 0);
			Psa(CIGAR, "%u%c", n, Fix(LastC));
			LastC = c;
			n = 1;
			}
		}
	if (n > 0)
		Psa(CIGAR, "%u%c", n, Fix(LastC));
	else
		asserta(false);

	void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);
	uint nm, nd, ni;
	GetPathCounts(m_Path, nm, nd, ni);
	uint na = m_LoA + nm + nd;
	uint nb = m_LoB + nm + ni;
	uint LA = m_ChainA->GetSeqLength();
	uint LB = m_ChainB->GetSeqLength();
	asserta(na <= LA);
	asserta(nb <= LB);
	if (na < LA)
		Psa(CIGAR, "%uS", LA - na);
	if (nb < LB)
		Psa(CIGAR, "%uT", LB - nb);
	}

void DSSAligner::ValidatePath() const
	{
	if (m_Path.empty())
		return;
	uint M, D, I;
	GetPathCounts(m_Path, M, D, I);
	uint LA = m_ChainA->GetSeqLength();
	uint LB = m_ChainB->GetSeqLength();
	asserta(m_LoA + M + D <= LA);
	asserta(m_LoB + M + I <= LB);
	}

void DSSAligner::AlignQueryTarget()
	{
	if (optset_label1 && optset_label2)
		{
		const string Label1 = string(opt(label1));
		const string Label2 = string(opt(label2));
		bool DoTrace = (m_ChainA->m_Label == Label1 && m_ChainB->m_Label == Label2) ||
			(m_ChainA->m_Label == Label2 && m_ChainB->m_Label == Label1);
		if (DoTrace)
			{
			ProgressLog("Trace %s, %s\n", Label1.c_str(), Label2.c_str());
			AlignQueryTarget_Trace();
			return;
			}
		}

	ClearAlign();

	if (DoMKF())
		{
		AlignMKF();
		return;
		}

	++m_AlnCount;

	if (m_Params->m_Omega > 0)
		{
		++m_MuFilterInputCount;
		bool MuFilterOk = MuDPFilter();
		if (!MuFilterOk)
			{
			++m_MuFilterDiscardCount;
			return;
			}
		}

	Align_NoAccel();
	}

void DSSAligner::CalcEvalue()
	{
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

	float Qual = 0;
	if (logE < -20)
		Qual = 1;
	else
		{
		float x = powf(10, logE/10);
		Qual = 1/(1 + x/2);
		}

	float E_scop = expf(logE)/SCOP40_DBSIZE;
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
	m_NewTestStatisticA = -FLT_MAX;
	m_NewTestStatisticB = -FLT_MAX;
	m_AlnFwdScore = 0;

	m_GlobalScore = -9999;
	m_GlobalPath.clear();
	}

void DSSAligner::Align_NoAccel()
	{
	ClearAlign();
	SetSMx_NoRev(*m_Params, *m_ProfileA, *m_ProfileB);

	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	StartTimer(SWFwd);
	uint Leni, Lenj;
	m_AlnFwdScore = SWFast(m_Mem, GetSMxData(), LA, LB,
	  m_Params->m_GapOpen, m_Params->m_GapExt,
	  m_LoA, m_LoB, Leni, Lenj, m_Path);
	EndTimer(SWFwd);

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

void DSSAligner::ToFasta2(FILE *f, bool Global, bool aUp) const
	{
	if (f == 0)
		return;

	const bool Up = !aUp;
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
	if (opt(noself) && m_ChainA->m_Label == m_ChainB->m_Label)
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

void DSSAligner::Stats()
	{
	uint Satn = m_ParasailSaturateCount;
	uint Disn = m_MuFilterDiscardCount;
	uint Inn = m_MuFilterInputCount;
	Log("DSSAligner::Stats() alns %s, mufil %u/%u %.1f%% (sat %u)",
	  FloatToStr(m_AlnCount), Inn, Disn, GetPct(Disn, Inn), Satn);
	Log(" xfil %.1f%%\n",
				GetPct(m_XDropDiscardCount, m_XDropAlnCount+1));
	MuKmerFilter::Stats();
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

float DSSAligner::GetKabsch(float t[3], float u[3][3], bool Up) const
	{
	if (Up)
		return Kabsch(*m_ChainA, *m_ChainB, m_LoA, m_LoB, m_Path, t, u);
	else
		{
		string Path;
		InvertPath(m_Path, Path);
		return Kabsch(*m_ChainB, *m_ChainA, m_LoB, m_LoA, Path, t, u);
		}
	}

void DSSAligner::AlignMKF()
	{
	ClearAlign();
	m_MKF.Align(*m_MuLettersB, *m_MuKmersB);
	PostAlignMKF();
	}

void DSSAligner::PostAlignMKF()
	{
	if (m_MKF.m_BestChainScore <= 0)
		return;

	++m_XDropAlnCount;
	float MegaHSPTotal = 0;
	const uint M = SIZE(m_MKF.m_ChainHSPLois);
	float BestMegaScore = 0;
	uint BestMegaIdx = 0;
	for (uint Idx = 0; Idx < M; ++Idx)
		{
		uint Loi = (uint) m_MKF.m_ChainHSPLois[Idx];
		uint Loj = (uint) m_MKF.m_ChainHSPLojs[Idx];
		uint Len = (uint) m_MKF.m_ChainHSPLens[Idx];
		float MegaScore = GetMegaHSPScore(Loi, Loj, Len);
		if (MegaScore > BestMegaScore)
			{
			BestMegaScore = MegaScore;
			BestMegaIdx = Idx;
			}
		MegaHSPTotal += MegaScore;
		}
	if (MegaHSPTotal < m_Params->m_MKF_MinMegaHSPScore)
		return;

	uint HSPLoA = (uint) m_MKF.m_ChainHSPLois[BestMegaIdx];
	uint HSPLoB = (uint) m_MKF.m_ChainHSPLojs[BestMegaIdx];
	uint Len = (uint) m_MKF.m_ChainHSPLens[BestMegaIdx];
	m_XDropScore = XDropHSP(HSPLoA, HSPLoB, Len,
							m_LoA, m_LoB, m_HiA, m_HiB);
	m_AlnFwdScore = m_XDropScore;
	m_Path = m_XDropPath;
	uint nM, nD, nI;
	GetPathCounts(m_Path, nM, nD, nI);
	m_HiA = m_LoA + nM + nD - 1;
	m_HiB = m_LoB + nM + nI - 1;
	CalcEvalue();
	}

const float * const *DSSAligner::GetSMxData() const
	{
	return m_SMx_Data;
	}

float **DSSAligner::GetSMxData()
	{
	return m_SMx_Data;
	}

void DSSAligner::AllocSMxData(uint LA, uint LB)
	{
	if (LA <= (2*m_SMx_Rows)/3 && LB <= (2*m_SMx_Cols)/3)
		return;
	FreeSMxData();
	
	size_t n = size_t(LA)*size_t(LB);
	uint un = uint(n);
	if (size_t(un) != n)
		Die("AllocSMxData(%u, %u) overflow", LA, LB);

	m_SMx_Buffer = (float *) malloc(un*sizeof(float));
	m_SMx_Data = (float **) malloc(LA*sizeof(float *));
	for (uint i = 0; i < LA; ++i)
		m_SMx_Data[i] = m_SMx_Buffer + i*LB;

	m_SMx_Rows = LA;
	m_SMx_Cols = LB;
	m_SMx_BufferSize = n;
	}

void DSSAligner::FreeSMxData()
	{
	if (m_SMx_Rows == 0)
		{
		asserta(m_SMx_BufferSize == 0);
		asserta(m_SMx_Buffer == 0);
		asserta(m_SMx_Data == 0);
		return;
		}

	free(m_SMx_Data);
	free(m_SMx_Buffer);

	m_SMx_Data = 0;
	m_SMx_Buffer = 0;
	m_SMx_BufferSize = 0;
	m_SMx_Rows = 0;
	m_SMx_Cols = 0;
	}
