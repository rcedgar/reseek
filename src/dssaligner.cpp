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
#include "statsig.h"
#include "alncounts.h"
#include <thread>
#include <set>
#include <mutex>

mutex DSSAligner::m_OutputLock;

void LogAln(const char *A, const char *B, const char *Path, unsigned ColCount);
float SWFast(XDPMem &Mem, const float * const *SMxData, uint LA, uint LB,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path);
void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);

atomic<uint> DSSAligner::m_AlnCount;
atomic<uint> DSSAligner::m_XDropAlnCount;
atomic<uint> DSSAligner::m_XDropDiscardCount;
atomic<uint> DSSAligner::m_SWCount;
atomic<uint> DSSAligner::m_MuFilterInputCount;
atomic<uint> DSSAligner::m_MuFilterDiscardCount;
atomic<uint> DSSAligner::m_ParasailSaturateCount;

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
	if (m_ProfPara8 != 0)	parasail_profile_free((parasail_profile_t *) m_ProfPara8);
	if (m_ProfPara16 != 0)	parasail_profile_free((parasail_profile_t *) m_ProfPara16);
	if (m_ProfParaRev8 != 0)	parasail_profile_free((parasail_profile_t *) m_ProfParaRev8);
	if (m_ProfParaRev16 != 0)	parasail_profile_free((parasail_profile_t *) m_ProfParaRev16);
	if (m_DProw != 0)
		myfree(m_DProw);
	FreeSMxData();
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
		m_UFs.push_back(UF_pvalue);
		}
	}

// GetDPScorePath calculates AlnScore which is optimized by SWFast.
float DSSAligner::GetDPScorePath(const vector<vector<byte> > &ProfileA,
  const vector<vector<byte> > &ProfileB, uint LoA, uint LoB,
  const string &Path) const
	{
	float Sum = 0;
	uint PosA = LoA;
	uint PosB = LoB;
	const float Open = DSSParams::m_GapOpen;
	const float Ext = DSSParams::m_GapExt;
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
			if (opt(tracedpscorepath))
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
				if (opt(tracedpscorepath))
					Log("De %5u  %5u  %10.3g  %10.3g\n",
					  PosA, PosB, Ext, Sum);
				}
			else
				{
				Sum += Open;
				if (opt(tracedpscorepath))
					Log("Do %5u  %5u  %10.3g  %10.3g\n",
					  PosA, PosB, Open, Sum);
				}
			++PosA;
			break;

		case 'I':
			if (Col != 0 && Path[Col-1] == 'I')
				{
				Sum += Ext;
				if (opt(tracedpscorepath))
					Log("Ie %5u  %5u  %10.3g  %10.3g\n",
					  PosA, PosB, Ext, Sum);
				}
			else
				{
				Sum += Open;
				if (opt(tracedpscorepath))
					Log("Io %5u  %5u  %10.3g  %10.3g\n",
					  PosA, PosB, Open, Sum);
				}
			++PosB;
			break;

		default:
			asserta(false);
			}
		}
	if (opt(tracedpscorepath))
		Log("Total score %.3g\n", Sum);
	return Sum;
	}

float DSSAligner::GetScorePosPair(const vector<vector<byte> > &ProfileA,
  const vector<vector<byte> > &ProfileB, uint PosA, uint PosB) const
	{
	const uint FeatureCount = DSSParams::GetFeatureCount();
	asserta(SIZE(ProfileA) == FeatureCount);
	asserta(SIZE(ProfileB) == FeatureCount);
	float MatchScore = 0;
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		asserta(PosA < SIZE(ProfileA[FeatureIdx]));
		asserta(PosB < SIZE(ProfileB[FeatureIdx]));
		FEATURE F = DSSParams::m_Features[FeatureIdx];
		uint AlphaSize = g_AlphaSizes2[F];
		float **ScoreMx = DSSParams::m_ScoreMxs[F];
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

void DSSAligner::SetSMx_QRev()
	{
	//const vector<vector<byte> > &ProfileA = *m_ProfileA;
	uint LA = m_ChainA->GetSeqLength();
	if (LA > 500)
		LA = 500;//FIXME

	//Mx<float> &SMx = m_SMx;
	//m_SMx.Alloc(LA, LA, __FILE__, __LINE__);
	AllocSMxData(LA, LA);
	StartTimer(SetSMx_QRev);
	float **Sim = GetSMxData();

	const uint FeatureCount = DSSParams::GetFeatureCount();
	asserta(SIZE((*m_ProfileA)) == FeatureCount);

// Special case first feature because = not += and
	FEATURE F0 = DSSParams::m_Features[0];
	uint AlphaSize0 = g_AlphaSizes2[F0];
	float **ScoreMx0 = DSSParams::m_ScoreMxs[F0];
	//const vector<byte> &ProfRowA = (*m_ProfileA)[0];
	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		byte ia = (*m_ProfileA)[0][PosA];
		float *SimRow = Sim[PosA];
		const float *ScoreMxRow = ScoreMx0[ia];

		for (uint PosB = 0; PosB < LA; ++PosB)
			{
			byte ib = (*m_ProfileA)[0][LA-1-PosB];
			assert(ia < AlphaSize0 && ib < AlphaSize0);
			float MatchScore = ScoreMxRow[ib];
			SimRow[PosB] = MatchScore;
			}
		}

	for (uint FeatureIdx = 1; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = DSSParams::m_Features[FeatureIdx];
		uint AlphaSize = g_AlphaSizes2[F];
		float **ScoreMx = DSSParams::m_ScoreMxs[F];
		//const vector<byte> &ProfRowA = (*m_ProfileA)[FeatureIdx];
		for (uint PosA = 0; PosA < LA; ++PosA)
			{
			byte ia = (*m_ProfileA)[FeatureIdx][PosA];
			float *SimRow = Sim[PosA];
			const float *ScoreMxRow = ScoreMx[ia];

			for (uint PosB = 0; PosB < LA; ++PosB)
				{
				byte ib = (*m_ProfileA)[FeatureIdx][LA-1-PosB];
				float MatchScore = 0;
				assert(ia < AlphaSize && ib < AlphaSize);
				MatchScore = ScoreMxRow[ib];
				SimRow[PosB] += MatchScore;
				}
			}
		}
#if DEBUG
	{
	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		for (uint PosB = 0; PosB < LA; ++PosB)
			{
			float MatchScore = Sim[PosA][PosB];
			float MatchScore2 = GetScorePosPair(*m_ProfileA, *m_ProfileA, PosA, LA-1-PosB);
			asserta(feq(MatchScore2, MatchScore));
			}
		}
	}
#endif
	EndTimer(SetSMx_QRev);
	}
//
//void DSSAligner::SetMuQP()
//	{
//	StartTimer(SetMuQP);
//	const vector<byte> &MuLettersA = *m_MuLettersA;
//	uint LA = SIZE(MuLettersA);
//	uint n = SIZE(m_ProfMu);
//	if (LA > n)
//		{
//		m_ProfMu.resize(LA);
//		m_ProfMuRev.resize(LA);
//		}
//	const uint AS = 36;
//	asserta(AS == 36);
//	for (uint PosA = 0; PosA < LA; ++PosA)
//		{
//		byte LetterA = MuLettersA[PosA];
//		asserta(LetterA < AS);
//		const float *MuMxRow = musubstmx[LetterA];
//		m_ProfMu[PosA] = MuMxRow;
//		m_ProfMuRev[LA-PosA-1] = MuMxRow;
//		}
//	EndTimer(SetMuQP);
//	}
//
//void DSSAligner::SetMuQPi()
//	{
//	const vector<byte> &MuLettersA = *m_MuLettersA;
//	uint LA = SIZE(MuLettersA);
//	uint n = SIZE(m_ProfMui);
//	if (LA > n)
//		{
//		m_ProfMui.resize(LA);
//		m_ProfMuRevi.resize(LA);
//		}
//	const uint AS = 36;
//	asserta(AS == 36);
//	for (uint PosA = 0; PosA < LA; ++PosA)
//		{
//		byte LetterA = MuLettersA[PosA];
//		asserta(LetterA < AS);
//		const int8_t *MuMxRow = IntScoreMx_Mu[LetterA];
//		m_ProfMui[PosA] = MuMxRow;
//		m_ProfMuRevi[LA-PosA-1] = MuMxRow;
//		}
//	}

void DSSAligner::LogHSP(uint Lo_i, uint Lo_j, uint Len) const
	{
	const char *SeqA = m_ChainA->m_Seq.c_str();
	const char *SeqB = m_ChainB->m_Seq.c_str();
	Log("\nDSSAligner::LogHSP(Loi=%u, Loj=%u, Len=%u)\n", Lo_i, Lo_j, Len);
	for (uint Col = 0; Col < Len; ++Col)
		Log("%c", SeqA[Lo_i+Col]);
	Log("\n");
	for (uint Col = 0; Col < Len; ++Col)
		Log("%c", SeqB[Lo_j+Col]);
	Log("\n");
	}

float DSSAligner::GetMegaHSPScore(uint Lo_i, uint Lo_j, uint Len)
	{
	StartTimer(GetMegaHSPScore);
	const vector<vector<byte> > &ProfileA = *m_ProfileA;
	const vector<vector<byte> > &ProfileB = *m_ProfileB;
	const uint FeatureCount = DSSParams::GetFeatureCount();
	asserta(SIZE(ProfileA) == FeatureCount);
	asserta(SIZE(ProfileB) == FeatureCount);

// Special case first feature because = not += and
	FEATURE F0 = DSSParams::m_Features[0];
	uint AlphaSize0 = g_AlphaSizes2[F0];
	float **ScoreMx0 = DSSParams::m_ScoreMxs[F0];
	const vector<byte> &ProfRowA = ProfileA[0];
	const vector<byte> &ProfRowB = ProfileB[0];
	float Total = 0;
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = DSSParams::m_Features[FeatureIdx];
		uint AlphaSize = g_AlphaSizes2[F];
		float **ScoreMx = DSSParams::m_ScoreMxs[F];
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

void DSSAligner::SetSMx_NoRev(const vector<vector<byte> > &ProfileA,
					  const vector<vector<byte> > &ProfileB)
	{
	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

// Memory blows up with grow-only strategy due to tail of long chains
	//m_SMx.Clear();
	//Mx<float> &SMx = m_SMx;
	//m_SMx.Alloc(LA, LB, __FILE__, __LINE__);
	AllocSMxData(LA, LB);
	StartTimer(SetSMx_NoRev);
	float **Sim = GetSMxData();

	const uint FeatureCount = DSSParams::GetFeatureCount();
	asserta(SIZE(ProfileA) == FeatureCount);
	asserta(SIZE(ProfileB) == FeatureCount);

// Special case first feature because = not += and
	FEATURE F0 = DSSParams::m_Features[0];
	uint AlphaSize0 = g_AlphaSizes2[F0];
	float **ScoreMx0 = DSSParams::m_ScoreMxs[F0];
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
		FEATURE F = DSSParams::m_Features[FeatureIdx];
		uint AlphaSize = g_AlphaSizes2[F];
		float **ScoreMx = DSSParams::m_ScoreMxs[F];
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

void DSSAligner::SetMuScore()
	{
	AlignMuQP_xx(*m_MuLettersA, *m_MuLettersB);
	}

bool DSSAligner::MuFilter()
	{
	if (m_MuLettersA == 0 || m_MuLettersB == 0)
		return true;

	if (DSSParams::m_ParaBits == 2)
		{
		int Score8 = AlignMuQP_Para8();
		if (Score8 < DSSParams::m_Omega8)
			return false;
		int Score16 = AlignMuQP_Para16();
		if (Score16 < DSSParams::m_Omega16)
			return false;
		return true;
		}

	const int MinMuScore = DSSParams::GetOmega();
	if (MinMuScore <= 0)
		return true;
	SetMuScore();
	if (m_MuFwdMinusRevScore < MinMuScore)
		return false;
	return true;
	}

void DSSAligner::Align_MuFilter(
  const PDBChain &ChainA, const PDBChain &ChainB,
  const vector<byte> &MuLettersA, const vector<uint> &MuKmersA,
  const vector<byte> &MuLettersB,const vector<uint> &MuKmersB,
  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB)
	{
	SetQuery(ChainA, &ProfileA, &MuLettersA, &MuKmersA, FLT_MAX);
	SetTarget(ChainB, &ProfileB, &MuLettersB, &MuKmersB, FLT_MAX);

	//m_EvalueA = FLT_MAX;
	//m_EvalueB = FLT_MAX;
	//m_Path.clear();
	ClearAlign();

	++m_AlnCount;

	bool MuFilterOk = MuFilter();
	++m_MuFilterInputCount;
	if (!MuFilterOk)
		{
		++m_MuFilterDiscardCount;
		return;
		}
	Align_NoAccel();
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
	if (ptrMuLetters != 0 && DSSParams::GetNeedMuLetters())
		SetMuQP_Para_xx();
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
	if (m_MuLettersA->empty() || m_MuLettersB->empty())
		return false;
	if (m_MuKmersA->empty() || m_MuKmersB->empty())
		return false;
	uint LA = m_ChainA->GetSeqLength();
	uint LB = m_ChainB->GetSeqLength();
	if (LA >= DSSParams::m_MKFL)
		return true;
	if (LB >= DSSParams::m_MKFL)
		return true;
	return false;
	}

void DSSAligner::AlignQueryTarget_Trace()
	{
	Log("\n");
	Log("______________________________________\n");
	Log("A>%s(%u)\n", m_ChainA->m_Label.c_str(), m_ChainA->GetSeqLength());
	Log("B>%s(%u)\n", m_ChainB->m_Label.c_str(), m_ChainB->GetSeqLength());

	void ProfileToLines(const vector<vector<byte> > &Profile,
		vector<string> &Lines);
	Log("\n");
	Log("ProfileA ==>\n");
	vector<string> Lines;
	ProfileToLines(*m_ProfileA, Lines);
	for (uint i = 0; i < SIZE(Lines); ++i)
		Log("%10.10s  %s\n", FeatureToStr(DSSParams::m_Features[i]), Lines[i].c_str());

	Log("\n");
	Log("ProfileB ==>\n");
	ProfileToLines(*m_ProfileB, Lines);
	for (uint i = 0; i < SIZE(Lines); ++i)
		Log("%10.10s  %s\n", FeatureToStr(DSSParams::m_Features[i]), Lines[i].c_str());

	ClearAlign();

	if (DoMKF())
		{
		Log("DoMKF()=true\n");
		AlignMKF();
		Log("m_LoA=%u, m_LoB=%u, m_HiA=%u, m_HiB=%u\n", m_LoA, m_LoB, m_HiA, m_HiB);
		Log("m_MKF.m_BestChainScore=%d\n", m_MKF.m_BestChainScore);
		Log("m_XDropScore=%.6g, Fwd=%.6g, Bwd=%.6g\n",
			m_XDropScore, m_XDropScoreFwd, m_XDropScoreBwd);
		Log("AlnFwdScore=%.6g\n", m_AlnFwdScore);
		Log("EvalueA=%.6g\n", m_EvalueA);
		Log("Path=(%u)%.10s...\n", SIZE(m_Path), m_Path.c_str());
		return;
		}

	int Omega = DSSParams::GetOmega();
	if (Omega > 0)
		{
		Log("Omega > 0\n");
		++m_MuFilterInputCount;
		bool MuFilterOk = MuFilter();
		Log("MuFilterOk=%c\n", tof(MuFilterOk));
		if (!MuFilterOk)
			{
			++m_MuFilterDiscardCount;
			return;
			}
		}

	SetSMx_NoRev(*m_ProfileA, *m_ProfileB);
	if (m_DoTrace)
		{
		ProgressLog("Trace SWMx\n");
		for (uint i = 0; i < 10; ++i)
			{
			for (uint j = 0; j < 10; ++j)
				Log("  %8.3g", m_SMx_Data[i][j]);
			Log("\n");
			}
		}

	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	uint Leni, Lenj;
	m_AlnFwdScore = SWFast(m_Mem, GetSMxData(), LA, LB,
	  DSSParams::m_GapOpen, DSSParams::m_GapExt,
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

void DSSAligner::AlignQueryTarget()
	{
	incac(alnqts);
	if (optset_label1 && optset_label2)
		{
		const string Label1 = string(opt(label1));
		const string Label2 = string(opt(label2));
		m_DoTrace = (m_ChainA->m_Label == Label1 && m_ChainB->m_Label == Label2) ||
			(m_ChainA->m_Label == Label2 && m_ChainB->m_Label == Label1);
		if (m_DoTrace)
			{
			ProgressLog("Trace %s, %s\n", Label1.c_str(), Label2.c_str());
			AlignQueryTarget_Trace();
			return;
			}
		}

	ClearAlign();

	if (DoMKF())
		{
		incac(alnmkfs);
		AlignMKF();
		return;
		}

	++m_AlnCount;

	if (DSSParams::GetDoMuFilter())
		{
		incac(mufilters);
		++m_MuFilterInputCount;
		bool MuFilterOk = MuFilter();
		if (!MuFilterOk)
			{
			incac(mufilterdiscards);
			++m_MuFilterDiscardCount;
			return;
			}
		}

	Align_NoAccel();
	}

void DSSAligner::CalcEvalue()
	{
// Threshold enables small speedup by avoiding LDDT and self-rev
	if (m_AlnFwdScore < DSSParams::m_MinFwdScore)
		return;

	uint M, D, I;
	GetPathCounts(m_Path, M, D, I);
	m_HiA = m_LoA + M + D - 1;
	m_HiB = m_LoB + M + I - 1;
	m_Ids = M;
	m_Gaps = D + I;

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

	m_NewTestStatisticA = DSSParams::m_lddtw*LDDT;
	m_NewTestStatisticA += (DSSParams::m_dpw*m_AlnFwdScore - DSSParams::m_revtsw*RevDPScore)/(L + DSSParams::m_ladd);

	m_NewTestStatisticB = m_NewTestStatisticA;

	float Pval = (float) StatSig::GetPvalue(m_NewTestStatisticA);
	float Qual = (float) StatSig::GetQual(m_NewTestStatisticA);
	float E = FLT_MAX;
	if (StatSig::m_DBSize != UINT_MAX)
		E = (float) StatSig::GetEvalue(m_NewTestStatisticA);

	m_QualityA = Qual;
	m_QualityB = Qual;
	m_PvalueA = Pval;
	m_PvalueB = Pval;
	m_EvalueA = E;
	m_EvalueB = E;
	EndTimer(CalcEvalue)
	}

void DSSAligner::ClearAlign()
	{
	m_Path.clear();
	m_XDropPath.clear();
	m_LoA = UINT_MAX;
	m_LoB = UINT_MAX;
	m_HiA = UINT_MAX;
	m_HiB = UINT_MAX;
	m_Ids = UINT_MAX;
	m_Gaps = UINT_MAX;
	m_PvalueA = FLT_MAX;
	m_PvalueB = FLT_MAX;
	m_EvalueA = FLT_MAX;
	m_EvalueB = FLT_MAX;
	m_QualityA = FLT_MAX;
	m_QualityB = FLT_MAX;
	m_TestStatisticA = -FLT_MAX;
	m_TestStatisticB = -FLT_MAX;
	m_NewTestStatisticA = -FLT_MAX;
	m_NewTestStatisticB = -FLT_MAX;
	m_XDropScore = 0;
	m_AlnFwdScore = 0;
	m_MuFwdScore8 = 0;
	m_MuRevScore8 = 0;
	m_MuFwdScore16 = 0;
	m_MuRevScore16 = 0;
	m_MuFwdMinusRevScore = 0;
	}

void DSSAligner::ClearAlign_ExceptMu()
	{
	m_Path.clear();
	m_XDropPath.clear();
	m_LoA = UINT_MAX;
	m_LoB = UINT_MAX;
	m_HiA = UINT_MAX;
	m_HiB = UINT_MAX;
	m_Ids = UINT_MAX;
	m_Gaps = UINT_MAX;
	m_PvalueA = FLT_MAX;
	m_PvalueB = FLT_MAX;
	m_EvalueA = FLT_MAX;
	m_EvalueB = FLT_MAX;
	m_QualityA = FLT_MAX;
	m_QualityB = FLT_MAX;
	m_TestStatisticA = -FLT_MAX;
	m_TestStatisticB = -FLT_MAX;
	m_NewTestStatisticA = -FLT_MAX;
	m_NewTestStatisticB = -FLT_MAX;
	m_XDropScore = 0;
	m_AlnFwdScore = 0;
	}

void DSSAligner::Align_NoAccel()
	{
	incac(alnnoaccels);
	ClearAlign_ExceptMu();
	SetSMx_NoRev(*m_ProfileA, *m_ProfileB);

	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	StartTimer(SWFwd);
	uint Leni, Lenj;
	if (m_DoTrace)
		{
		ProgressLog("Trace SWMx");
		for (uint i = 0; i < 10; ++i)
			{
			for (uint j = 0; j < 10; ++j)
				Log("  %8.3g", m_SMx_Data[i][j]);
			Log("\n");
			}
		}
	m_AlnFwdScore = SWFast(m_Mem, GetSMxData(), LA, LB,
	  DSSParams::m_GapOpen, DSSParams::m_GapExt,
	  m_LoA, m_LoB, Leni, Lenj, m_Path);
	EndTimer(SWFwd);

	CalcEvalue();
	}

void DSSAligner::Align_QRev()
	{
	ClearAlign();
	SetSMx_QRev();
	m_AlnFwdScore = 0;
	return;//FIXME
	uint LA = m_ChainA->GetSeqLength();
	if (LA > 500)
		LA = 500;//FIXME

	StartTimer(SWFwd);
	uint Leni, Lenj;
	m_AlnFwdScore = SWFast(m_Mem, GetSMxData(), LA, LA,
	  DSSParams::m_GapOpen, DSSParams::m_GapExt,
	  m_LoA, m_LoB, Leni, Lenj, m_Path);
	EndTimer(SWFwd);
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

int DSSAligner::AlignMuQP_Para_xx()
	{
	if (DSSParams::m_ParaBits == 8)
		return AlignMuQP_Para8();
	if (DSSParams::m_ParaBits == 16)
		return AlignMuQP_Para16();
	Die("AlignMuQP_Para_xx");
	return 0;
	}

int DSSAligner::AlignMuParaBags_xx(const ChainBag &BagA, const ChainBag &BagB)
	{
	if (DSSParams::m_ParaBits == 8)
		return AlignMuParaBags8(BagA, BagB);
	if (DSSParams::m_ParaBits == 16)
		return AlignMuParaBags16(BagA, BagB);
	Die("AlignMuParaBags_xx");
	return 0;
	}

void DSSAligner::SetMuQP_Para_xx()
	{
	if (DSSParams::m_ParaBits == 8)
		SetMuQP_Para8();
	else if (DSSParams::m_ParaBits == 16)
		SetMuQP_Para16();
	else if (DSSParams::m_ParaBits == 2)
		{
		SetMuQP_Para8();
		SetMuQP_Para16();
		}
	else
		Die("SetMuQP_Para_xx");
	}

int DSSAligner::AlignMuQP_xx(const vector<byte> &LettersA,
  const vector<byte> &LettersB)
	{
	m_MuLettersA = &LettersA;
	m_MuLettersB = &LettersB;
	return AlignMuQP_Para_xx();
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
	double GetLDDT_mu_fast(const PDBChain &Q, const PDBChain &T,
	  const vector<uint> &PosQs, const vector<uint> &PosTs);
	vector<uint> PosAs;
	vector<uint> PosBs;
	GetPosABs(PosAs, PosBs);
	double LDDT = 
	  GetLDDT_mu_fast(*m_ChainA, *m_ChainB, PosAs, PosBs);
	return (float) LDDT;
	}

float DSSAligner::GetBitScore() const
	{
	// S = Smith-Waterman score
	// S' = (lambda S - ln K)/ ln 2
	TODO
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
	m_MKF.m_DA = this;
	m_MKF.Align(*m_MuLettersB, *m_MuKmersB);
	PostAlignMKF();
	}

void DSSAligner::PostAlignMKF()
	{
	incac(postalignmkfs);
	if (m_MKF.m_BestChainScore <= 0)
		return;

	incac(postaligntryxdrops);
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
	if (MegaHSPTotal < DSSParams::m_MKF_MinMegaHSPScore)
		{
		incac(postalignlohsps);
		return;
		}

	incac(postalignxdrops);
	uint HSPLoA = (uint) m_MKF.m_ChainHSPLois[BestMegaIdx];
	uint HSPLoB = (uint) m_MKF.m_ChainHSPLojs[BestMegaIdx];
	uint Len = (uint) m_MKF.m_ChainHSPLens[BestMegaIdx];
#if TRACE_XDROP
	LogHSP(HSPLoA, HSPLoB, Len);
#endif
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

float DSSAligner::GetSBScore(SBSCORE SBS, bool Up) const
	{
	switch (SBS)
		{
	case SBS_Evalue:	return GetEvalue(Up);
	case SBS_TS:		return GetNewTestStatistic(Up);
	case SBS_Raw:		return m_AlnFwdScore;
		}
	Die("DSSAligner::GetSBScore(%d)", int(SBS));
	return FLT_MAX;
	}
