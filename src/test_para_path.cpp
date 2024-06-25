#include "myutils.h"
#include "dssaligner.h"
#include "xdpmem.h"
#include "trainer.h"
#include <set>

#define	TEST_XDrop 1
#define	TEST_GetDPScorePath 0
#define	TEST_GetComboDPScorePath 0

int ParasailAlign(const vector<byte> &Q, const vector<byte> &T,
  int GapOpen, int GapExt, uint &LoQ, uint &LoT, string &Path);

static DSSAligner DA;
static DSS D;

static double CmpVecs(const vector<uint> &PosQs, const vector<uint> &PosRs,
  const vector<uint> &TrPosQs, const vector<uint> &TrPosRs)
	{
	const uint N = SIZE(PosQs);
	asserta(SIZE(PosRs) == N);
	set<pair<uint, uint> > PairSet;
	for (uint i = 0; i < N; ++i)
		PairSet.insert(pair<uint, uint>(PosQs[i], PosRs[i]));

	const uint TrN = SIZE(TrPosQs);
	asserta(SIZE(TrPosRs) == TrN);
	uint n = 0;
	for (uint i = 0; i < TrN; ++i)
		{
		pair<uint, uint> TrPair(TrPosQs[i], TrPosRs[i]);
		if (PairSet.find(TrPair) != PairSet.end())
			++n;
		}
	return double(n)/N;
	}

static void GetPosVecs(uint LoQ, uint LoR, const string &Path,
  vector<uint> &PosQs, vector<uint> &PosRs)
	{
	PosQs.clear();
	PosRs.clear();
	uint PosQ = LoQ;
	uint PosR = LoR;
	const uint ColCount = SIZE(Path);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			{
			PosQs.push_back(PosQ++);
			PosRs.push_back(PosR++);
			break;
			}
		case 'D':
			{
			++PosQ;
			break;
			}
		case 'I':
			{
			++PosR;
			break;
			}
		default:
			asserta(false);
			}
		}
	}

static void OnPair(const Trainer &Tr, uint ChainIdxQ, uint ChainIdxR,
  const vector<uint> &TrPosQs, const vector<uint> &TrPosRs)
	{
	asserta(SIZE(TrPosQs) == SIZE(TrPosRs));

	const uint ChainCount = SIZE(Tr.m_Chains);
	asserta(ChainIdxQ < ChainCount);
	asserta(ChainIdxR < ChainCount);
	const PDBChain &ChainQ = *Tr.m_Chains[ChainIdxQ];
	const PDBChain &ChainR = *Tr.m_Chains[ChainIdxR];

	vector<byte> ComboLettersQ;
	vector<byte> ComboLettersR;

	vector<vector<byte> > ProfileQ;
	vector<vector<byte> > ProfileR;

	D.Init(ChainQ);
	D.GetComboLetters(ComboLettersQ);
	D.GetProfile(ProfileQ);

	D.Init(ChainR);
	D.GetComboLetters(ComboLettersR);
	D.GetProfile(ProfileR);

	DA.Align_Test(ChainQ, ChainR, ComboLettersQ, ComboLettersR,
	  ProfileQ, ProfileR);
	}

static uint g_LA;
static uint g_LB;
static const vector<vector<byte> > *g_ProfileA;
static const vector<vector<byte> > *g_ProfileB;
static DSSAligner *g_DA;
static float SubFn(void *UserData, uint PosA, uint PosB)
	{
	float Score = g_DA->GetScorePosPair(*g_ProfileA, *g_ProfileB, PosA, PosB);
	return Score;
	}

void DSSAligner::Align_Test(
  const PDBChain &ChainA, const PDBChain &ChainB,
  const vector<byte> &ComboLettersA, const vector<byte> &ComboLettersB,
  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB)
	{
	SetQuery(ChainA, &ProfileA, 0, &ComboLettersA, FLT_MAX, FLT_MAX);
	SetTarget(ChainB, &ProfileB, 0, &ComboLettersB, FLT_MAX, FLT_MAX);

	m_EvalueA = FLT_MAX;
	m_EvalueB = FLT_MAX;
	m_PathA.clear();

	bool ComboFilterOk = ComboFilter();
	if (!ComboFilterOk)
		return;

#if TEST_GetComboDPScorePath
	SetSMx_Combo();
	uint LA = SIZE(ComboLettersA);
	uint LB = SIZE(ComboLettersB);
	float GapOpen = -(float) m_Params->m_ParaComboGapOpen;
	float GapExt = -(float) m_Params->m_ParaComboGapExt;
	uint LoA, LoB, LenA, LenB;
	string Path;
	float FwdScore = SWFast(m_Mem, m_SMx, LA, LB, GapOpen, GapExt,
	  LoA, LoB, LenA, LenB, Path);
	float FwdScore2 = GetComboDPScorePath(ComboLettersA, ComboLettersB,
	  LoA, LoB, GapOpen, GapExt, Path);
	asserta(feq(FwdScore, FwdScore2));
#endif

#if TEST_GetDPScorePath
	m_PathA.clear();
	m_LoA = UINT_MAX;
	m_LoB = UINT_MAX;

	SetSMx_NoRev();

	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	uint LenA, LenB;
	float FwdScore = SWFast(m_Mem, m_SMx, LA, LB,
	  m_Params->m_GapOpen, m_Params->m_GapExt,
	  m_LoA, m_LoB, LenA, LenB, m_PathA);
	float FwdScore2 = GetDPScorePath(*m_ProfileA, *m_ProfileB,
	  m_LoA, m_LoB, m_PathA);
	Log("%.1f %.1f\n", FwdScore, FwdScore2);
#endif

#if TEST_XDrop
	SetSMx_Combo();
	uint LA = SIZE(ComboLettersA);
	uint LB = SIZE(ComboLettersB);
	float GapOpen = -(float) m_Params->m_ParaComboGapOpen;
	float GapExt = -(float) m_Params->m_ParaComboGapExt;
	uint LoA, LoB, LenA, LenB;
	string Path;
	float FwdScore = SWFast(m_Mem, m_SMx, LA, LB, GapOpen, GapExt,
	  LoA, LoB, LenA, LenB, Path);
	const float X = 8;

	string FwdPathX;
	string BwdPathX;
	g_DA = this;
	g_ProfileA = m_ProfileA;
	g_ProfileB = m_ProfileB;
	float FwdScoreX = XDropFwd(m_Mem, X, GapOpen, GapExt, SubFn, 0, LoA, LA, LoB, LB, FwdPathX);
	float FwdScoreX2 = GetDPScorePath(*m_ProfileA, *m_ProfileB, LoA, LoB, FwdPathX);

	float BwdScoreX = XDropBwd(m_Mem, X, GapOpen, GapExt, SubFn, 0, LoA, LA, LoB, LB, BwdPathX);
	float BwdScoreX2 = GetDPScorePath(*m_ProfileA, *m_ProfileB, LoA, LoB, BwdPathX);
	Log("%.1f %.1f %.1f %.1f %.1f\n", FwdScore, FwdScoreX, FwdScoreX2, BwdScoreX, BwdScoreX2);
#endif
	}

void cmd_test_para_path()
	{
	DSSParams Params;
	Params.SetFromCmdLine(10000);
	DA.m_Params = &Params;
	D.m_Params = &Params;

	Trainer Tr;
	Tr.Init(g_Arg1, opt_train_cal);
	Tr.Scan(OnPair, 0);
	}
