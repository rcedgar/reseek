#include "myutils.h"
#include "dssaligner.h"
#include "xdpmem.h"
#include "trainer.h"
#include <set>

#define	TEST_XDrop 1
#define	TEST_GetDPScorePath 0
#define	TEST_GetMuDPScorePath 0

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

	vector<byte> MuLettersQ;
	vector<byte> MuLettersR;

	vector<vector<byte> > ProfileQ;
	vector<vector<byte> > ProfileR;

	D.Init(ChainQ);
	D.GetMuLetters(MuLettersQ);
	D.GetProfile(ProfileQ);

	D.Init(ChainR);
	D.GetMuLetters(MuLettersR);
	D.GetProfile(ProfileR);

	DA.Align_Test(ChainQ, ChainR, MuLettersQ, MuLettersR,
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
  const vector<byte> &MuLettersA, const vector<byte> &MuLettersB,
  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB)
	{
	SetQuery(ChainA, &ProfileA, &MuLettersA, 0, FLT_MAX);
	SetTarget(ChainB, &ProfileB, &MuLettersB, 0, FLT_MAX);

	m_EvalueA = FLT_MAX;
	m_EvalueB = FLT_MAX;
	m_Path.clear();

	bool MuFilterOk = MuFilter();
	if (!MuFilterOk)
		return;

#if TEST_GetMuDPScorePath
	SetSMx_Mu();
	uint LA = SIZE(MuLettersA);
	uint LB = SIZE(MuLettersB);
	float GapOpen = -(float) m_Params->m_ParaMuGapOpen;
	float GapExt = -(float) m_Params->m_ParaMuGapExt;
	uint LoA, LoB, LenA, LenB;
	string Path;
	float FwdScore = SWFast(m_Mem, m_SMx, LA, LB, GapOpen, GapExt,
	  LoA, LoB, LenA, LenB, Path);
	float FwdScore2 = GetMuDPScorePath(MuLettersA, MuLettersB,
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

#if 0 // TEST_XDrop
	//SetSMx_Mu(); 
	Die("TODO");
	uint LA = SIZE(MuLettersA);
	uint LB = SIZE(MuLettersB);
	float GapOpen = -(float) m_Params->m_ParaMuGapOpen;
	float GapExt = -(float) m_Params->m_ParaMuGapExt;
	uint LoA, LoB, LenA, LenB;
	string Path;
	//float FwdScore = SWFast(m_Mem, m_SMx.GetData(), LA, LB, GapOpen, GapExt,
	//  LoA, LoB, LenA, LenB, Path);
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
	DA.SetParams(Params);
	D.SetParams(Params);

	Trainer Tr;
	Tr.Init(g_Arg1, opt_train_cal);
	Tr.Scan(OnPair, 0);
	}
