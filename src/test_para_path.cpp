#include "myutils.h"
#include "dssaligner.h"
#include "trainer.h"
#include <set>

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

static vector<uint> g_Counts(11);
static uint g_Count;

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

	D.Init(ChainQ);
	D.GetComboLetters(ComboLettersQ);

	D.Init(ChainR);
	D.GetComboLetters(ComboLettersR);

	int GapOpen = D.m_Params->m_ParaComboGapOpen;
	int GapExt = D.m_Params->m_ParaComboGapExt;

	string Path;
	uint LoQ, LoR;
#if 0
	int Score = ParasailAlign(
	  ComboLettersQ, ComboLettersR, GapOpen, GapExt, LoQ, LoR, Path);
	if (Score < 0)
		return;
#else
	float Score = DA.AlignCombo(ComboLettersQ, ComboLettersR, LoQ, LoR, Path);
	vector<uint> PosQs;
	vector<uint> PosRs;
	GetPosVecs(LoQ, LoR, Path, PosQs, PosRs);
#endif
	++g_Count;
	double Cm = CmpVecs(PosQs, PosRs, TrPosQs, TrPosRs);
	uint Bin = uint(Cm*10);
	asserta(Bin >= 0 && Bin <= 10);
	g_Counts[Bin] += 1;
	}

void cmd_test_para_path()
	{
	DSSParams Params;
	Params.SetFromCmdLine(true);
	DA.m_Params = &Params;
	D.m_Params = &Params;

	Trainer Tr;
	Tr.Init(g_Arg1, opt_train_cal);
	Tr.Scan(OnPair, 0);

	double SumPct = 0;
	for (int Bin = 10; Bin >= 0; --Bin)
		{
		uint n = g_Counts[Bin];
		double Pct = GetPct(n, g_Count);
		SumPct += Pct;
		ProgressLog("%.4f  %5.1f%%  %7u\n", Bin/10.0, SumPct, n);
		}
	}
