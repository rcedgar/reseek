#include "myutils.h"
#include "seedtrainer.h"
#include "scop40bench.h"

bool Scop40_IsTP_SF(const string &Label1, const string &Label2)
	{
	string Dom1, Cls1, Fold1, SF1, Fam1;
	string Dom2, Cls2, Fold2, SF2, Fam2;
	SCOP40Bench::ParseScopLabel(Label1, Dom1, Cls1, Fold1, SF1, Fam1);
	SCOP40Bench::ParseScopLabel(Label2, Dom2, Cls2, Fold2, SF2, Fam2);
	bool t = (SF1 == SF2);
	return t;
	}

static uint s_NT;
static uint s_NF;

void OnPairT(Trainer &Tr, uint PairIdx)
	{
	++s_NT;
	}

void OnPairF(Trainer &Tr, uint ChainIdxQ, uint ChainIdxR)
	{
	++s_NF;
	}

void cmd_train_seeds()
	{
	opt_fast = true;
	optset_fast = true;

	DSSParams Params;
	Params.SetFromCmdLine(10000);

	SeedTrainer Tr;
	Tr.Init(Params, g_Arg1, opt_train_cal);
	Tr.EnumChainPairsT(OnPairT);
	Tr.EnumChainPairsF(Scop40_IsTP_SF, OnPairF);
	ProgressLog("NT %u, NF %u\n", s_NT, s_NF);
	}
