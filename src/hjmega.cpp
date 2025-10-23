#include "myutils.h"
#include "statsig.h"
#include "scop40bench.h"
#include "peaker.h"

static SCOP40Bench *s_SB;
static Peaker *s_Peaker;

/***
mindy=0.001
maxdy=0.1
minh=0.0001
latin=yes
sigfig=3
var=AA	min=0	max=1	delta=0.05	bins=64	init=0.40
var=NENDist	min=0	max=1	delta=0.05	bins=64	init=0.13
var=Conf	min=0	max=1	delta=0.05	bins=64	init=0.20
var=NENConf	min=0	max=1	delta=0.05	bins=64	init=0.15
var=RENDist	min=0	max=1	delta=0.05	bins=64	init=0.10
var=DstNxtHlx	min=0	max=1	delta=0.05	bins=64	init=0
var=StrandDens	min=0	max=1	delta=0.05	bins=64	init=0.02
var=NormDens	min=0	max=1	delta=0.05	bins=64	init=0
var=open	min=0	max=4	delta=0.2	bins=64	init=0.7
var=ext	min=0	max=1	delta=0.04	bins=64	init=0.05
***/
static double EvalArea(const vector<double> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);

	vector<float> Weights;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		float Value = (float) xv[VarIdx];
		const string VarName = string(s_Peaker->GetVarName(VarIdx));
		if (VarName == "open")
			{
			if (Value < 0)
				Value = 0;
			DSSParams::m_GapOpen = -Value;
			}
		else if (VarName == "ext")
			{
			if (Value < 0)
				Value = 0;
			DSSParams::m_GapExt = -Value;
			}
			
#define F(x)	else if (VarName == string(#x)) { Weights.push_back(Value); }
#include "featurelist.h"

		else
			Die("VarName=%s", VarName.c_str());
		}

	DSSParams::UpdateWeights(Weights);
	s_SB->ClearHits();
	s_SB->RunSelf();
	s_SB->m_Level = "sf";
	s_SB->SetStats(0.005f);
	return s_SB->m_Area;
	}

void cmd_hjmega()
	{
	const string SpecFN = g_Arg1;
	asserta(optset_db);
	const string &DBFN = opt(db);

	if (!optset_evalue)
		{
		opt_evalue = 9999;
		optset_evalue = true;
		}

	if (!optset_mints)
		{
		opt_mints = -99999;
		optset_mints = true;
		}

	DSSParams::Init(DM_UseCommandLineOption);

	s_SB = new SCOP40Bench;
	s_SB->LoadDB(DBFN);
	StatSig::Init(s_SB->GetDBSize());
	s_SB->Setup();
	s_SB->m_QuerySelf = true;
	s_SB->m_ScoresAreEvalues = true;
	if (opt(scores_are_not_evalues))
		s_SB->m_ScoresAreEvalues = false;

	Log("SpecFN=%s\n", SpecFN.c_str());
	vector<string> SpecLines;
	ReadLinesFromFile(SpecFN, SpecLines);
	for (uint i = 0; i < SIZE(SpecLines); ++i)
		Log("%s\n", SpecLines[i].c_str());

	Peaker P;
	s_Peaker = &P;
	P.Init(SpecLines, EvalArea);
	P.Run();
	}
