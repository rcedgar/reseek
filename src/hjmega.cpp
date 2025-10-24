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
	//AddFeature(FEATURE_AA,			0.398145f); // 0
	//AddFeature(FEATURE_NENDist,		0.129367f);	// 1
	//AddFeature(FEATURE_Conf,		0.202354f);		// 2
	//AddFeature(FEATURE_NENConf,		0.149383f);	// 3
	//AddFeature(FEATURE_RENDist,		0.0937677f);// 4
	//AddFeature(FEATURE_DstNxtHlx,	0.00475462f);	// 5
	//AddFeature(FEATURE_StrandDens,	0.0183853f);// 6
	//AddFeature(FEATURE_NormDens,	0.00384384f);	// 7

static double EvalArea(const vector<double> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);

	asserta(SIZE(DSSParams::m_Features) == 8);
	asserta(DSSParams::m_Features[0] == FEATURE_AA);
	asserta(DSSParams::m_Features[1] == FEATURE_NENDist);
	asserta(DSSParams::m_Features[2] == FEATURE_Conf);
	asserta(DSSParams::m_Features[3] == FEATURE_NENConf);
	asserta(DSSParams::m_Features[4] == FEATURE_RENDist);
	asserta(DSSParams::m_Features[5] == FEATURE_DstNxtHlx);
	asserta(DSSParams::m_Features[6] == FEATURE_StrandDens);
	asserta(DSSParams::m_Features[7] == FEATURE_NormDens);

	vector<float> Weights(8);
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		float Value = (float) xv[VarIdx];
		if (Value < 0)
			return DBL_MAX;
		const string VarName = string(s_Peaker->GetVarName(VarIdx));
		if (VarName == "open")				DSSParams::m_GapOpen = -Value;
		else if (VarName == "ext")			DSSParams::m_GapExt = -Value;
		else if (VarName == "AA")			Weights[0] = Value;
		else if (VarName == "NENDist")		Weights[1] = Value;
		else if (VarName == "Conf")			Weights[2] = Value;
		else if (VarName == "NENConf")		Weights[3] = Value;
		else if (VarName == "RENDist")		Weights[4] = Value;
		else if (VarName == "DstNxtHlx")	Weights[5] = Value;
		else if (VarName == "StrandDens")	Weights[6] = Value;
		else if (VarName == "NormDens")		Weights[7] = Value;
		else
			Die("VarName=%s", VarName.c_str());
		}

	DSSParams::UpdateWeights(Weights);
	s_SB->ClearHits();
	s_SB->RunSelf();
	s_SB->m_Level = "sf";
	s_SB->SetStats(0.005f);
	s_SB->WriteSummary();
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

void cmd_evalarea()
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
	vector<float> Weights(8);

	vector<string> SpecLines;
	vector<string> Fields;
	ReadLinesFromFile(SpecFN, SpecLines);
	for (uint i = 0; i < SIZE(SpecLines); ++i)
		{
		const string &Line = SpecLines[i].c_str();
		if (StartsWith(Line, "#"))
			continue;
		Log("%s\n", Line.c_str());
		Split(Line, Fields, '=');
		asserta(SIZE(Fields) == 2);
		const string &VarName = Fields[0];
		float Value = StrToFloatf(Fields[1]);
		if (VarName == "open")				DSSParams::m_GapOpen = -Value;
		else if (VarName == "ext")			DSSParams::m_GapExt = -Value;
		else if (VarName == "AA")			Weights[0] = Value;
		else if (VarName == "NENDist")		Weights[1] = Value;
		else if (VarName == "Conf")			Weights[2] = Value;
		else if (VarName == "NENConf")		Weights[3] = Value;
		else if (VarName == "RENDist")		Weights[4] = Value;
		else if (VarName == "DstNxtHlx")	Weights[5] = Value;
		else if (VarName == "StrandDens")	Weights[6] = Value;
		else if (VarName == "NormDens")		Weights[7] = Value;
		else
			Die("VarName=%s", VarName.c_str());
		}

	DSSParams::UpdateWeights(Weights);
	s_SB->ClearHits();
	s_SB->RunSelf();
	s_SB->m_Level = "sf";
	s_SB->SetStats(0.005f);
	s_SB->WriteSummary();
	ProgressLog("spec=%s;", SpecFN.c_str());
	ProgressLog("db=%s;", DBFN.c_str());
	ProgressLog("area=%.5g;", s_SB->m_Area);
	ProgressLog("\n");
	}
