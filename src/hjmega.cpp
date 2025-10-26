#include "myutils.h"
#include "statsig.h"
#include "scop40bench.h"
#include "peaker.h"

static SCOP40Bench *s_SB;
static Peaker *s_Peaker;

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
		Weights[0] = 1; // AA defaults to 1
		float Value = (float) xv[VarIdx];
		if (Value < 0)
			return DBL_MAX;
		const string VarName = string(s_Peaker->GetVarName(VarIdx));
		if (VarName == "open")				DSSParams::m_GapOpen = -Value;
		else if (VarName == "ext")			DSSParams::m_GapExt = -Value;
		else if (VarName == "gap")			{ DSSParams::m_GapOpen = DSSParams::m_GapExt = -Value; }
		else if (VarName == "gap2")			{ DSSParams::m_GapOpen = DSSParams::m_GapExt = -Value*0.1f; }
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
	s_SB->ClearHitsAndResults();
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

	DSSParams::Init(DM_UseCommandLineOption);
	string ParamStr;
	DSSParams::GetParamStr(ParamStr);
	ProgressLog("ParamStr:%s\n", ParamStr.c_str());

	s_SB = new SCOP40Bench;
	s_SB->LoadDB(DBFN);
	StatSig::Init(s_SB->GetDBSize());
	s_SB->Setup();
	s_SB->m_QuerySelf = true;

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
