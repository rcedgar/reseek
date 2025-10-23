#include "myutils.h"
#include "statsig.h"
#include "scop40bench.h"
#include "peaker.h"

static SCOP40Bench *s_SB;

// SpecLines.push_back("var=AA\tmin=0\tmax=1\tdelta=0.05\tbins=8\tinit=0.2");			// 0
// SpecLines.push_back("var=NENDist\tmin=0\tmax=1\tdelta=0.05\tbins=8\tinit=0.2");		// 1
// SpecLines.push_back("var=Conf\tmin=0\tmax=1\tdelta=0.05\tbins=8\tinit=0.2");			// 2
// SpecLines.push_back("var=NENConf\tmin=0\tmax=1\tdelta=0.05\tbins=8\tinit=0.2");		// 3
// SpecLines.push_back("var=RENDist\tmin=0\tmax=1\tdelta=0.05\tbins=8\tinit=0.2");		// 4
// SpecLines.push_back("var=DstNxtHlx\tmin=0\tmax=1\tdelta=0.05\tbins=8\tinit=0.2");	// 5 
// SpecLines.push_back("var=StrandDens\tmin=0\tmax=1\tdelta=0.05\tbins=8\tinit=0.2");	// 6
// SpecLines.push_back("var=NormDens\tmin=0\tmax=1\tdelta=0.05\tbins=8\tinit=0.2");		// 7
// SpecLines.push_back("var=open\tmin=0\tmax=4\tdelta=0.2\tbins=8\tinit=1");			// 8
// SpecLines.push_back("var=ext\tmin=0\tmax=1\tdelta=0.04\tbins=8\tinit=0.1");			// 9
static double EvalArea(const vector<double> &xv)
	{
	asserta(SIZE(xv) == 10);

	float Open = float(xv[8]);
	float Ext = float(xv[9]);
	if (Open < 0)
		Open = 0;
	if (Ext < 0)
		Ext = 0;
	asserta(s_SB != 0);
	vector<float> Weights;
	for (uint i = 0; i < 8; ++i)
		Weights.push_back((float) xv[i]);
	DSSParams::UpdateWeights(Weights);
	DSSParams::m_GapOpen = -Open;
	DSSParams::m_GapExt = -Ext;
	s_SB->ClearHits();
	s_SB->RunSelf();
	s_SB->m_Level = "sf";
	s_SB->SetStats(0.005f);
	return s_SB->m_Area;
	}

void cmd_hjmega()
	{
	string CalFN = g_Arg1;

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
	s_SB->LoadDB(CalFN);
	StatSig::Init(s_SB->GetDBSize());
	s_SB->Setup();
	s_SB->m_QuerySelf = true;
	s_SB->m_ScoresAreEvalues = true;
	if (opt(scores_are_not_evalues))
		s_SB->m_ScoresAreEvalues = false;

	//vector<double> xv;
	//xv.push_back(0.685533);
	//xv.push_back(0.051881);
	//double Area = EvalArea(xv);
	//ProgressLog("Area=%.3g\n", s_SB->m_Area);

	vector<string> SpecLines;

	// y is Area, in range 0 to 1
	SpecLines.push_back("mindy=0.001");
	SpecLines.push_back("maxdy=0.1");
	SpecLines.push_back("minh=0.0005");
	SpecLines.push_back("latin=yes");
	SpecLines.push_back("sigfig=3");
	SpecLines.push_back("var=AA\tmin=0\tmax=1\tdelta=0.05\tbins=128\tinit=0.2");			// 0
	SpecLines.push_back("var=NENDist\tmin=0\tmax=1\tdelta=0.05\tbins=128\tinit=0.2");		// 1
	SpecLines.push_back("var=Conf\tmin=0\tmax=1\tdelta=0.05\tbins=128\tinit=0.2");		// 2
	SpecLines.push_back("var=NENConf\tmin=0\tmax=1\tdelta=0.05\tbins=128\tinit=0.2");		// 3
	SpecLines.push_back("var=RENDist\tmin=0\tmax=1\tdelta=0.05\tbins=128\tinit=0.2");		// 4
	SpecLines.push_back("var=DstNxtHlx\tmin=0\tmax=1\tdelta=0.05\tbins=128\tinit=0.2");	// 5 
	SpecLines.push_back("var=StrandDens\tmin=0\tmax=1\tdelta=0.05\tbins=128\tinit=0.2");	// 6
	SpecLines.push_back("var=NormDens\tmin=0\tmax=1\tdelta=0.05\tbins=128\tinit=0.2");	// 7
	SpecLines.push_back("var=open\tmin=0\tmax=4\tdelta=0.2\tbins=128\tinit=2");			// 8
	SpecLines.push_back("var=ext\tmin=0\tmax=1\tdelta=0.04\tbins=128\tinit=0.1");			// 9

	Peaker P;
	P.Init(SpecLines, EvalArea);
	P.Run();
	}
