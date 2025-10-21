#include "myutils.h"
#include "featuretrainer2.h"
#include "peaker.h"
#include "sort.h"

static const vector<float> *s_AlnSubstScores;
static const vector<uint> *s_AlnColCountVec;
static const vector<uint> *s_AlnOpenVec;
static const vector<uint> *s_AlnExtVec;
static const vector<bool> *s_TPs;
uint s_AlnCount;
static float s_Area;
static Peaker s_Peaker;

static vector<float> s_AlnScores;

static double EvalArea_OpenOnly(const vector<double> &xv)
	{
	asserta(SIZE(xv) == 1);
	float OpenPenalty = float(xv[0]);
	float ExtPenalty = OpenPenalty*0.1f;
	float Bias = 0; // float(xv[2]);

	FeatureTrainer2::GetAlnScores(*s_AlnSubstScores, *s_AlnColCountVec,
		*s_AlnOpenVec, *s_AlnExtVec, OpenPenalty, ExtPenalty, Bias,
		s_AlnScores);

	float Area = FeatureTrainer2::CalcArea(s_AlnScores, *s_TPs);
	return Area;
	}

static double EvalArea_OPenExt(const vector<double> &xv)
	{
	asserta(SIZE(xv) == 2);
	float OpenPenalty = float(xv[0]);
	float ExtPenalty = float(xv[1]);
	float Bias = 0; // float(xv[2]);

	FeatureTrainer2::GetAlnScores(*s_AlnSubstScores, *s_AlnColCountVec,
		*s_AlnOpenVec, *s_AlnExtVec, OpenPenalty, ExtPenalty, Bias,
		s_AlnScores);

	float Area = FeatureTrainer2::CalcArea(s_AlnScores, *s_TPs);
	return Area;
	}

void FeatureTrainer2::OptimizeArea(
	const vector<float> &AlnSubstScores,
	const vector<uint> &AlnColCountVec,
	const vector<uint> &AlnOpenVec,
	const vector<uint> &AlnExtVec,
	const vector<bool> &TPs,
	float &BestOpenPenalty,
	float &BestExtPenalty,
	float &BestBias,
	float &BestArea,
	uint Iters)
	{
	BestOpenPenalty = FLT_MAX;
	BestExtPenalty = FLT_MAX;
	BestBias = FLT_MAX;
	BestArea = FLT_MAX;

	s_AlnCount = SIZE(AlnSubstScores);
	asserta(SIZE(AlnColCountVec) == s_AlnCount);
	asserta(SIZE(AlnOpenVec) == s_AlnCount);
	asserta(SIZE(AlnExtVec) == s_AlnCount);
	asserta(SIZE(TPs) == s_AlnCount);

	s_AlnSubstScores = &AlnSubstScores;
	s_AlnColCountVec= &AlnColCountVec;
	s_AlnOpenVec = &AlnOpenVec;
	s_AlnExtVec = &AlnExtVec;
	s_TPs = &TPs;

	vector<string> SpecLines;

	// y is Area, in range 0 to 1
	SpecLines.push_back("mindy=0.001");
	SpecLines.push_back("maxdy=0.1");
	SpecLines.push_back("minh=0.0005");
	SpecLines.push_back("latin=yes");
	SpecLines.push_back("sigfig=3");
	SpecLines.push_back("var=open\tmin=0\tmax=4\tdelta=0.1\tbins=32\tinit=0.5");
	//SpecLines.push_back("var=ext\tmin=0\tmax=1\tdelta=0.01\tbins=32\tinit=0.1");
	//SpecLines.push_back("var=bias\tmin=-1\tmax=1\tdelta=0.01\tbins=32\tinit=0.0");

	vector<float> Areas;
	vector<float> Opens;
	vector<float> Exts;
	vector<float> Biases;

	vector<double> xv0 = { 0 };
	float Area0 = (float) EvalArea_OpenOnly(xv0);
	Psa(m_FevStr, "area0=%.3g;", Area0);
	Areas.push_back(Area0);
	Opens.push_back(0);
	Exts.push_back(0);
	Biases.push_back(0);

	BestArea = 0;
	for (uint Iter = 0; Iter < Iters; ++Iter)
		{
		s_Peaker.m_Progress = false;
		s_Peaker.Init(SpecLines, EvalArea_OpenOnly);
		s_Peaker.Run();
		const vector<double> &xv = s_Peaker.GetBestParams();
		asserta(SIZE(xv) == 1);
		float OpenPenalty = (float) xv[0];
		float ExtPenalty = OpenPenalty*0.1f; // (float) xv[1];
		float Bias = 0; // (float) xv[2];
		float Area = (float) s_Peaker.m_Best_y;

		Areas.push_back(Area);
		Opens.push_back(OpenPenalty);
		Exts.push_back(ExtPenalty);
		Biases.push_back(Bias);

		if (Area > BestArea)
			{
			BestArea = Area;
			BestOpenPenalty = OpenPenalty;
			BestExtPenalty = ExtPenalty;
			BestBias = Bias;
			}
		}
	vector<uint> Order(Iters);
	QuickSortOrder(Areas.data(), Iters, Order.data());
	ProgressLog("%s\n", m_FevStr.c_str());
	ProgressLog("  Area    Diff    Open     Ext    Bias\n");
	//   123456  123456  123456  123456  123456
	for (uint k = 0; k < Iters; ++k)
		{
		uint Iter = Order[k];
		float Area = Areas[Iter];
		float Open = Opens[Iter];
		float Ext = Exts[Iter];
		float Bias = Biases[Iter];
		float Diff = BestArea - Area;
		ProgressLog("%6.3f", Area);
		if (Diff < 0.001)
			ProgressLog("  %6.6s", "");
		else
			ProgressLog("  %6.3f", Diff);
		ProgressLog("  %6.3f", Open);
		ProgressLog("  %6.3f", Ext);
		ProgressLog("  %6.3f", Bias);
		ProgressLog("\n");
		}
	}
