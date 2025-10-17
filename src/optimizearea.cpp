#include "myutils.h"
#include "featuretrainer2.h"
#include "peaker.h"

static const vector<float> *s_AlnSubstScores;
static const vector<uint> *s_AlnColCountVec;
static const vector<uint> *s_AlnOpenVec;
static const vector<uint> *s_AlnExtVec;
static const vector<bool> *s_TPs;
uint s_AlnCount;
static float s_Area;
static Peaker s_Peaker;

static vector<float> s_AlnScores;

static double EvalArea(const vector<double> &xv)
	{
	asserta(SIZE(xv) == 2);
	float OpenPenalty = float(xv[0]);
	float ExtPenalty = float(xv[1]);
	float Bias = 0;

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
	float &OpenPenalty,
	float &ExtPenalty,
	float &Bias,
	float &Area)
	{
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
	SpecLines.push_back("var=open\tmin=-2\tmax=4\tdelta=0.2\tbins=32\tinit=1");
	SpecLines.push_back("var=ext\tmin=-1\tmax=1\tdelta=0.04\tbins=32\tinit=0.1");

	s_Peaker.Init(SpecLines, EvalArea);
	s_Peaker.Run();
	const vector<double> &xv = s_Peaker.GetBestParams();
	OpenPenalty = (float) xv[0];
	ExtPenalty = (float) xv[1];
	Bias = 0;
	Area = (float) s_Peaker.m_Best_y;
	}
