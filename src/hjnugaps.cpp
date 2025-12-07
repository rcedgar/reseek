#include "myutils.h"
#include "featuretrainer2.h"
#include "parasearch.h"
#include <list>

static ParaSearch *s_PS;
static double s_BestSum3 = -999;
static int s_BestOpen = -999;
static int s_BestExt = -999;

static list<int> s_PendingOpens;
static list<int> s_PendingExts;

static vector<int> s_DoneOpens;
static vector<int> s_DoneExts;
static vector<double> s_DoneSum3s;

static vector<FEATURE> s_Fs;
static vector<float> s_Weights;

static double EvalSum3(int IntOpen, int IntExt)
	{
	s_PS->SetGapParams(IntOpen, IntExt);
	s_PS->ClearHitsAndResults();
	s_PS->Search("para", false);
	s_PS->Bench();
	return s_PS->m_Sum3;
	}

static void GetFeatures(const string &s,
	vector<FEATURE> &Fs, vector<float> &Weights)
	{
	Fs.clear();
	Weights.clear();
	vector<string> Fields;
	Split(s, Fields, ';');
	const uint n = SIZE(Fields);
	for (uint Idx = 0; Idx < n; ++Idx)
		{
		vector<string> Fields2;
		Split(Fields[Idx], Fields2, '=');
		asserta(SIZE(Fields2) == 2);
		FEATURE F = StrToFeature(Fields2[0].c_str());
		float Weight = StrToFloatf(Fields2[1]);
		Fs.push_back(F);
		Weights.push_back(Weight);
		}
	}

static bool ExtStalled(int Open, int Ext)
	{
	double Score3_Minus1 = DBL_MAX;
	double Score3_Minus2 = DBL_MAX;
	double Score3_Minus3 = DBL_MAX;
	const uint N = SIZE(s_DoneOpens);
	asserta(SIZE(s_DoneExts) == N);
	for (uint i = 0; i < N; ++i)
		{
		int Open2 = s_DoneOpens[i];
		int Ext2 = s_DoneExts[i];
		if (Open2 == Open && Ext2 == Ext - 1)
			Score3_Minus1 = s_DoneSum3s[i];
		if (Open2 == Open && Ext2 == Ext - 2)
			Score3_Minus2 = s_DoneSum3s[i];
		if (Open2 == Open && Ext2 == Ext - 3)
			Score3_Minus3 = s_DoneSum3s[i];
		}
	if (Score3_Minus1 == DBL_MAX ||
		Score3_Minus2 == DBL_MAX ||
		Score3_Minus3 == DBL_MAX)
		return false;
	return	feq(Score3_Minus1, Score3_Minus2) &&
			feq(Score3_Minus2, Score3_Minus3);
	}

static bool AddPendingIfOk(int Open, int Ext)
	{
	if (Open < 0 || Ext < 0)
		return false;

	bool OpenFound = (std::find(s_DoneOpens.begin(), s_DoneOpens.end(), Open) != s_DoneOpens.end());
	bool ExtFound = (std::find(s_DoneExts.begin(), s_DoneExts.end(), Ext) != s_DoneExts.end());
	if (OpenFound && ExtFound)
		return false;

	if (ExtStalled(Open, Ext))
		{
		ProgressLog("\nSTALLED %d/%d\n\n", Open, Ext);
		return false;
		}

	s_PendingOpens.push_back(Open);
	s_PendingExts.push_back(Ext);

	return true;
	}

static void Optimize(int ScaleFactor, int FirstOpen, int FirstExt)
	{
	s_BestSum3 = -999;
	s_BestOpen = -999;
	s_BestExt = -999;

	s_PendingOpens.clear();
	s_PendingExts.clear();

	s_DoneOpens.clear();
	s_DoneExts.clear();
	s_DoneSum3s.clear();

	Paralign::SetCompoundMx(s_Fs, s_Weights, ScaleFactor, FirstOpen, FirstExt, 777);

	AddPendingIfOk(FirstOpen, FirstExt);
	for (uint Iter = 0; Iter < 100; ++Iter)
		{
		uint PendingCount = SIZE(s_PendingOpens);
		asserta(SIZE(s_PendingExts) == PendingCount);
		if (PendingCount == 0)
			break;
		int Open = s_PendingOpens.front();
		int Ext = s_PendingExts.front();
		s_PendingOpens.pop_front();
		s_PendingExts.pop_front();
		double Sum3 = EvalSum3(Open, Ext);
		s_DoneOpens.push_back(Open);
		s_DoneExts.push_back(Ext);
		s_DoneSum3s.push_back(Sum3);

		bool Better = false;
		if (Sum3 > s_BestSum3)
			{
			Better = true;
			s_BestSum3 = Sum3;
			s_BestOpen = Open;
			s_BestExt = Ext;

			AddPendingIfOk(Open+1, Ext);
			AddPendingIfOk(Open-1, Ext);
			AddPendingIfOk(Open, Ext+1);
			AddPendingIfOk(Open, Ext-1);

			ProgressLog("\n >>> [%.3f]   open/ext %d/%d scale %d\n\n",
				s_BestSum3, Open, Ext, ScaleFactor);
			}
		else
			ProgressLog("\n ... [%.3f]   open/ext %d/%d %.3f scale %d\n\n",
				s_BestSum3, Open, Ext, Sum3, ScaleFactor);

		if (!Better && feq(Sum3, s_BestSum3))
			{
			AddPendingIfOk(Open+1, Ext);
			AddPendingIfOk(Open-1, Ext);
			AddPendingIfOk(Open, Ext+1);
			AddPendingIfOk(Open, Ext-1);
			}
		}

	if (s_PendingOpens.empty())
		ProgressLog("Optimize converged\n");
	else
		Warning("Failed to converge\n");
	}

// Optimize gap penalies on SCOP40 for Nu given scale factor and 
// fixed (pre-optimized) weights for components
// H-J-like hack for integers
void cmd_hjnugaps()
	{
	GetFeatures(g_Arg1, s_Fs, s_Weights);

	s_PS = new ParaSearch;
	s_PS->m_NuFs = s_Fs;
	s_PS->GetByteSeqs(opt(db), "nuletters");
	s_PS->SetLookupFromLabels();

	double BestSum3 = -999;
	int BestOpen = -999;
	int BestExt = -999;
	int BestScaleFactor = -999;
	int MaxScaleFactor = 3;
	int ScaleFactor = 1;
	for (;;)
		{
		int FirstOpen = ScaleFactor*3;
		int FirstExt = ScaleFactor*3;
		ProgressLog("Optimize open=%d, ext=%d, scale=%d\n",
			FirstOpen, FirstExt, ScaleFactor);
		Optimize(ScaleFactor, FirstOpen, FirstExt);

		if (s_BestSum3 > BestSum3)
			{
			BestSum3 = s_BestSum3;
			BestOpen = s_BestOpen;
			BestExt = s_BestExt;
			BestScaleFactor = ScaleFactor;
			if (ScaleFactor >= 3)
				MaxScaleFactor++;
			}
		if (ScaleFactor == MaxScaleFactor)
			break;
		++ScaleFactor;
		}

	ProgressLog("\n");
	ProgressLog("Best Sum3=%.3f open=%d ext=%d scale=%d bits=%d\n",
		BestSum3, BestOpen, BestExt, BestScaleFactor, Paralign::m_Bits);
	ProgressLog("\n");
	}
