#include "myutils.h"
#include "binner.h"

// y = mx + b
void LinearFit(const vector<float> &xs, const vector<float> &ys,
  float &m, float &b);

static float g_MinScore;
static float g_MaxScore;
static uint g_BinCount;

static void FitBins(const string &Label,
  const vector<uint> &AccumRevBins,
  float *ptr_m, float *ptr_b)
	{
	*ptr_m = 0;
	*ptr_b = 0;
	asserta(SIZE(AccumRevBins) == g_BinCount);
	const uint N = AccumRevBins[0];
	const uint Max_n = N/2;
	const uint Min_n = Max_n/100;;
	const float BinSize = (g_MaxScore - g_MinScore)/g_BinCount;
	vector<uint> ns;
	vector<float> xs;
	vector<float> ys;
	for (uint Bin = 0; Bin < g_BinCount; ++Bin)
		{
		uint n = AccumRevBins[Bin];
		if (n > Max_n)
			continue;
		if (n < Min_n)
			break;
		float x = g_MinScore + BinSize*Bin + BinSize/2;
		float y = log10f(float(n));
		ns.push_back(n);
		xs.push_back(x);
		ys.push_back(y);
		}
	if (SIZE(xs) < 3)
		return;

	float m, b;
	LinearFit(xs, ys, m, b);
	*ptr_m = m;
	*ptr_b = b;

#if 0
	const uint K = SIZE(xs);
	Log("%s m=%.3g b=%.3g ", Label.c_str(), m, b);
	for (uint k = 0; k < K; ++k)
		{
		uint n = ns[k];
		float x = xs[k];
		float y = ys[k];
		float yhat = m*x + b;
		uint nhat = uint(powf(10, yhat));
		Log(" %u,%u", n, nhat);
		}
	Log("\n");
#endif
	}

/***
* Second step following calibrate5.
***/

void cmd_calibrate5a()
	{
	g_MinScore = (float) opt_minscore;
	g_MaxScore = (float) opt_maxscore;
	g_BinCount = opt_n;
	FILE *fIn = OpenStdioFile(g_Arg1);
	FILE *fOut = CreateStdioFile(opt_output);
	asserta(fOut != 0);

	string Line;
	vector<string> Fields;
	vector<float> Scores;
	Scores.reserve(100000000);
	uint LineCount = 0;
	bool HdrDone = false;
	while (ReadLineStdioFile(fIn, Line))
		{
		if (++LineCount%100 == 0)
			Progress("%u\r", LineCount);
		Split(Line, Fields, '\t');
		const uint n = SIZE(Fields);
		const string &Label = Fields[0];
		for (uint i = 1; i < n; ++i)
			{
			float Score = (float) StrToFloat(Fields[i]);
			Scores.push_back(Score);
			}
		Binner<float> B(Scores, g_BinCount, g_MinScore, g_MaxScore);
		vector<uint> AccumRevBins;
		B.GetAccumBinsReverse(AccumRevBins);
		asserta(SIZE(AccumRevBins) == g_BinCount);
		float m, b;
		FitBins(Label, AccumRevBins, &m, &b);
		if (!HdrDone)
			{
			fprintf(fOut, "Bin\tm\tb");
			for (uint i = 0; i < g_BinCount; ++i)
				fprintf(fOut, "\t%.3g", B.GetBinMid(i));
			fprintf(fOut, "\n");
			HdrDone = true;
			}
		fprintf(fOut, "%s\t%.3g\t%.3g", Label.c_str(), m, b);
		for (uint i = 0; i < g_BinCount; ++i)
			fprintf(fOut, "\t%u", AccumRevBins[i]);
		fprintf(fOut, "\n");
		Scores.clear();
		}
	Progress("%u\n", ++LineCount);
	CloseStdioFile(fIn);
	CloseStdioFile(fOut);
	}
