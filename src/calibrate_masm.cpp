#include "myutils.h"
#include "binner.h"

/***
Input:
	smith \
		-searchs ../out/masm.files \
		-input ../out/test_reverse.mega \
		-minscore 1 \
		-calibrate calibrate.tsv

# head -n3 calibrate.tsv | less -S
a.1.1.msta      4693    15.4    15.2    15.2    15.1    14.8    14.7    14.7    14.6   >
a.1.2.msta      6191    7.3     7.3     6.9     6.4     6.1     6.0     5.9     5.8    >
a.100.1.msta    2881    6.2     6.0     5.7     5.5     5.4     5.4     5.3     5.2    >
***/

// y = mx + b
void LinearFit(const vector<float> &xs, const vector<float> &ys,
  float &m, float &b);

static float g_MinScore = 1;
static float g_MaxScore = 30;
static uint g_BinCount = 32;
static bool g_DoLog2 = false;

static void FitBins(const string &Label, 
  const vector<uint> &AccumRevBins, float *ptr_m, float *ptr_b)
	{
	asserta(SIZE(AccumRevBins) == g_BinCount);
	const uint N = AccumRevBins[0];
	const uint Min_n = N/100 + 1;
	const float BinSize = (g_MaxScore - g_MinScore)/g_BinCount;
	vector<float> xs;
	vector<float> ys;
	for (uint Bin = 0; Bin < g_BinCount; ++Bin)
		{
		uint n = AccumRevBins[Bin];
		if (n < Min_n)
			break;
		float x = g_MinScore + BinSize*Bin + BinSize/2;
		float y = log10f(float(n));
		xs.push_back(x);
		ys.push_back(y);
		}

	float m, b;
	LinearFit(xs, ys, m, b);
	*ptr_m = m;
	*ptr_b = b;

#if 0
	const uint K = SIZE(xs);
	Log("%s m=%.3g b=%.3g ", Label.c_str(), m, b);
	for (uint k = 0; k < K; ++k)
		{
		uint n = AccumRevBins[k];
		float x = xs[k];
		float y = ys[k];
		float yhat = m*x + b;
		uint nhat = uint(powf(10, yhat));
		Log(" %u,%u", n, nhat);
		}
	Log("\n");
#endif
	}

void cmd_calibrate_masm()
	{
	if (optset_log2) g_DoLog2 = opt_log2;
	if (optset_minscore) g_MinScore = float(opt_minscore);
	if (optset_maxscore) g_MaxScore = float(opt_maxscore);
	if (optset_bins) g_BinCount = opt_bins;

	string Line;
	vector<string> Fields;
	FILE *fIn = OpenStdioFile(g_Arg1);
	FILE *fOut = CreateStdioFile(opt_output);
	FILE *fOut2 = CreateStdioFile(opt_output2);

	bool HdrDone = false;
	while (ReadLineStdioFile(fIn, Line))
		{
		Split(Line, Fields, '\t');
		const uint FieldCount = SIZE(Fields);
		string MasmLabel = Fields[0];
		size_t dotpos = MasmLabel.find(".masm");
		if (dotpos != string::npos)
			MasmLabel = MasmLabel.substr(0, dotpos);

		const uint n = StrToUint(Fields[1]);
		asserta(FieldCount == n + 2);
		vector<float> Scores;
		for (uint i = 0; i < n; ++i)
			{
			float Score = (float) StrToFloat(Fields[i+2]);
			if (g_DoLog2)
				Score = log2f(Score);
			Scores.push_back(Score);
			}

		Binner<float> *B =
		  new Binner(Scores, g_BinCount, g_MinScore, g_MaxScore);

		if (!HdrDone)
			{
			if (fOut != 0)
				{
				fprintf(fOut, "Bin");
				for (uint Bin = 0; Bin < g_BinCount; ++Bin)
					{
					float Lo = B->GetBinMid(Bin);
					fprintf(fOut, "\t%.1f", Lo);
					}
				fprintf(fOut, "\n");
				}
			if (fOut2 != 0)
				{
				fprintf(fOut2, "AccRevBin\tm\tb");
				for (uint Bin = 0; Bin < g_BinCount; ++Bin)
					{
					float Lo = B->GetBinMid(Bin);
					fprintf(fOut2, "\t%.1f", Lo);
					}
				fprintf(fOut2, "\n");
				}
			HdrDone = true;
			}

		const vector<uint32_t> &Bins = B->GetBins();
		if (fOut != 0)
			{
			fprintf(fOut, "%s", MasmLabel.c_str());
			for (uint i = 0; i < SIZE(Bins); ++i)
				fprintf(fOut, "\t%u", Bins[i]);
			fprintf(fOut, "\n");
			}

		vector<uint32_t> AccRevBins;
		B->GetAccumBinsReverse(AccRevBins);
		float m = FLT_MAX;
		float b = FLT_MAX;
		FitBins(MasmLabel, AccRevBins, &m, &b);
		if (fOut2 != 0)
			{
			if (isnan(m))
				m = 0;
			if (isnan(b))
				b = 0;
			fprintf(fOut2, "%s\t%.3g\t%.3g", 
			  MasmLabel.c_str(), m, b);
			for (uint i = 0; i < SIZE(AccRevBins); ++i)
				fprintf(fOut2, "\t%u", AccRevBins[i]);
			fprintf(fOut2, "\n");
			}

		delete B;
		}
	CloseStdioFile(fIn);
	CloseStdioFile(fOut);
	CloseStdioFile(fOut2);
	}
