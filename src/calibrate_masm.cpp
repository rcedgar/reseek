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

void cmd_calibrate_masm()
	{
	bool DoLog2 = optset_log2;
	float MinScore = 5; if (optset_minscore) MinScore = float(opt_minscore);
	float MaxScore = 30;	if (optset_maxscore) MaxScore = float(opt_maxscore);
	uint BinCount = 16; if (optset_bins) BinCount = opt_bins;

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
		const string &MasmLabel = Fields[0];
		const uint n = StrToUint(Fields[1]);
		asserta(FieldCount == n + 2);
		vector<float> Scores;
		for (uint i = 0; i < n; ++i)
			{
			float Score = (float) StrToFloat(Fields[i+2]);
			if (DoLog2)
				Score = log2f(Score);
			Scores.push_back(Score);
			}

		Binner<float> *B =
		  new Binner(Scores, BinCount, MinScore, MaxScore);

		if (!HdrDone)
			{
			if (fOut != 0)
				{
				fprintf(fOut, "Bin");
				for (uint Bin = 0; Bin < BinCount; ++Bin)
					{
					float Lo = B->GetBinMid(Bin);
					fprintf(fOut, "\t%.1f", Lo);
					}
				fprintf(fOut, "\n");
				}
			if (fOut2 != 0)
				{
				fprintf(fOut2, "AccRevBin");
				for (uint Bin = 0; Bin < BinCount; ++Bin)
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
		if (fOut2 != 0)
			{
			fprintf(fOut2, "%s", MasmLabel.c_str());
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
