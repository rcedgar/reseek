#include "myutils.h"
#include "binner.h"
#include "quarts.h"

void cmd_binner()
	{
	FILE *fIn = OpenStdioFile(g_Arg1);

	uint Bins = 32;
	uint FieldNr = 0;
	if (optset_fieldnr)
		{
		asserta(opt(fieldnr) > 0);
		FieldNr = opt(fieldnr) - 1;
		}
	if (optset_bins)
		Bins = opt(bins);
	float MinValue = -FLT_MAX;
	float MaxValue = FLT_MAX;
	if (optset_minval)
		MinValue = float(opt(minval));
	if (optset_maxval)
		MaxValue = float(opt(maxval));
	asserta(MinValue < MaxValue);

	vector<float> Values;
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(fIn, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) > FieldNr);
		float v = StrToFloatf(Fields[FieldNr]);
		if (opt(log10))
			{
			if (v < 1e-20)
				v = -20;
			else
				v = log10f(v);
			}
		Values.push_back(v);
		}

	FILE *fOut = CreateStdioFile(opt(output));
	FILE *fAccum = CreateStdioFile(opt(accum));
	FILE *fAccumRev = CreateStdioFile(opt(accumrev));
	QuartsFloat QF;
	GetQuartsFloat(Values, QF);
	QF.LogMe();
	QF.WriteMe(g_fLog);
	QF.WriteMe(stderr);

	if (MinValue == -FLT_MAX || MaxValue == FLT_MAX)
		{
		Binner<float> B(Values, Bins);
		B.ToTsv(fOut);
		B.AccumToTsv(fAccum);
		B.AccumToTsvReverse(fAccumRev);
		}
	else
		{
		Binner<float> B(Values, Bins, MinValue, MaxValue);
		B.ToTsv(fOut);
		B.AccumToTsv(fAccum);
		B.AccumToTsvReverse(fAccumRev);
		}
	CloseStdioFile(fOut);
	}
