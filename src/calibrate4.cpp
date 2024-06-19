#include "myutils.h"
#include "scop40bench.h"
#include <map>

/***
* -calibrate4
* Unsupervised per-chain linear fit of -log(P) to TS
* P is probability randomly chosen chain will have score >= TS,
*   estimated by finding linear-like region of accumulated
*   count bins for TS values, fitting -log(P) = m*TS + b.
* Per-chain m and b written to -output slopes.tsv
* For search, slopes.tsv loaded into DBSearcher::m_ms, m_bs
* Gives much worse results for E-value.
* See -calibrate2 for supervised linear fit (good results).
* See -calibrate for unsupervised Gumbel fit (bad results).
***/

#define LOG_FIT	0

void LinearFit(const vector<float> &xs, const vector<float> &ys,
  float &m, float &b);

/***
-calibrate3 . -log calibrate3.log -calib_output calibrate3.tsv
less -S calibrate3.tsv
TS      -0.0938 -0.0812 -0.0688 -0.0562 -0.0437 -0.0312 -0.0187 -0.00625        0.00625 0.0188  0.0>
d1vkya_/e.53.1.1        11205   11205   11204   11197   11181   11134   11041   10745   10055   779>
d3nfka_/b.36.1.1        11205   11205   11205   11205   11205   11205   11202   11191   11158   109>

-scop40bench . -sens1fp_report sens1pf_report.tsv -log bench.log
# head sens1fp_report.tsv | columns.py
 d1vkya_/e.53.1                .      .      .   d2jdia2/b.49.1   0.096   86.2
 d3nfka_/b.36.1   d2i6va1/b.36.1  0.182   14.9    d2gj6d1/b.1.1   0.164   21.4
 d1t6ca2/c.55.1   d3ezwa1/c.55.1  0.128   44.6   d3kb2a_/c.37.1   0.128   44.9
d1v33a_/d.264.1  d1zt2a1/d.264.1  0.403   0.16   d1s12a_/d.64.2  0.0732    138
 d2gtlm1/b.61.7   d1jmxa5/b.61.4  0.362  0.367   d2a13a1/b.60.1   0.346  0.515
***/

// Bins are descending accumulated counts
// Find 50% top / 5% bottom percentiles
// this should look like log-linear FPs for
// fitting calibrated E-value for this DB chain.
static void GetBinLoHi(const vector<uint> &Bins, uint &BinLo, uint &BinHi)
	{
	BinLo = UINT_MAX;
	BinHi = UINT_MAX;
	const uint NBINS = SIZE(Bins);
	uint Max = Bins[0];
	for (uint Bin = 0; Bin < NBINS; ++Bin)
		{
		uint n = Bins[Bin];
		if (Bin > 0)
			asserta(n <= Bins[Bin-1]);
		if (BinLo == UINT_MAX && (n*100.0)/Max <= 50.0)
			BinLo = Bin;
		if (BinHi == UINT_MAX && (n*100.0)/Max <= 5.0)
			BinHi = Bin;
		}
	}

// Per-chain unsupervised log-linear fit
// Results are terrible, not sure why.
void cmd_calibrate4()
	{
	const string &FN1 = g_Arg1;
	const string &FN3 = opt_input3;

	string Line;
	vector<string> Fields;

	vector<string> Doms;
	map<string, uint> DomToIdx;
	vector<float> TSFP1s;

	FILE *f1 = OpenStdioFile(FN1);
	while (ReadLineStdioFile(f1, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 7);
		const string &Label = Fields[0];
		string Dom;
		SCOP40Bench::GetDomFromLabel(Label, Dom);
		const string &Fld5 = Fields[5];
		float TSFP1 = (Fld5 == "." ? FLT_MAX : (float) StrToFloat(Fields[5]));
		uint DomIdx = SIZE(Doms);
		Doms.push_back(Dom);
		TSFP1s.push_back(TSFP1);
		DomToIdx[Dom] = DomIdx;
		}
	CloseStdioFile(f1);
	const uint DomCount = SIZE(Doms);
	f1 = 0;

	FILE *f3 = OpenStdioFile(FN3);
	bool Ok = ReadLineStdioFile(f3, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) > 2);
	asserta(Fields[0] == "TS");
	const uint NBINS = SIZE(Fields) - 1;
	vector<float> TSHdr;
	for (uint i = 0; i < NBINS; ++i)
		{
		float TS = (float) StrToFloat(Fields[i+1]);
		TSHdr.push_back(TS);
		}

	vector<vector<uint> > BinsVec(DomCount);
	while (ReadLineStdioFile(f3, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == NBINS + 1);
		const string &Label = Fields[0];
		string Dom;
		SCOP40Bench::GetDomFromLabel(Label, Dom);

		map<string, uint>::const_iterator iter = DomToIdx.find(Dom);
		asserta(iter != DomToIdx.end());
		uint DomIdx = iter->second;
		asserta(DomIdx < DomCount);
		vector<uint> &Bins = BinsVec[DomIdx];
		for (uint i = 0; i < NBINS; ++i)
			{
			const string &Fld = Fields[i+1];
			uint n = (Fld == "" ? 0 : StrToUint(Fld));
			Bins.push_back(n);
			}
		}
	CloseStdioFile(f3);
	f3 = 0;

	vector<uint> BinLos;
	vector<uint> BinHis;
	vector<float> TSLos;
	vector<float> TSHis;
	vector<float> ms;
	vector<float> bs;
	float Sum_m = 0;
	float Sum_b = 0;
	uint Count = 0;
	for (uint i = 0; i < DomCount; ++i)
		{
		const string &Dom = Doms[i];
		float TSFP1 = TSFP1s[i];
		const vector<uint> Bins = BinsVec[i];
		uint BinLo, BinHi;
		GetBinLoHi(Bins, BinLo, BinHi);
		BinLos.push_back(BinLo);
		BinHis.push_back(BinHi);

		asserta(BinLo < SIZE(TSHdr));
		asserta(BinHi < SIZE(TSHdr));
		float TSLo = TSHdr[BinLo];
		float TSHi = TSHdr[BinHi];
		TSLos.push_back(TSLo);
		TSHis.push_back(TSHi);
		float m = FLT_MAX;
		float b = FLT_MAX;
		if (BinHi > BinLo)
			{
			uint Max = Bins[0];
			vector<float> TSs;
			vector<float> MinusLogPs;
			for (uint Bin = BinLo; Bin <= BinHi; ++Bin)
				{
				uint n = Bins[Bin];
				float P = float(n)/Max;
				float MinusLogP = (n == 0 ? 99.0f : -logf(P));
				float TS = TSLo + (Bin - BinLo)*(TSHi - TSLo)/(BinHi - BinLo);
				TSs.push_back(TS);
				MinusLogPs.push_back(MinusLogP);
				}
			LinearFit(TSs, MinusLogPs, m, b);
			++Count;
			Sum_m += m;
			Sum_b += b;
#if LOG_FIT
			{
			Log(">%s m=%.3g b=%.3g\n", Doms[i].c_str(), m, b);
			for (uint Bin = BinLo; Bin <= BinHi; ++Bin)
				{
				float TS = TSLo + (Bin - BinLo)*(TSHi - TSLo)/(BinHi - BinLo);
				uint n = Bins[Bin];
				float P = float(n)/Max;
				float PredMinusLogP = m*TS + b;
				Log(" TS=%.3g P=%.3g MinusLogP=%.3g PredMinusLogP=%.3g\n",
				  TSs[Bin-BinLo], P, MinusLogPs[Bin-BinLo], PredMinusLogP);
				}
			}
#endif
			}
		else
			{
#if LOG_FIT
			Log(">%s (no lo-hi)\n", Doms[i].c_str());
#endif
			}
		ms.push_back(m);
		bs.push_back(b);
		}
	ProgressLog("Mean m %.3g, b %.3g\n",
	  Sum_m/Count, Sum_b/Count);

	if (optset_report)
		{
		FILE *fOut = CreateStdioFile(opt_report);
		fprintf(fOut, "Dom\tTSFP1\tTSLo\tTSHi\tm\tb");
		for (uint Bin = 0; Bin < NBINS; ++Bin)
			fprintf(fOut, "\t%.3g", TSHdr[Bin]);
		fprintf(fOut, "\n");

		for (uint i = 0; i < DomCount; ++i)
			{
			const string &Dom = Doms[i];
			float TSFP1 = TSFP1s[i];
			const vector<uint> Bins = BinsVec[i];
			uint BinLo = BinLos[i];
			uint BinHi = BinHis[i];
			float TSLo = TSLos[i];
			float TSHi = TSHis[i];
			float m = ms[i];
			float b = bs[i];
			fprintf(fOut, "%s\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g",
			  Dom.c_str(), TSFP1, TSLo, TSHi, m, b);
			for (uint Bin = 0; Bin < NBINS; ++Bin)
				fprintf(fOut, "\t%u", Bins[Bin]);
			fprintf(fOut, "\n");
			}
		CloseStdioFile(fOut);
		}

	if (optset_output)
		{
		FILE *fOut = CreateStdioFile(opt_output);
		fprintf(fOut, "Dom\tm\tb\n");

		for (uint i = 0; i < DomCount; ++i)
			{
			const string &Dom = Doms[i];
			float m = ms[i];
			float b = bs[i];
			fprintf(fOut, "%s\t%.3g\t%.3g\n", Dom.c_str(), m, b);
			}
		CloseStdioFile(fOut);
		}

	if (optset_output2)
		{
		FILE *fOut = CreateStdioFile(opt_output2);
		fprintf(fOut, "Dom\tm\tb");
		for (uint Bin = 0; Bin < NBINS; ++Bin)
			fprintf(fOut, "\t%.3g", TSHdr[Bin]);
		fprintf(fOut, "\n");

		for (uint i = 0; i < DomCount; ++i)
			{
			const string &Dom = Doms[i];
			float m = ms[i];
			float b = bs[i];
			fprintf(fOut, "%s\t%.3g\t%.3g", Dom.c_str(), m, b);
			for (uint Bin = 0; Bin < NBINS; ++Bin)
				{
				float TS = TSHdr[Bin];
				float MinusLogP = m*TS + b;
				float P = expf(-MinusLogP);
				fprintf(fOut, "\t%.3g", P);
				}
			fprintf(fOut, "\n");
			}
		CloseStdioFile(fOut);
		}
	}
