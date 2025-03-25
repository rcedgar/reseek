#include "myutils.h"
#include "scop40bench.h"
#include "binner.h"
#include "output.h"

/***
* Calibrate distribution of FP errors on all-vs-all SCOP40 by
* linear fit of TS to -log(P) in range NFP=NQ/100 .. NFP=NQ*100.
* This works well, but similar fit for de novo per-chain
* calibration does not (calibrate3&4).
* 
* Superfamily: Linear fit to -log(P) m=20.5 b=2.89
* Fold:        Linear fit to -log(P) m=26.6 b=2.42
***/

static const uint BinCount = 100;

// y = mx + b
void LinearFit(const vector<float> &xs, const vector<float> &ys,
  float &m, float &b)
	{
	const uint N = SIZE(xs);
	float sumx = 0;
	float sumx2 = 0;
	float sumy = 0;
	float sumxy = 0;
	for (uint i = 0; i < N; ++i)
		{
		float x = xs[i];
		float y = ys[i];
		sumx += x;
		sumx2 += x*x;
		sumy += y;
		sumxy += x*y;
		}
	float meanx = sumx/N;
	float meany = sumy/N;

	float sumxd = 0;
	float sumxd2 = 0;
	float sumyd = 0;
	for (uint i = 0; i < N; ++i)
		{
		float xd = (meanx - xs[i]);
		float yd = (meany - ys[i]);
		sumxd += xd;
		sumxd2 += xd*xd;
		sumyd += yd;
		}
	m = (N*sumxy - sumx*sumy)/(N*sumx2 - sumx*sumx);
	b = meany - m*meanx;
	}

void cmd_calibrate2()
	{
	string CalFN;
	if (g_Arg1 == ".")
#ifdef _MSC_VER
		CalFN = "c:/src/reseek_scop40/reseek_db/scop40_family.cal";
#else
		CalFN = "/c/src/reseek_scop40/reseek_db/scop40_family.cal";
#endif
	else
		CalFN = g_Arg1;

	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast, SCOP40_DBSIZE);
	SCOP40Bench SB;
	SB.m_Params = &Params;
	SB.LoadDB(CalFN);

	asserta(SB.m_Params == &Params);
	Params.m_DBSize = (float) SB.GetDBSize();

	SB.Setup();
	
	float MaxFPR = 0.005f;
	if (optset_maxfpr)
		MaxFPR = (float) opt(maxfpr);

	OpenOutputFiles();

	SB.m_QuerySelf = true;
	SB.m_ScoresAreEvalues = true;
	if (opt(scores_are_not_evalues))
		SB.m_ScoresAreEvalues = false;	SB.RunSelf();
	SB.SetTFs();
	SB.SetStats(MaxFPR, true);

	const vector<float> &Scores = SB.m_SmoothScores;
	const vector<uint> &NTPs = SB.m_SmoothNTPs;
	const vector<uint> &NFPs = SB.m_SmoothNFPs;
	const vector<float> &TPRs = SB.m_SmoothTPRs;
	const vector<float> &FPRs = SB.m_SmoothFPRs;
	const uint N = SIZE(Scores);
	asserta(SIZE(TPRs) == N);
	asserta(SIZE(FPRs) == N);
	asserta(SIZE(NTPs) == N);
	asserta(SIZE(NFPs) == N);
	uint NQ = SB.GetDBChainCount();
	uint TotalAlns = NQ*NQ;

	vector<float> TSs;
	vector<float> Ps;
	vector<float> MinusLogPs;
	for (uint i = 0; i < N; ++i)
		{
		float TS = Scores[i];
		uint NFP = NFPs[i];
		if (NFP < NQ/100)
			continue;
		if (NFP > NQ*100)
			break;
	// P is probability hit is FP if TS >= thisTS
		float P = float(NFP)/TotalAlns;
		TSs.push_back(TS);
		Ps.push_back(P);
		float MinusLogP = -logf(P);
		MinusLogPs.push_back(MinusLogP);
		}
	const uint M = SIZE(TSs);
	float m, b;
	LinearFit(TSs, MinusLogPs, m, b);
	Log("Linear fit to -log(P) m=%.3g b=%.3g\n", m, b);

	if (optset_output)
		{
		FILE *f = CreateStdioFile(opt(output));
		fprintf(f, "TS\tP\tMinusLogP\tMinusLogP_fit\tP_fit\n");
		for (uint i = 0; i < M; ++i)
			{
			float TS = TSs[i];
			float P = Ps[i];
			float MinusLogP = MinusLogPs[i];
			float MinusLogP_fit = m*TS + b;
			float P_fit = expf(-MinusLogP_fit);
			fprintf(f, "%.4g\t%.4g\t%.4g\t%.4g\t%.4g\n",
			  TS, P, MinusLogP, MinusLogP_fit, P_fit);
			}
		}
	}
