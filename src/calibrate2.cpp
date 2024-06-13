#include "myutils.h"
#include "scop40bench.h"
#include "binner.h"

static const uint BinCount = 100;

void cmd_calibrate2()
	{
	asserta(optset_benchlevel);

	string CalFN;
	if (g_Arg1 == ".")
#ifdef _MSC_VER
		CalFN = "c:/src/reseek_scop40/reseek_db/scop40_family.cal";
#else
		CalFN = "/c/src/reseek_scop40/reseek_db/scop40_family.cal";
#endif
	else
		CalFN = g_Arg1;


	SCOP40Bench SB;
	SB.ReadChains(CalFN, "");

	DSSParams Params;
	Params.SetFromCmdLine(SB.GetDBSize());

	SB.Setup(Params);

	SB.m_QuerySelf = true;
	SB.m_ScoresAreEvalues = true;
	SB.Run();
	SB.SetTFs();
	const float MaxFPR = 0.1;
	SB.SetStats(MaxFPR, true);
	//ROCStepsToTsv(opt_roc, m_SmoothScores, m_SmoothNTPs, m_SmoothNFPs,
	//  m_SmoothTPRs, m_SmoothFPRs);
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
	uint NQ = SB.m_QueryChainCount;
	uint TotalAlns = NQ*NQ;

	if (optset_output)
		{
		FILE *f = CreateStdioFile(opt_output);
		fprintf(f, "TS\tP\n");
		for (uint i = 0; i < N; ++i)
			{
			float TS = Scores[i];
			uint NFP = NFPs[i];
			if (NFP < NQ/100)
				continue;
			if (NFP > NQ*100)
				break;
		// P is probability hit is FP if TS >= thisTS
			double P = double(NFP)/TotalAlns;
			fprintf(f, "%.4g\t%.4g\n", TS, P);
			}
		}
	}
