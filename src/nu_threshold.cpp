#include "myutils.h"
#include "fastbench.h"

void cmd_nu_threshold()
	{
	asserta(optset_lookup);
	const string &nufn = g_Arg1;
	const string &vsfn = opt(input2);

	FastBench FB_nu, FB_vs, FB_t;
	FB_nu.ReadLookup(opt(lookup));
	FB_vs.ReadLookup(opt(lookup));
	FB_t.ReadLookup(opt(lookup));
	FB_t.Alloc();

	FB_nu.m_scores_are_evalues = false;
	FB_vs.m_scores_are_evalues = true;
	FB_t.m_scores_are_evalues = true;

	FB_nu.ReadBits(nufn);
	FB_vs.ReadBits(vsfn);

	ProgressLog("nu bench:\n");
	FB_nu.SetScoreOrder();
	FB_nu.Bench();
	ProgressLog("vs bench:\n");
	FB_vs.SetScoreOrder();
	FB_vs.Bench();

	const uint npair = FB_vs.m_PairCount;
	const vector<float> ts =
		{ 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
	for (auto t : ts)
		{
		uint n = 0;
		for (uint i = 0; i < npair; ++i)
			{
			if (FB_nu.m_Scores[i] >= t)
				{
				++n;
				FB_t.m_Scores[i] = FB_vs.m_Scores[i];
				}
			else
				FB_t.m_Scores[i] = FLT_MAX;
			}
		FB_t.SetScoreOrder();
		FB_t.Bench();
		ProgressLog("t=%.0f;", t);
		ProgressLog("Sum3=%.4f;", FB_t.m_Sum3);
		ProgressLog("n=%u;", n);
		ProgressLog("pct=%.3g%%;", GetPct(n, npair));
		ProgressLog("\n");
		}
	}
