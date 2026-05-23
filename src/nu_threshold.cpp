#include "myutils.h"
#include "fastbench.h"

/***
Optimized Nu alphabet with weights
----------------------------------
	C:\src\2025-10_reseek_tune\bash\2026-04-01_hjnumega_parasail.bash
	aa4=0.481;pm2=0.301;sec32=0.219;

Threshold tuning
----------------
	2025-10_reseek_tune/2026-04-03_nu_threshold
	bash/2026-04-03_nu_threshold.bash

SEPQ0.1=0.296 SEPQ1=0.407 SEPQ10=0.491 Sum3=1.694 | verysensitive v2.7 without filter
SEPQ0.1=0.294 SEPQ1=0.406 SEPQ10=0.545 Sum3=1.742 | verysensitive v2.7 with nu filter

hresh    Sum3  Passed    %pass  Speedup   %sum3+
     1  1.7416    7.1M   20.32%     0.27     6.5%  |  SEPQ0.1=0.296 SEPQ1=0.411 SEPQ10=0.532 Sum3=1.742
    10  1.7417    4.9M   14.12%     0.38     6.5%  |  SEPQ0.1=0.296 SEPQ1=0.410 SEPQ10=0.534 Sum3=1.742
    20  1.7476    3.2M    9.30%     0.58     6.8%  |  SEPQ0.1=0.296 SEPQ1=0.414 SEPQ10=0.534 Sum3=1.748
    30  1.7416    2.1M    6.07%     0.89     6.5%  |  SEPQ0.1=0.296 SEPQ1=0.412 SEPQ10=0.532 Sum3=1.742
    40  1.7393    1.4M    3.98%     1.36     6.3%  |  SEPQ0.1=0.296 SEPQ1=0.413 SEPQ10=0.528 Sum3=1.739
    50  1.7259  920.8k    2.65%     2.05     5.5%  |  SEPQ0.1=0.295 SEPQ1=0.409 SEPQ10=0.521 Sum3=1.726
    60  1.7171  623.2k    1.79%     3.03     5.0%  |  SEPQ0.1=0.295 SEPQ1=0.409 SEPQ10=0.514 Sum3=1.717
    70  1.7616  432.6k    1.24%     4.37     7.7%  |  SEPQ0.1=0.294 SEPQ1=0.407 SEPQ10=0.563 Sum3=1.762
    72  1.7534  403.7k    1.16%     4.68     7.2%  |  SEPQ0.1=0.294 SEPQ1=0.405 SEPQ10=0.557 Sum3=1.753
    74  1.7498  377.5k    1.09%     5.00     7.0%  |  SEPQ0.1=0.294 SEPQ1=0.407 SEPQ10=0.551 Sum3=1.750
    76  1.7419  353.3k    1.02%     5.35     6.5%  |  SEPQ0.1=0.294 SEPQ1=0.406 SEPQ10=0.545 Sum3=1.742
    78  1.7343  331.1k    0.95%     5.70     6.0%  |  SEPQ0.1=0.294 SEPQ1=0.404 SEPQ10=0.539 Sum3=1.734
    80  1.7305  310.8k    0.89%     6.08     5.8%  |  SEPQ0.1=0.294 SEPQ1=0.406 SEPQ10=0.533 Sum3=1.731
    90  1.6955  231.7k    0.67%     8.15     3.7%  |  SEPQ0.1=0.294 SEPQ1=0.401 SEPQ10=0.507 Sum3=1.696
   100  1.6663  179.0k    0.51%    10.55     1.9%  |  SEPQ0.1=0.295 SEPQ1=0.396 SEPQ10=0.483 Sum3=1.666
   110  1.6423  145.1k    0.42%    13.02     0.4%  |  SEPQ0.1=0.295 SEPQ1=0.394 SEPQ10=0.461 Sum3=1.642
   125  1.5957  111.0k    0.32%    17.01    -2.5%  |  SEPQ0.1=0.295 SEPQ1=0.384 SEPQ10=0.429 Sum3=1.596
   150  1.5151   81.1k    0.23%    23.28    -7.4%  |  SEPQ0.1=0.294 SEPQ1=0.364 SEPQ10=0.382 Sum3=1.515
   200  1.3189   56.2k    0.16%    33.59   -19.4%  |  SEPQ0.1=0.282 SEPQ1=0.302 SEPQ10=0.302 Sum3=1.319

  mufilter  1.6358  n=1890353   5.43%
***/

void cmd_mu_threshold()
	{
	asserta(optset_lookup);
	asserta(optset_dope);
	const string &vsfn = g_Arg1;

	FastBench FB_vs, FB_t;
	FB_vs.ReadLookup(opt(lookup));
	FB_vs.ReadBits(vsfn);
	FB_t.ReadLookup(opt(lookup));
	FB_t.ReadDope(opt(dope));
	FB_t.Alloc();
	FB_t.m_scores_are_evalues = true;
	const uint npair = FB_vs.m_npair;
	for (uint i = 0; i < npair; ++i)
		FB_t.m_Scores[i] = FLT_MAX;

	uint n = 0;
	for (uint hitidx = 0; hitidx < FB_t.m_dope_nhit; ++hitidx)
		{
		uint k = FB_t.m_dope_ks[hitidx];
		asserta(k < npair);
		float score = FB_vs.m_Scores[k];
		FB_t.m_Scores[k] = score;
		n += 1;
		}
	FB_t.SetScoreOrder();
	double Sum3 = FB_t.Bench();
	ProgressLog("t=mufilter;");
	ProgressLog("Sum3=%.4f;", Sum3);
	ProgressLog("n=%u;", n);
	ProgressLog("pct=%.3g%%;", GetPct(n, npair));
	ProgressLog("\n");
	}

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

	const uint npair = FB_vs.m_npair;
	const vector<float> ts =
		{ 1, 10, 20, 30, 40, 50, 60, 70, 72, 74, 76, 78, 80, 90, 100, 110, 125, 150, 200 };

	ProgressLog("\n");
	ProgressLog("Thresh");
	//           123456
	ProgressLog("    Sum3");
	//             123456
	ProgressLog("  Passed");
	//             123456
	ProgressLog("    %%pass");
	//             1234567
	ProgressLog("  Speedup");
	//             1234567
	ProgressLog("   %%sum3+");
	//             1234567
	ProgressLog("\n");

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
		double Sum3 = FB_t.Bench("noshow");

	//mufilter  1.6358  n=1890353   5.43%
		double pct = GetPct(n, npair);
		double sum3 = Sum3;
		double speedx = 5.43/pct;
		double sum3x = sum3/1.6358;

		ProgressLog("%6.0f", t);
		ProgressLog("  %6.3f", sum3);
		ProgressLog("  %6.6s", IntToStr(n));
		ProgressLog("  %6.2f%%", pct);
		ProgressLog("  %7.2f", speedx);
		ProgressLog("  %6.1f%%", (sum3x - 1)*100);
		ProgressLog("  |  SEPQ0.1=%.3f", FB_t.m_SEPQ0_1);
		ProgressLog(" SEPQ1=%.3f", FB_t.m_SEPQ1);
		ProgressLog(" SEPQ10=%.3f", FB_t.m_SEPQ10);
		ProgressLog(" Sum3=%.3f", FB_t.m_Sum3);
		ProgressLog("\n");
		}
	}
