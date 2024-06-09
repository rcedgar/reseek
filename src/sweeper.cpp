#include "myutils.h"
#include "sweeper.h"
#include "sort.h"

uint Sweeper::Run(const DSSParams &Params, const string &Why)
	{
	asserta(m_SB != 0);
	m_SB->m_Params = &Params;
	m_SB->ClearHits();
	m_SB->Run();
	uint Score = m_SB->GetSens1stFP();

	if (m_fFev != 0)
		{
		fprintf(m_fFev, "%4u\twhy=%s", Score, Why.c_str());
		m_SB->m_Params->ToFev(m_fFev, false);
		}

	if (m_FirstScore == UINT_MAX)
		{
		m_FirstScore = Score;
		m_BestScore = Score;
		Log("%4u\twhy=%s", Score, Why.c_str());
		Params.ToFev(g_fLog, false);
		return Score;
		}

	double PctImproveVsPrevBest = 0;
	double d_first = double(m_BestScore) - double(m_FirstScore);
	double d_prevbest = double(Score) - double(m_BestScore);

	double PctBetter_first = m_FirstScore == 0 ? 0 : 100.0*d_first/m_FirstScore;
	double PctBetter_prevbest = m_BestScore == 0 ? 0 : 100.0*d_prevbest/m_BestScore;

	if (Score > m_BestScore)
		{
		m_BestScore = Score;
		if (m_fFev != 0)
			fprintf(m_fFev, "\timprove=%.1f%%", PctBetter_first);
		}

	if (PctBetter_prevbest > 0.1)
		{
		ProgressLog(" +++ ");
		if (m_fFev != 0)
			fprintf(m_fFev, "\tplusdelta=+%.3g%%", PctBetter_prevbest);
		}

	ProgressLog("%s Sens %u [%+%.0f] (%+.2f%%)\n",
	  Why.c_str(), Score, d_prevbest, PctBetter_first);

	return Score;
	}
