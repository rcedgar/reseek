#include "myutils.h"
#include "sweeper.h"
#include "sort.h"

uint Sweeper::Run(const DSSParams &Params)
	{
	asserta(m_SB != 0);
	m_SB->m_Params = &Params;
	m_SB->ClearHits();
	m_SB->Run();
	uint Score = m_SB->GetSens1stFP();
	if (m_FirstScore == UINT_MAX)
		m_FirstScore = Score;
	int d = int(Score) - int(m_BestScore);
	bool Improve = (d > 0);
	if (Improve)
		m_BestScore = Score;
	int d2 = int(m_BestScore) - int(m_FirstScore);

	if (m_fFev != 0)
		{
		fprintf(m_fFev, "%4u\t", Score);
		m_SB->m_Params->ToFev(m_fFev);
		}

	Params.ToFev(g_fLog);
	ProgressLog(" Sens %u [%u %+d]", Score, m_BestScore, d);
	if (Improve)
		ProgressLog(" +++ ");
	ProgressLog(" (%+.2f%%)\n", GetPct(d2, m_FirstScore));

	//m_Scores.push_back(Score);
	return Score;
	}
