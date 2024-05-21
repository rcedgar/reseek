#include "myutils.h"
#include "scop40bit.h"

uint SCOP40Bit::GetSens1stFP(const vector<float> &Scores) const
	{
	const uint DomCount = GetDomCount();
	const uint HitCount = GetHitCount();
	vector<float> DomToScoreFirstFP;
	if (m_ScoresAreEvalues)
		DomToScoreFirstFP.resize(DomCount, FLT_MAX);
	else
		DomToScoreFirstFP.resize(DomCount, 0);

	for (uint i = 0; i < HitCount; ++i)
		{
		uint Dom1 = m_DomIdx1s[i];
		uint Dom2 = m_DomIdx2s[i];
		float Score = Scores[i];
		uint Fam1 = m_DomIdxToFamIdx[Dom1];
		uint Fam2 = m_DomIdxToFamIdx[Dom2];
		if (Fam1 != Fam2)
			{
			if (m_ScoresAreEvalues)
				DomToScoreFirstFP[Dom1] = min(DomToScoreFirstFP[Dom1], Score);
			else
				DomToScoreFirstFP[Dom1] = max(DomToScoreFirstFP[Dom1], Score);
			}
		}

	uint GoodCount = 0;
	for (uint i = 0; i < HitCount; ++i)
		{
		uint Dom1 = m_DomIdx1s[i];
		uint Dom2 = m_DomIdx2s[i];
		float Score = Scores[i];
		uint Fam1 = m_DomIdxToFamIdx[Dom1];
		uint Fam2 = m_DomIdxToFamIdx[Dom2];
		if (Fam1 == Fam2)
			{
			if (m_ScoresAreEvalues)
				{
				if (Score < DomToScoreFirstFP[Dom1])
					++GoodCount;
				}
			else
				{
				if (Score > DomToScoreFirstFP[Dom1])
					++GoodCount;
				}
			}
		}
	return GoodCount;
	}
