#include "myutils.h"
#include "seqdb.h"
#include "chainreader2.h"
#include "alpha.h"
#include "daliscorer.h"

double DALIScorer::GetDALIScore_ChainPair(uint ChainIdx1, uint ChainIdx2,
  const vector<uint> &Pos1s, vector<uint> &Pos2s) const
	{
	const uint Lali = SIZE(Pos1s);
	double OffDiag = GetDALIScore_OffDiag(ChainIdx1, ChainIdx2, Pos1s, Pos2s);
	double Score = OffDiag + Lali*m_DALI_Theta;
	return Score;
	}

double DALIScorer::GetDALIScore_OffDiag(uint ChainIdx1, uint ChainIdx2,
  const vector<uint> &Pos1s, const vector<uint> &Pos2s) const
	{
	const uint Lali = SIZE(Pos1s);
	asserta(SIZE(Pos2s) == Lali);

	const vector<vector<double> > &DistMx1 = GetDistMx(ChainIdx1);
	const vector<vector<double> > &DistMx2 = GetDistMx(ChainIdx2);

	double Sum = 0;
	for (uint i = 0; i < Lali; ++i)
		{
		uint Pos1i = Pos1s[i];
		uint Pos2i = Pos2s[i];
		const vector<double> &DistMxRow1i = DistMx1[Pos1i];
		const vector<double> &DistMxRow2i = DistMx2[Pos2i];
		for (uint j = 0; j < Lali; ++j)
			{
			if (i == j)
				continue;

			uint Pos1j = Pos1s[j];
			uint Pos2j = Pos2s[j];

			double dij_1 = DistMxRow1i[Pos1j];
			double dij_2 = DistMxRow2i[Pos2j];
			if (dij_1 > m_DALI_R0 || dij_2 > m_DALI_R0)
				continue;

			double x = DALI_dpscorefun(dij_1, dij_2);
			Sum += x;
			}
		}
	return Sum;
	}
