#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"

double DALIScorer::GetLDDT_muscle() const
	{
	const SeqDB &MSA = *m_MSA;
	const uint SeqCount = MSA.GetSeqCount();
	double SumLDDT = 0;
	uint PairCount = 0;
	for (uint SeqIdx1 = 0; SeqIdx1 < SeqCount; ++SeqIdx1)
		{
		const uint ChainIdx1 = m_SeqIdxToChainIdx[SeqIdx1];
		const char *Label1 = MSA.GetLabel(SeqIdx1).c_str();
		const vector<uint> &ColToPos1 = m_ColToPosVec[SeqIdx1];
		for (uint SeqIdx2 = SeqIdx1 + 1; SeqIdx2 < SeqCount; ++SeqIdx2)
			{
			const uint ChainIdx2 = m_SeqIdxToChainIdx[SeqIdx2];
			const vector<uint> &ColToPos2 = m_ColToPosVec[SeqIdx2];
			const char *Label2 = MSA.GetLabel(SeqIdx2).c_str();
			if (ChainIdx1 == UINT_MAX || ChainIdx2 == UINT_MAX)
				continue;
			double PairLDDT = 
			  GetLDDTChainPair_muscle(ChainIdx1, ChainIdx2, ColToPos1, ColToPos2);
			++PairCount;
			SumLDDT += PairLDDT;
			}
		}
	double LDDT = 0;
	if (PairCount > 0)
		LDDT = SumLDDT/PairCount;
	return LDDT;
	}

double DALIScorer::GetLDDTChainPair_muscle(uint ChainIdx1, uint ChainIdx2,
  const vector<uint> &col_to_pos1s, const vector<uint> &col_to_pos2s) const
	{
	const uint nr_cols = SIZE(col_to_pos1s);
	if (nr_cols == 0)
		return 0;

	const vector<vector<double> > &DistMx1 = GetDistMx(ChainIdx1);
	const vector<vector<double> > &DistMx2 = GetDistMx(ChainIdx2);

	const uint nr_thresholds = SIZE(m_LDDT_thresholds);
	asserta(SIZE(col_to_pos2s) == nr_cols);
	double total = 0;
	uint nr_cols_considered = 0;
	for (uint coli = 0; coli < nr_cols; ++coli)
		{
		if (m_DoCore && !m_ColIsCore[coli])
			continue;
		uint pos1i = col_to_pos1s[coli];
		uint pos2i = col_to_pos2s[coli];
		if (pos1i == UINT_MAX || pos2i == UINT_MAX)
			continue;

		const vector<double> &DistMxRow1i = DistMx1[pos1i];
		const vector<double> &DistMxRow2i = DistMx2[pos2i];

		++nr_cols_considered;
		uint nr_considered = 0;
		uint nr_preserved = 0;
		for (uint colj = 0; colj < nr_cols; ++colj)
			{
			if (coli == colj)
				continue;
			if (m_DoCore && !m_ColIsCore[colj])
				continue;
			uint pos1j = col_to_pos1s[colj];
			uint pos2j = col_to_pos2s[colj];
			if (pos1j == UINT_MAX || pos2j == UINT_MAX)
				continue;

			//double d1_old = GetDist(ChainIdx1, pos1i, pos1j);
			//double d2_old = GetDist(ChainIdx2, pos2i, pos2j);
			double d1 = DistMxRow1i[pos1j];
			double d2 = DistMxRow2i[pos2j];

			if (d1 > m_LDDT_R0)
				continue;

			for (uint k = 0; k < nr_thresholds; ++k)
				{
				double t = m_LDDT_thresholds[k];
				nr_considered += 1;
				double diff = abs(d1 - d2);
				if (diff <= t)
					nr_preserved += 1;
				}
			}
		double score = 0;
		if (nr_considered > 0)
			score = double(nr_preserved)/nr_considered;
		total += score;
		}

	if (nr_cols_considered == 0)
		return 0;
	double avg = total/nr_cols_considered;
	return avg;
	}
