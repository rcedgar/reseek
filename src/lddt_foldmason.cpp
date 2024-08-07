#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"

double DALIScorer::GetLDDT_foldmason() const
	{
	uint nr_cols = GetColCount();
	double total_col_scores = 0;
	for (uint col = 0; col < nr_cols; ++col)
		{
		double col_score = GetLDDTColScore_foldmason(col);
		total_col_scores += col_score;
		}
	double mean_col_score = 0;
	if (nr_cols > 0)
		mean_col_score = total_col_scores/nr_cols;
	return mean_col_score;
	}

double DALIScorer::GetLDDTColScore_foldmason(uint col) const
	{
	double total_pair_scores = 0;
	uint nr_seq_pairs = 0;
	uint nr_seqs = GetSeqCount();
	uint nr_cols = GetColCount();
	for (uint seq_idxi = 0; seq_idxi < nr_seqs; ++seq_idxi)
		{
		uint posi = m_ColToPosVec[seq_idxi][col];
		if (posi == UINT_MAX)
			continue;
		asserta(seq_idxi < SIZE(m_DistMxVec));
		const vector<vector<double> > &dist_mxi = m_DistMxVec[seq_idxi];
		for (uint seq_idxj = seq_idxi+1; seq_idxj < nr_seqs; ++seq_idxj)
			{
			uint posj = m_ColToPosVec[seq_idxj][col];
			if (posj == UINT_MAX)
				continue;
			nr_seq_pairs += 1;
			const vector<vector<double> > &dist_mxj = m_DistMxVec[seq_idxj];
			double pair_score = 0;

			uint nr_pairs = 0;
			double score = 0;
			for (uint col2 = 0; col2 < nr_cols; ++col2)
				{
				if (col2 == col)
					continue;
				uint posi2 = m_ColToPosVec[seq_idxi][col2];
				uint posj2 = m_ColToPosVec[seq_idxj][col2];
				if (posi2 == UINT_MAX || posj2 == UINT_MAX)
					continue;
				
				double di = dist_mxi[posi][posi2];
				double dj = dist_mxj[posj][posj2];
				if (m_LDDT_symm == SYMM_First)
					{
					if (di > m_LDDT_R0)
						continue;
					}
				else if (m_LDDT_symm == SYMM_Both)
					{
					if (di > m_LDDT_R0 && dj > m_LDDT_R0)
						continue;
					}
				else if (m_LDDT_symm == SYMM_Either)
					{
					if (di > m_LDDT_R0 || dj > m_LDDT_R0)
						continue;
					}
				else
					asserta(false);

				nr_pairs += 1;
				double d_l = abs(di - dj);
				int isum = int(d_l < 0.5) + int(d_l < 1.0) + int(d_l < 2.0) + int(d_l < 4.0);
				score += double(isum)/4;
				}
			if (nr_pairs > 0)
				{
				pair_score = score/nr_pairs;
				total_pair_scores += pair_score;
				}
			}
		}
	if (nr_seq_pairs == 0)
		return 0;
	double col_score = total_pair_scores/nr_seq_pairs;
	return col_score;
	}
