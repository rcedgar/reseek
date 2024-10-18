#include "myutils.h"
#include "pdbchain.h"

static const double g_LDDT_R0 = 15;
static const double g_LDDT_thresholds[4] = { 0.5, 1, 2, 4 };
static const uint g_nr_thresholds = 4;

double GetLDDT_mu(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  bool DaliScorerCompatible)
	{
	const uint nr_cols = SIZE(PosQs);
	if (nr_cols == 0)
		return 0;
	asserta(SIZE(PosTs) == nr_cols);
	double total = 0;
	uint nr_cols_considered = 0;
	for (uint coli = 0; coli < nr_cols; ++coli)
		{
		uint pos1i = PosQs[coli];
		uint pos2i = PosTs[coli];
		if (pos1i == UINT_MAX || pos2i == UINT_MAX)
			continue;

		++nr_cols_considered;
		uint nr_considered = 0;
		uint nr_preserved = 0;
		for (uint colj = 0; colj < nr_cols; ++colj)
			{
			if (coli == colj)
				continue;
			uint pos1j = PosQs[colj];
			uint pos2j = PosTs[colj];
			if (pos1j == UINT_MAX || pos2j == UINT_MAX)
				continue;

			double d1 = Q.GetDist(pos1i, pos1j);
			double d2 = T.GetDist(pos2i, pos2j);
			if (DaliScorerCompatible)
				{
				if (d1 > g_LDDT_R0)
					continue;
				}
			else
				{
				if (d1 > g_LDDT_R0 && d2 > g_LDDT_R0)
					continue;
				}

			for (uint k = 0; k < g_nr_thresholds; ++k)
				{
				double t = g_LDDT_thresholds[k];
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

#if 0
#include "seqdb.h"
#include "daliscorer.h"

void cmd_test()
	{
	asserta(optset_input);

	string Name;
	GetStemName(g_Arg1, Name);

	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);

	FILE* fOut = CreateStdioFile(opt_output);
	const bool MissingSeqOk = opt_missingtestseqok;

	DALIScorer DS;
	DS.LoadChains(opt_input);
	bool Ok = DS.SetMSA(Name, MSA, false, MissingSeqOk);
	if (!Ok)
		Die("SetMSA failed");

	const uint SeqCount = MSA.GetSeqCount();
	double Sum_Z = 0;
	double Sum_Z15 = 0;
	double Sum_LDDT_mu = 0;
	double Sum_LDDT_fm = 0;

	uint PairCount = 0;
	for (uint SeqIdx1 = 0; SeqIdx1 < SeqCount; ++SeqIdx1)
		{
		const uint ChainIdx1 = DS.m_SeqIdxToChainIdx[SeqIdx1];
		const char *Label1 = MSA.GetLabel(SeqIdx1).c_str();
		const vector<uint> &ColToPos1 = DS.m_ColToPosVec[SeqIdx1];
		for (uint SeqIdx2 = SeqIdx1 + 1; SeqIdx2 < SeqCount; ++SeqIdx2)
			{
			const uint ChainIdx2 = DS.m_SeqIdxToChainIdx[SeqIdx2];
			const vector<uint> &ColToPos2 = DS.m_ColToPosVec[SeqIdx2];
			const char *Label2 = MSA.GetLabel(SeqIdx2).c_str();
			if (ChainIdx1 == UINT_MAX || ChainIdx2 == UINT_MAX)
				{
				Log("%s\t%s\tERROR_structure_not_found\n", Label1, Label2);
				continue;
				}
			++PairCount;

			double LDDT_mu = DS.GetLDDTChainPair_muscle(
			  ChainIdx1, ChainIdx2, ColToPos1, ColToPos2);

			const PDBChain &Chain1 = *DS.m_Chains[ChainIdx1];
			const PDBChain &Chain2 = *DS.m_Chains[ChainIdx2];
			double LDDT_mu2 = GetLDDT_mu(
			  Chain1, Chain2, ColToPos1, ColToPos2, true);
			double LDDT_mu3 = GetLDDT_mu(
			  Chain1, Chain2, ColToPos1, ColToPos2, false);

			Log("label1=%s\tlabel2=%s\tLDDT_mu=%.4f\tLDDT_mu2=%.4f\tLDDT_mu3=%.4f\n",
			  Label1, Label2, LDDT_mu, LDDT_mu2, LDDT_mu3);
			
			asserta(feq(LDDT_mu, LDDT_mu2));
			}
		}
	}
#endif // 0
