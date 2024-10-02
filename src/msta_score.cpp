#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"

void cmd_msta_score()
	{
	asserta(optset_input);

	const bool DoCore = opt_core;

	string Name;
	GetStemName(g_Arg1, Name);

	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);

	FILE* fOut = CreateStdioFile(opt_output);
	const bool MissingSeqOk = opt_missingtestseqok;

	DALIScorer DS;
	DS.LoadChains(opt_input);
	bool Ok = DS.SetMSA(Name, MSA, DoCore, MissingSeqOk);
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
				Pf(fOut, "%s\t%s\tERROR_structure_not_found\n", Label1, Label2);
				continue;
				}
			++PairCount;

			double LDDT_mu = DS.GetLDDTChainPair_muscle(ChainIdx1, ChainIdx2, ColToPos1, ColToPos2);
			double Score_notused;
			double Z = DBL_MAX;
			double Z15 = DBL_MAX;

			DS.m_DALI_R0 = DBL_MAX;
			bool Ok = DS.GetDALIRowPair(SeqIdx1, SeqIdx2, Score_notused, Z);
			asserta(Ok);

			DS.m_DALI_R0 = 15;
			Ok = DS.GetDALIRowPair(SeqIdx1, SeqIdx2, Score_notused, Z15);
			asserta(Ok);
			DS.m_DALI_R0 = DBL_MAX;

			Sum_Z += Z;
			Sum_Z15 += Z15;
			Sum_LDDT_mu += LDDT_mu;
			Pf(fOut, "label1=%s\tlabel2=%s\tLDDT_mu=%.4f\tZ=%.3f\tZ15=%.3f\n",
			  Label1, Label2, LDDT_mu, Z, Z15);
			}
		}

	double LDDT_fm = DS.GetLDDT_foldmason();
	double Mean_Z = 0;
	double Mean_Z15 = 0;
	double Mean_LDDT_mu = 0;
	if (PairCount > 0)
		{
		Mean_Z = Sum_Z/PairCount;
		Mean_Z15 = Sum_Z15/PairCount;
		Mean_LDDT_mu = Sum_LDDT_mu/PairCount;
		}
	Pf(fOut, "MSA=%s", Name.c_str());
	Pf(fOut, "\tLDDT_fm=%.4f", LDDT_fm);
	Pf(fOut, "\tavg_LDDT_mu=%.4f", Mean_LDDT_mu);
	Pf(fOut, "\tavg_Z=%.3f", Mean_Z);
	Pf(fOut, "\tavg_Z15=%.3f", Mean_Z15);
	Pf(fOut, "\n");
	CloseStdioFile(fOut);

	ProgressLog("MSA=%s", Name.c_str());
	ProgressLog("\tLDDT_fm=%.4f", LDDT_fm);
	ProgressLog("\tavg_LDDT_mu=%.4f", Mean_LDDT_mu);
	ProgressLog("\tavg_Z=%.3f", Mean_Z);
	ProgressLog("\tavg_Z15=%.3f", Mean_Z15);
	ProgressLog("\n");
	}
