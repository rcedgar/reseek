#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"

void cmd_lddt_msa()
	{
	asserta(optset_input);
	const bool MissingTestSeqOk = !opt_missingtestseqok;

	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);
	FILE *fOut = CreateStdioFile(opt_output);

	string Name;
	GetStemName(g_Arg1, Name);

	bool DoCore = opt_core;

	DALIScorer DS;
	DS.LoadChains(opt_input);
	DS.SetMSA(Name, MSA, DoCore, MissingTestSeqOk);

	const uint SeqCount = MSA.GetSeqCount();
	double SumLDDT = 0;
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
			double PairLDDT = DS.GetLDDTPair_muscle(ChainIdx1, ChainIdx2, ColToPos1, ColToPos2);
			++PairCount;
			SumLDDT += PairLDDT;
			Pf(fOut, "%s\t%s\t%.4f\n", Label1, Label2, PairLDDT);
			}
		}

	double LDDT = 0;
	if (PairCount > 0)
		LDDT = SumLDDT/PairCount;
	ProgressLog("LDDT=%.4f MSA=%s\n", LDDT, Name.c_str());
	Pf(fOut, "LDDT=%.4f\tMSA=%s\n", LDDT, Name.c_str());
	CloseStdioFile(fOut);
	}
