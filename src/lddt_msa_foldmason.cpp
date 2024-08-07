#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"

void cmd_lddt_msa_foldmason()
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
	double LDDT = DS.GetLDDT_foldmason();
	ProgressLog("LDDT_fm=%.4f MSA=%s\n", LDDT, Name.c_str());
	Pf(fOut, "LDDT_fm=%.4f\tMSA=%s\n", LDDT, Name.c_str());
	CloseStdioFile(fOut);
	}
