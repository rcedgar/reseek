#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"

void cmd_msta_scores()
	{
	asserta(optset_input);
	asserta(optset_testdir);
	FILE* fOut = CreateStdioFile(opt_output);
	const bool MissingSeqOk = !opt_missingtestseqok;

	string TestDir = string(opt_testdir);
	Dirize(TestDir);

	DALIScorer DS;
	DS.LoadChains(opt_input);

	vector<string> FNs;
	ReadLinesFromFile(g_Arg1, FNs);

	const bool DoCore = opt_core;

	const uint N = SIZE(FNs);
	double SumZ = 0;
	for (uint i = 0; i < N; ++i)
		{
		const string &FN = FNs[i];

		SeqDB MSA;
		MSA.FromFasta(TestDir + FN, true);

		DS.SetMSA(FN, MSA, DoCore, MissingSeqOk);
		ProgressStep(i, N, "%s", DS.m_Name.c_str());

		double Z = DS.GetZ();
		double LDDT_muscle = DS.GetLDDT_muscle();
		double LDDT_foldmason = DS.GetLDDT_foldmason();
		uint CoreColCount = DS.m_CoreColCount;

		SumZ += Z;
		if (fOut != 0)
			{
			fprintf(fOut, "aln=%s", FN.c_str());
			fprintf(fOut, "\tZ=%.3f", Z);
			fprintf(fOut, "\tLDDT_mu=%.4f", LDDT_muscle);
			fprintf(fOut, "\tLDDT_fm=%.4f", LDDT_foldmason);
			if (DoCore)
				fprintf(fOut, "\tnr_core_cols=%u", CoreColCount);
			fprintf(fOut, "\n");
			}
		}

	double MeanZ = 0;
	if (N > 0)
		MeanZ = SumZ/N;

	if (fOut != 0)
		{
		fprintf(fOut, "testdir=%s", TestDir.c_str());
		fprintf(fOut, "\tZ=%.1f", MeanZ);
		fprintf(fOut, "\n");
		}
	  
	CloseStdioFile(fOut);
	}
