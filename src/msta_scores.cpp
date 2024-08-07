#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"

void cmd_msta_scores()
	{
	asserta(optset_input);
	asserta(optset_testdir);
	FILE* fOut = CreateStdioFile(opt_output);
	const bool MissingSeqOk = opt_missingtestseqok;

	string TestDir = string(opt_testdir);
	Dirize(TestDir);

	DALIScorer DS;
	DS.LoadChains(opt_input);

	vector<string> FNs;
	ReadLinesFromFile(g_Arg1, FNs);

	const bool DoCore = opt_core;

	const uint N = SIZE(FNs);
	double Sum_Z = 0;
	double Sum_Z15 = 0;
	double Sum_LDDT_mu = 0;
	double Sum_LDDT_fm = 0;
	for (uint i = 0; i < N; ++i)
		{
		const string &FN = FNs[i];

		SeqDB MSA;
		MSA.FromFasta(TestDir + FN, true);

		DS.SetMSA(FN, MSA, DoCore, MissingSeqOk);
		ProgressStep(i, N, "%s", DS.m_Name.c_str());

		double Z = DS.GetZ();
		double LDDT_mu = DS.GetLDDT_muscle();
		double LDDT_fm = DS.GetLDDT_foldmason();

		DS.m_DALI_R0 = 15;
		double Z15 = DS.GetZ();
		DS.m_DALI_R0 = DBL_MAX;

		uint CoreColCount = DS.m_CoreColCount;

		Sum_Z += Z;
		Sum_Z15 += Z15;
		Sum_LDDT_mu += LDDT_mu;
		Sum_LDDT_fm += LDDT_fm;
		if (fOut != 0)
			{
			fprintf(fOut, "aln=%s", FN.c_str());
			fprintf(fOut, "\tZ=%.3f", Z);
			fprintf(fOut, "\tZ15=%.3f", Z15);
			fprintf(fOut, "\tLDDT_mu=%.4f", LDDT_mu);
			fprintf(fOut, "\tLDDT_fm=%.4f", LDDT_fm);
			if (DoCore)
				fprintf(fOut, "\tnr_core_cols=%u", CoreColCount);
			fprintf(fOut, "\n");
			}
		}

	double Mean_Z = 0;
	double Mean_Z15 = 0;
	double Mean_LDDT_mu = 0;
	double Mean_LDDT_fm = 0;
	if (N > 0)
		{
		Mean_Z = Sum_Z/N;
		Mean_LDDT_mu = Sum_LDDT_mu/N;
		Mean_LDDT_fm = Sum_LDDT_fm/N;
		}

	if (fOut != 0)
		{
		fprintf(fOut, "testdir=%s", TestDir.c_str());
		fprintf(fOut, "\tavg_Z=%.1f", Mean_Z);
		fprintf(fOut, "\tavg_Z15=%.1f", Mean_Z15);
		fprintf(fOut, "\tavg_LDDT_mu=%.1f", Mean_LDDT_mu);
		fprintf(fOut, "\tavg_LDDT_fm=%.1f", Mean_LDDT_fm);
		fprintf(fOut, "\n");
		}
	  
	CloseStdioFile(fOut);
	}
