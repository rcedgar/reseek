#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"

void cmd_daliscore_msa()
	{
	asserta(optset_input);
	const bool MissingTestSeqOk = !opt_missingtestseqok;

	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);
	FILE* fOut = CreateStdioFile(opt_output);

	string Name;
	GetStemName(g_Arg1, Name);

	DALIScorer DS;
	DS.LoadChains(opt_input);
	DS.SetMSA(Name, MSA, true, MissingTestSeqOk);
	double Z = DS.GetZ();
	double Z2 = DS.GetZ_Rows();
	asserta(feq(Z, Z2));

	DS.SetMSA(Name, MSA, false, MissingTestSeqOk);
	double Z_core = DS.GetZ();
	double Z_core2 = DS.GetZ_Rows();
	asserta(feq(Z_core, Z_core2));

	ProgressLog("Z %.1f Z_core %.1f %s\n", Z, Z_core, Name.c_str());
	}

void cmd_daliscore_msas()
	{
	asserta(optset_input);
	FILE* fOut = CreateStdioFile(opt_output);
	const bool MissingSeqOk = !opt_missingtestseqok;

	string TestDir = ".";
	if (optset_testdir)
		TestDir = string(opt_testdir);
	if (!EndsWith(TestDir, "/") && !EndsWith(TestDir, "\\"))
		TestDir += "/";

	DALIScorer DS;
	DS.LoadChains(opt_input);

	vector<string> FNs;
	ReadLinesFromFile(g_Arg1, FNs);

	const uint N = SIZE(FNs);
	double SumZ = 0;
	double SumZ_core = 0;
	double MeanZ = 0;
	double MeanZ_core = 0;
	for (uint i = 0; i < N; ++i)
		{
		const string &FN = FNs[i];

		SeqDB MSA;
		MSA.FromFasta(TestDir + FN, true);

		DS.SetMSA(FN, MSA, true, MissingSeqOk);
		double Z = DS.GetZ();
		double Z2 = DS.GetZ_Rows();

		DS.SetMSA(FN, MSA, false, MissingSeqOk);
		double Z_core = DS.GetZ();
		double Z_core2 = DS.GetZ_Rows();
		asserta(feq(Z, Z2));
		asserta(feq(Z_core, Z_core2));

		SumZ += Z;
		SumZ_core += Z_core;
		MeanZ = SumZ/(i+1);
		MeanZ_core = SumZ_core/(i+1);
		ProgressStep(i, N, "n %u Z %.2f Z_core %.2f", i+1, MeanZ, MeanZ_core);
		if (fOut != 0)
			{
			fprintf(fOut, "aln=%s", FN.c_str());
			fprintf(fOut, "\tZ=%.1f", Z);
			fprintf(fOut, "\tZ_core=%.1f", Z_core);
			fprintf(fOut, "\n");
			}
		}

	if (fOut != 0)
		{
		fprintf(fOut, "testdir=%s", TestDir.c_str());
		fprintf(fOut, "\tZ=%.1f", MeanZ);
		fprintf(fOut, "\tZ_core=%.1f", MeanZ_core);
		fprintf(fOut, "\n");
		}
	  
	CloseStdioFile(g_fLog);
	CloseStdioFile(fOut);
	}
