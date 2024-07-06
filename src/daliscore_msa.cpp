#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"

void cmd_daliscore_msa()
	{
	asserta(optset_input);

	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);
	FILE* fOut = CreateStdioFile(opt_output);

	DALIScorer DS;
	DS.LoadChains(opt_input);
	DS.Run(g_Arg1, MSA);
	DS.ToFev(g_fLog);
	DS.RunCols(g_Arg1, MSA, true);
	}

void cmd_daliscore_msas()
	{
	asserta(optset_input);
	FILE* fOut = CreateStdioFile(opt_output);

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
	double SumMeanScore = 0;
	double SumMeanScore_core = 0;
	double SumMeanZ_core = 0;
	double SumMeanZ = 0;
	uint n = 0;
	double MeanScore = 0;
	double MeanScore_core = 0;
	double MeanZ_core = 0;
	double MeanZ = 0;
	for (uint i = 0; i < N; ++i)
		{
		const string &FN = FNs[i];
		ProgressStep(i, N, "n %u score %.1f score_core %.1f Z %.2f Z_core %.2f",
		  n, MeanScore, MeanScore_core, MeanZ, MeanZ_core);

		SeqDB MSA;
		MSA.FromFasta(TestDir + FN, true);
		DS.Run(FN, MSA);
		DS.ToFev(fOut);
		if (DS.m_MeanScore == DBL_MAX)
			continue;
		n += 1;
		SumMeanScore += DS.m_MeanScore;
		SumMeanScore_core += DS.m_MeanScore_core;
		SumMeanZ_core += DS.m_MeanZ_core;
		SumMeanZ += DS.m_MeanZ;
		MeanScore = SumMeanScore/n;
		MeanScore_core = SumMeanScore_core/n;
		MeanZ_core = SumMeanZ_core/n;
		MeanZ = SumMeanZ/n;
		}
	ProgressLog("\nTotal n %u score %.1f score_core %.1f Z %.2f Z_core %.2f",
	  n, MeanScore, MeanScore_core, MeanZ, MeanZ_core);
	if (fOut != 0)
		{
		fprintf(fOut, "dir=%s", TestDir.c_str());
		ProgressLog("\nTotal n %u score %.1f score_core %.2f Z %.1f Z_core %.1f",
		  n, MeanScore, MeanScore_core, MeanZ, MeanZ_core);
		fprintf(fOut, "\tn=%u", n);
		fprintf(fOut, "\tfinal_score=%1.f", MeanScore);
		fprintf(fOut, "\tfinal_score_core=%1.f", MeanScore_core);
		fprintf(fOut, "\tfinal_Z=%.2f", MeanZ);
		fprintf(fOut, "\tfinal_Z_core=%.2f", MeanZ_core);
		fprintf(fOut, "\n");
		}
	  
	CloseStdioFile(g_fLog);
	CloseStdioFile(fOut);
	}
