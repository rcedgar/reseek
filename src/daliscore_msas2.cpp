#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"

void cmd_daliscore_msas2()
	{
	asserta(optset_input);
	asserta(optset_testdir);
	asserta(optset_testdir2);

	const bool MissingSeqOk = !opt_missingtestseqok;
	const bool DoCore = opt_core;
	FILE* fOut = CreateStdioFile(opt_output);

	string TestDir1 = string(opt_testdir);
	string TestDir2 = string(opt_testdir2);
	Dirize(TestDir1);
	Dirize(TestDir2);

	DALIScorer DS;
	DS.LoadChains(opt_input);

	vector<string> FNs;
	ReadLinesFromFile(g_Arg1, FNs);

	const uint N = SIZE(FNs);
	double SumZ = 0;
	double SumZ_core = 0;
	double MeanZ = 0;
	double MeanZ_core = 0;
	uint N1 = 0;
	uint N2 = 0;
	uint Ntie = 0;
	double SumZ1 = 0;
	double SumZ2 = 0;
	double Sum1 = 0;
	double Sum2 = 0;
	for (uint i = 0; i < N; ++i)
		{
		const string &FN = FNs[i];

		SeqDB MSA1;
		SeqDB MSA2;
		MSA1.FromFasta(TestDir1 + FN, true);
		MSA2.FromFasta(TestDir2 + FN, true);
		const uint SeqCount = MSA1.GetSeqCount();
		asserta(MSA2.GetSeqCount() == SeqCount);

		DS.SetMSA(FN, MSA1, DoCore, MissingSeqOk);
		double Score1 = DS.GetSumScore_Rows();
		double Z1 = DS.GetZ();
		SumZ1 += Z1;

		DS.SetMSA(FN, MSA2, DoCore, MissingSeqOk);
		double Score2 = DS.GetSumScore_Rows();
		double Z2 = DS.GetZ();
		SumZ2 += Z2;

		if (Score1 == Score2)
			++Ntie;
		else if (Score1 > Score2)
			{
			// gt/lt not preserved because of averaging
			//if (Z1 < Z2)
			//	Die("%.1f %.1f %.4f %.4f", Score1, Score2, Z1, Z2);
			++N1;
			}
		else if (Score2 > Score1)
			{
			// gt/lt not preserved because of averaging
			//if (Z2 < Z1)
			//	Die("%.1f %.1f %.4f %.4f", Score1, Score2, Z1, Z2);
			++N2;
			}
		else
			asserta(false);
		
		if (Score1 < 0)
			{
			Score1 = 0;
			Score2 -= Score1;
			}
		if (Score2 < 0)
			{
			Score2 = 0;
			Score1 -= Score1;
			}
		double Norm1 = Score1/(Score1 + Score2 + 1);
		double Norm2 = Score2/(Score1 + Score2 + 1);
		Sum1 += Norm1;
		Sum2 += Norm2;

		Log("%12.1f %12.1f", Score1, Score2);
		Log(" %8.3g %8.3g\n",Norm1, Norm2);
		if (fOut != 0)
			{
			fprintf(fOut, "aln=%s", FN.c_str());
			fprintf(fOut, "\tscore1=%.1f", Score1);
			fprintf(fOut, "\tscore2=%.1f", Score2);
			fprintf(fOut, "\tz1=%.1f", Z1);
			fprintf(fOut, "\tz2=%.1f", Z2);
			fprintf(fOut, "\tz2=%.1f", Z2);
			fprintf(fOut, "\tnorm1=%.1f", Score1);
			fprintf(fOut, "\tnorm2=%.1f", Score2);
			fprintf(fOut, "\n");
			}
		ProgressStep(i, N, "%s  N1 %u N2 %u tie %u avg1 %.4f avg2 %.4g",
		  FN.c_str(), N1, N2, Ntie, Sum1/(i+1), Sum2/(i+1));
		}
	  
	if (fOut != 0)
		{
		fprintf(fOut, "testdir1=%s", TestDir1.c_str());
		fprintf(fOut, "\ttestdir2=%s", TestDir2.c_str());
		fprintf(fOut, "\tn1better=%u", N1);
		fprintf(fOut, "\tn2better=%u", N2);
		fprintf(fOut, "\tntie=%u", Ntie);
		fprintf(fOut, "\tavg1=%.8f", Sum1/N);
		fprintf(fOut, "\tavg2=%.8f", Sum2/N);
		fprintf(fOut, "\tZ1=%.2f", SumZ1/N);
		fprintf(fOut, "\tZ2=%.2f", SumZ2/N);
		fprintf(fOut, "\n");
		CloseStdioFile(fOut);
		}
	}
