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

	vector<string> Accs;
	ReadLinesFromFile(g_Arg1, Accs);

	const bool DoCore = opt_core;

	string SeqsDir;
	if (optset_seqsdir)
		{
		SeqsDir = string(opt_seqsdir);
		Dirize(SeqsDir);
		}

	const uint N = SIZE(Accs);
	double Sum_Z = 0;
	double Sum_LDDT_mu = 0;
	uint FoundFileCount = 0;
	for (uint i = 0; i < N; ++i)
		{
		const string &Acc = Accs[i];
		ProgressStep(i, N, "%s", Acc.c_str());

		const string &FN = TestDir + Acc;
		if (!StdioFileExists(FN))
			{
			fprintf(fOut, "missing_aln=%s\n", FN.c_str());
			continue;
			}

		SeqDB MSA;
		if (optset_seqsdir)
			{
			string SeqsFN = SeqsDir + Acc;
			SeqDB EvalSeqs;
			EvalSeqs.FromFasta(SeqsFN, false);
			MSA.FromFasta_Seqs(FN, EvalSeqs, true);
			}
		else
			MSA.FromFasta(FN, true);
		if (MSA.GetSeqCount() == 0)
			{
			fprintf(fOut, "empty_aln=%s\n", FN.c_str());
			continue;
			}

		++FoundFileCount;

		bool Ok = DS.SetMSA(Acc, MSA, DoCore, MissingSeqOk);
		if (!Ok)
			continue;

		double Z = DS.GetZ();
		double LDDT_mu = DS.GetLDDT_muscle();

		uint SeqCount = DS.GetSeqCount();
		uint CoreColCount = DS.m_CoreColCount;

		Sum_Z += Z;
		Sum_LDDT_mu += LDDT_mu;
		if (fOut != 0)
			{
			fprintf(fOut, "aln=%s", FN.c_str());
			fprintf(fOut, "\tseqs=%u", SeqCount);
			fprintf(fOut, "\tZ=%.3f", Z);
			fprintf(fOut, "\tLDDT_mu=%.4f", LDDT_mu);
			if (DoCore)
				fprintf(fOut, "\tnr_core_cols=%u", CoreColCount);
			fprintf(fOut, "\n");
			}
		}

	double Mean_Z = 0;
	double Mean_LDDT_mu = 0;
	if (FoundFileCount > 0)
		{
		Mean_Z = Sum_Z/FoundFileCount;
		Mean_LDDT_mu = Sum_LDDT_mu/FoundFileCount;
		}

	if (fOut != 0)
		{
		fprintf(fOut, "testdir=%s", TestDir.c_str());
		fprintf(fOut, "\tavg_Z=%.4f", Mean_Z);
		fprintf(fOut, "\tavg_LDDT_mu=%.4f", Mean_LDDT_mu);
		fprintf(fOut, "\n");
		}
	  
	CloseStdioFile(fOut);

	ProgressLog("MSAs=%u/%u", FoundFileCount, N);
	ProgressLog(" Z=%.3f", Mean_Z);
	//ProgressLog(" Z15=%.3f", Mean_Z15);
	ProgressLog(" LDDT_mu=%.4f", Mean_LDDT_mu);
	//ProgressLog(" LDDT_fm=%.4f", Mean_LDDT_fm);
	ProgressLog("\n");
	}
