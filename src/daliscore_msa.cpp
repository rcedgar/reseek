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
	FILE *fOut = CreateStdioFile(opt_output);

	string Name;
	GetStemName(g_Arg1, Name);

	bool DoCore = opt_core;

	DALIScorer DS;
	DS.LoadChains(opt_input);
	DS.SetMSA(Name, MSA, DoCore, MissingTestSeqOk);
	double Z = DS.GetZ();
	double Z2 = DS.GetZ_Rows();

	double Score = DS.GetSumScore_Rows();
	double Score2 = DS.GetSumScore_Cols();

	const uint SeqCount = MSA.GetSeqCount();
	double SumScore = 0;
	for (uint SeqIdx1 = 0; SeqIdx1 < SeqCount; ++SeqIdx1)
		{
		const char *Label1 = MSA.GetLabel(SeqIdx1).c_str();
		for (uint SeqIdx2 = SeqIdx1 + 1; SeqIdx2 < SeqCount; ++SeqIdx2)
			{
			const char *Label2 = MSA.GetLabel(SeqIdx2).c_str();
			double Score, Z;
			bool Ok = DS.GetDALIRowPair(SeqIdx1, SeqIdx2, Score, Z);
			if (Ok)
				Pf(fOut, "%s\t%s\t%.3g\t%.1f\n", Label1, Label2, Score, Z);
			else
				Pf(fOut, "%s\t%s\tERROR\n", Label1, Label2);
			}
		}

	ProgressLog("Z=%.1f Score=%.1f MSA=%s\n", Z, Score, Name.c_str());
	Pf(fOut, "Z=%.1f\tScore=%.1f\tMSA=%s\n", Z, Score, Name.c_str());

	asserta(feq(Z, Z2));
	asserta(feq(Score, Score2));
	CloseStdioFile(fOut);
	}
