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

	DALIScorer DS;
	DS.LoadChains(opt_input);
	DS.SetMSA(Name, MSA, true, MissingTestSeqOk);
	double Z = DS.GetZ();
	double Z2 = DS.GetZ_Rows();

	double Score = DS.GetSumScore_Rows();
	double Score2 = DS.GetSumScore_Cols();

	DS.SetMSA(Name, MSA, false, MissingTestSeqOk);
	double Z_core = DS.GetZ();
	double Z_core2 = DS.GetZ_Rows();

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
				Pf(fOut, "%s\t%s\t%.3g\n", Label1, Label2, Score);
			else
				Pf(fOut, "%s\t%s\tERROR\n", Label1, Label2);
			}
		}

	double Score_core = DS.GetSumScore_Rows();
	double Score_core2 = DS.GetSumScore_Cols();

	ProgressLog("Z %.1f Z_core %.1f Score %.1f Score_core %.1f %s\n", 
	  Z, Z_core, Score, Score_core, Name.c_str());

	asserta(feq(Z, Z2));
	asserta(feq(Score, Score2));
	asserta(feq(Z_core, Z_core2));
	asserta(feq(Score_core, Score_core2));
	CloseStdioFile(fOut);
	}
