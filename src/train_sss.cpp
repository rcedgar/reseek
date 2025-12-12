#include "myutils.h"
#include "featuretrainer2.h"
#include "seqdb.h"
#include "sort.h"
#include "hexintseq.h"

void cmd_train_sss()
	{
	asserta(optset_background_style);
	asserta(optset_alpha_size);

	const string &SeqsFN = g_Arg1;
	const string &TrainTPAlnFN = opt(traintps); // "src/2025-10_reseek_tune/big_fa2/tp.mints05.maxts25.fa2";
	FILE *fOut = CreateStdioFile(opt(output));
	const uint AS = opt(alpha_size);

	FeatureTrainer2::m_AlphaSize = AS;
	FeatureTrainer2::m_FevStr.clear();

///////////////////////////////////////////////////////////////////////////////////////
// Structure chains, must include all Train and Eval chains, may include others
///////////////////////////////////////////////////////////////////////////////////////
	vector<vector<uint> > IntSeqs;
	map<string, uint> LabelToSeqIdx;
	Progress("Reading int seqs... ");
	ReadHexIntSeqs(AS, SeqsFN, IntSeqs, LabelToSeqIdx);
	Progress("done.\n");

	vector<bool> TrainsTPs_notused;
	vector<string> TrainRows;
	vector<string> TrainLabels;
	vector<uint> TrainSeqIdxs;
	FeatureTrainer2::AppendAlns("traintps", TrainTPAlnFN, LabelToSeqIdx, true,
		TrainRows, TrainLabels, TrainSeqIdxs, TrainsTPs_notused);

	vector<vector<float> > ScoreMx;
	BACKGROUND_STYLE BS = FeatureTrainer2::StrToBS(opt(background_style));
	FeatureTrainer2::TrainSSS(IntSeqs, TrainRows, TrainLabels,
		TrainSeqIdxs, ScoreMx, BS);

	FeatureTrainer2::ScoreMxToTsv(fOut, ScoreMx);
	Log("@FEV@ %s\n", FeatureTrainer2::m_FevStr.c_str());
	if (fOut != 0)
		{
		string CmdLine;
		GetCmdLine(CmdLine);
		fprintf(fOut, "cmd\t%s\n", CmdLine.c_str());
		fprintf(fOut, "fev\t%s\n", FeatureTrainer2::m_FevStr.c_str());
		fprintf(fOut, "git\t%s\n", g_GitVer);
		CloseStdioFile(fOut);
		}
	}
