#include "myutils.h"
#include "featuretrainer2.h"
#include "parasearch.h"
#include "peaker.h"

static ParaSearch *s_PS;
static Peaker *s_Peaker;

static double EvalSum3(const vector<string> &xv)
	{
	asserta(SIZE(xv) == 2);

	const float Open = StrToFloatf(xv[0]);
	const float Ext = StrToFloatf(xv[1]);

	const int IntOpen = int(round(Open));
	const int IntExt= int(round(Ext));

	s_PS->SetGapParams(IntOpen, IntExt);
	s_PS->ClearHitsAndResults();
	s_PS->Search("para");
	s_PS->Bench();
	return s_PS->m_SB.m_Sum3;
	}

static void TrainMx(vector<vector<int> > &IntScoreMx, double ScaleFactor)
	{
	const string &ChainFN = opt(db);			// "src/reseek/test_data/scop40.bca";
	const string &TrainTPAlnFN = opt(traintps); // "src/2025-10_reseek_tune/big_fa2/tp.mints05.maxts25.fa2";

	FeatureTrainer2::m_FevStr.clear();

	vector<PDBChain *> Chains;
	map<string, uint> LabelToChainIdx;
	FeatureTrainer2::ReadChains(ChainFN, Chains, LabelToChainIdx);

///////////////////////////////////////////////////////////////////////////////////////
// Training alignments (TP only)
///////////////////////////////////////////////////////////////////////////////////////
	vector<bool> TrainsTPs_notused;
	vector<string> TrainRows;
	vector<string> TrainLabels;
	vector<uint> TrainChainIdxs;
	FeatureTrainer2::AppendAlns("traintps", TrainTPAlnFN, LabelToChainIdx, true,
		TrainRows, TrainLabels, TrainChainIdxs, TrainsTPs_notused);

///////////////////////////////////////////////////////////////////////////////////////
// Eval alignments not used, empty vectors needed for TrainDSSFeature
///////////////////////////////////////////////////////////////////////////////////////
	vector<bool> _notused_EvalTPs;
	vector<string> _notused_EvalRows;
	vector<string> _notused_EvalLabels;
	vector<uint> _notused_EvalRowChainIdxs;
	vector<uint> _notused_EvalAlnColCountVec, _notused_EvalAlnOpenVec, _notused_EvalAlnExtVec;

	BACKGROUND_STYLE BS = BS_UniqueChains;
	if (optset_background_style)
		FeatureTrainer2::StrToBS(opt(background_style));

	FILE *fSteps = 0;
	vector<vector<float> > ScoreMx;
	FeatureTrainer2::TrainDSSFeature(FEATURE_Mu, Chains, LabelToChainIdx,
		TrainRows, TrainLabels, TrainChainIdxs,
		_notused_EvalTPs, _notused_EvalRows, _notused_EvalLabels, _notused_EvalRowChainIdxs,
		_notused_EvalAlnColCountVec, _notused_EvalAlnOpenVec, _notused_EvalAlnExtVec,
		ScoreMx, BS, 0);

	const uint AS = SIZE(ScoreMx);
	IntScoreMx.clear();
	IntScoreMx.resize(AS);

	asserta(AS == 36);
	for (uint i = 0; i < AS; ++i)
		{
		IntScoreMx[i].resize(AS);
		for (uint j = 0; j < AS; ++j)
			IntScoreMx[i][j] = int(round(ScaleFactor*ScoreMx[i][j]));
		}
	}

static void Optimize(
	const vector<string> &SpecLines,
	double &Best_y,
	vector<string> &Best_xv)
	{
	uint HJCount = 1;

	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);
	uint LatinBinCount = Peaker::SpecGetInt(GlobalSpec, "latin", UINT_MAX);
	asserta(LatinBinCount != UINT_MAX);

	Peaker &P = *new Peaker(0, "mumx");
	P.Init(SpecLines, EvalSum3);
	s_Peaker = &P;

	ProgressLog("==============\n");
	ProgressLog("Latin (%u)\n", LatinBinCount);
	ProgressLog("==============\n");
	asserta(LatinBinCount > 0);
	P.RunLatin(LatinBinCount);

	vector<uint> TopEvalIdxs;
	P.GetTopEvalIdxs(HJCount, TopEvalIdxs);
	const uint n = SIZE(TopEvalIdxs);
	if (n == 0)
		Die("No evals");

	uint EvalIdx = TopEvalIdxs[0];
	Peaker *Child = P.MakeChild("HJ");
	double y = P.m_ys[EvalIdx];
	const vector<string> &xv = P.m_xvs[EvalIdx];
	Child->AppendResult(xv, y, "HJstart");
	Child->HJ_RunHookeJeeves();
	P.AppendChildResults(*Child);
	delete Child;
	ProgressLog("============\n");
	ProgressLog("HJ converged\n");
	ProgressLog("============\n");

	Best_y = P.m_Best_y;
	Best_xv = P.m_Best_xv;
	}

// Train matrix and optimize gap penalies on SCOP40
// to find best subst mx for Mu
void cmd_hjmumx()
	{
	DSSParams::Init(DM_DefaultSensitive);

	double ScaleFactor = 1.0;
	if (optset_scale)
		ScaleFactor = opt(scale);

	vector<vector<int> > ScoreMx;
	if (optset_mxname)
		{
		ProgressLog("SetSubstMx(%s)\n", opt(mxname));
		Paralign::SetSubstMx(opt(mxname));
		}
	else
		{
		ProgressLog("Training, scale = %.1f\n", ScaleFactor);
		TrainMx(ScoreMx, ScaleFactor);
		ProgressLog("Training complete\n");
		ProgressLog("%s\n", FeatureTrainer2::m_FevStr.c_str());
		Paralign::SetMatrix(ScoreMx, 0, 0, 9999);
		Paralign::LogMatrix();
		}

	const string SpecFN = g_Arg1;
	Log("SpecFN=%s\n", SpecFN.c_str());
	vector<string> SpecLines;
	ReadLinesFromFile(SpecFN, SpecLines);

	string SeqsMethod = "muletters";
	if (optset_seqsmethod)
		SeqsMethod = opt(seqsmethod);
	s_PS = new ParaSearch;
	s_PS->GetByteSeqs(opt(input2), SeqsMethod);
	s_PS->m_SB.ReadLookup(opt(lookup));

	double Best_y;
	vector<string> Best_xv;
	Optimize(SpecLines, Best_y, Best_xv);

	ProgressLog("\n");
	ProgressLog("Best Sum3=%.3f\n", Best_y);
	ProgressLog("\n");
	}
