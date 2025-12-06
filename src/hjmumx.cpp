#include "myutils.h"
#include "featuretrainer2.h"
#include "parasearch.h"
#include <list>

static ParaSearch *s_PS;
static double s_BestSum3 = DBL_MAX;
static int s_BestOpen = -999;
static int s_BestExt = -999;

static list<int> s_PendingOpens;
static list<int> s_PendingExts;

static vector<int> s_DoneOpens;
static vector<int> s_DoneExts;
static vector<double> s_DoneSum3s;

static double EvalSum3(int IntOpen, int IntExt)
	{
	s_PS->SetGapParams(IntOpen, IntExt);
	s_PS->ClearHitsAndResults();
	s_PS->Search("para");
	s_PS->Bench();
	return s_PS->m_SB.m_Sum3;
	}

static void TrainMx(const string &ChainFN,
	vector<vector<int> > &IntScoreMx, double ScaleFactor)
	{
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

static bool ExtStalled(int Open, int Ext)
	{
	double Score3_Minus1 = DBL_MAX;
	double Score3_Minus2 = DBL_MAX;
	double Score3_Minus3 = DBL_MAX;
	const uint N = SIZE(s_DoneOpens);
	asserta(SIZE(s_DoneExts) == N);
	for (uint i = 0; i < N; ++i)
		{
		int Open2 = s_DoneOpens[i];
		int Ext2 = s_DoneExts[i];
		if (Open2 == Open && Ext2 == Ext - 1)
			Score3_Minus1 = s_DoneSum3s[i];
		if (Open2 == Open && Ext2 == Ext - 2)
			Score3_Minus2 = s_DoneSum3s[i];
		if (Open2 == Open && Ext2 == Ext - 3)
			Score3_Minus3 = s_DoneSum3s[i];
		}
	if (Score3_Minus1 == DBL_MAX ||
		Score3_Minus2 == DBL_MAX ||
		Score3_Minus3 == DBL_MAX)
		return false;
	return	feq(Score3_Minus1, Score3_Minus2) &&
			feq(Score3_Minus2, Score3_Minus3);
	}

static bool AddPendingIfOk(int Open, int Ext)
	{
	if (Open < 0 || Ext < 0)
		return false;

	bool OpenFound = (std::find(s_DoneOpens.begin(), s_DoneOpens.end(), Open) != s_DoneOpens.end());
	bool ExtFound = (std::find(s_DoneExts.begin(), s_DoneExts.end(), Ext) != s_DoneExts.end());
	if (OpenFound && ExtFound)
		return false;

	if (ExtStalled(Open, Ext))
		{
		ProgressLog("\nSTALLED %d/%d\n\n", Open, Ext);
		return false;
		}

	s_PendingOpens.push_back(Open);
	s_PendingExts.push_back(Ext);

	return true;
	}

static void Optimize(int ScaleFactor)
	{
	int FirstOpen = ScaleFactor*3;
	int FirstExt = ScaleFactor*3;

	AddPendingIfOk(FirstOpen, FirstExt);
	for (uint Iter = 0; Iter < 100; ++Iter)
		{
		uint PendingCount = SIZE(s_PendingOpens);
		asserta(SIZE(s_PendingExts) == PendingCount);
		if (PendingCount == 0)
			break;
		int Open = s_PendingOpens.front();
		int Ext = s_PendingExts.front();
		s_PendingOpens.pop_front();
		s_PendingExts.pop_front();
		double Sum3 = EvalSum3(Open, Ext);
		s_DoneOpens.push_back(Open);
		s_DoneExts.push_back(Ext);
		s_DoneSum3s.push_back(Sum3);

		bool Better = false;
		if (s_BestSum3 == DBL_MAX || Sum3 > s_BestSum3)
			{
			Better = true;
			s_BestSum3 = Sum3;
			s_BestOpen = Open;
			s_BestExt = Ext;

			AddPendingIfOk(Open+1, Ext);
			AddPendingIfOk(Open-1, Ext);
			AddPendingIfOk(Open, Ext+1);
			AddPendingIfOk(Open, Ext-1);

			ProgressLog("\n >>> [%.3f]   open/ext %d/%d\n\n", s_BestSum3, Open, Ext);
			}
		else
			ProgressLog("\n ... [%.3f]   open/ext %d/%d %.3f\n\n", s_BestSum3, Open, Ext, Sum3);

		if (!Better && feq(Sum3, s_BestSum3))
			{
			AddPendingIfOk(Open+1, Ext);
			AddPendingIfOk(Open-1, Ext);
			AddPendingIfOk(Open, Ext+1);
			AddPendingIfOk(Open, Ext-1);
			}
		}

	if (!s_PendingOpens.empty())
		Warning("Failed to converge\n");
	}

// Train matrix and optimize gap penalies on SCOP40
// to find best subst mx for Mu
// H-J-like hack for integers
void cmd_hjmumx()
	{
// -fixmubyteseq applies only when reading FASTA
	asserta(!optset_fixmubyteseq);

	DSSParams::Init(DM_DefaultSensitive);

	int ScaleFactor = 1;
	if (optset_scale)
		ScaleFactor = opt(scale);
	if (optset_parabits)
		Paralign::m_Bits = opt(parabits);

	vector<vector<int> > ScoreMx;
	if (optset_mxname)
		{
		ProgressLog("SetSubstMx(%s)\n", opt(mxname));
		Paralign::SetSubstMxByName(opt(mxname));
		}
	else
		{
		ProgressLog("Training, scale = %d\n", ScaleFactor);
		const string &ChainFN = g_Arg1;				// "src/reseek/test_data/scop40.bca";
		TrainMx(g_Arg1, ScoreMx, double(ScaleFactor));
		ProgressLog("Training complete\n");
		ProgressLog("%s\n", FeatureTrainer2::m_FevStr.c_str());
		Paralign::SetMatrix(ScoreMx, 0, 0, 9999);
		Paralign::LogMatrix();
		}

	string SeqsMethod = "muletters";
	if (optset_seqsmethod)
		SeqsMethod = opt(seqsmethod);
	s_PS = new ParaSearch;
	s_PS->GetByteSeqs(opt(input2), SeqsMethod);
	s_PS->m_SB.ReadLookup(opt(lookup));

	Optimize(ScaleFactor);

	ProgressLog("\n");
	ProgressLog("Best Sum3=%.3f open=%d ext=%d scale=%d bits=%d\n",
		s_BestSum3, s_BestOpen, s_BestExt, ScaleFactor, Paralign::m_Bits);
	ProgressLog("\n");
	}
