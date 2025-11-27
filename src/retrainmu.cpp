#include "myutils.h"
#include "featuretrainer2.h"
#include "sort.h"
#include "pearson_correl.h"

const float SCALE_FACTOR = 3.5f;

// Trying to reproduce this matrix (for parasail)
extern int8_t Mu_S_ij_i8[36][36];

// Float matrix for ... not sure what
extern float musubstmx[36][36];

static void IntMxToSrc(FILE *f, const string &Name, 
  const vector<vector<int> > &Mx)
	{
	if (f == 0)
		return;
	uint N = SIZE(Mx);

	fprintf(f, "\n");
	fprintf(f, "static int %s[%u][%u] = {\n", Name.c_str(), N, N);
	fprintf(f, "//");
	for (uint i = 0; i < N; ++i)
		fprintf(f, " %4d ", i);
	fprintf(f, "\n");
	for (uint i = 0; i < N; ++i)
		{
		asserta(SIZE(Mx[i]) == N);
		fprintf(f, " {");
		for (uint j = 0; j < N; ++j)
			{
			fprintf(f, " %4d", Mx[i][j]);
			if (j+1 != N)
				fprintf(f, ",");
			}
		fprintf(f, "},\n");
		}
	fprintf(f, "};\n");
	}

static double Pear_Corr(const vector<vector<float> > &Mx)
	{
	vector<float> v1;
	vector<float> v2;
	for (uint i = 0; i < 36; ++i)
		for (uint j = 0; j < 36; ++j)
			v1.push_back((float) Mu_S_ij_i8[i][j]);

	for (uint i = 0; i < 36; ++i)
		for (uint j = 0; j < 36; ++j)
			v2.push_back((float) Mx[i][j]);
	double R = pearson_correlation(v1, v2);
	return R;
	}

static void CmpMx(const string &Msg,
	const vector<vector<float> > &ScoreMx)
	{
	double R = Pear_Corr(ScoreMx);
	ProgressLog("%s R=%.4f\n", Msg.c_str(), R);
	}

static void Cmp_ScoreMx_Mu()
	{
	vector<vector<float> > ScoreMx(36);
	for (uint i = 0; i < 36; ++i)
		{
		ScoreMx[i].resize(36);
		for (uint j = 0; j < 36; ++j)
			ScoreMx[i][j] = musubstmx[i][j];
		}
	CmpMx("extern float musubstmx", ScoreMx);
	}

void cmd_retrainmu()
	{
	Cmp_ScoreMx_Mu();

	DSSParams::Init(DM_DefaultSensitive);

	const string &ChainFN = opt(db);			// "src/reseek/test_data/scop40.bca";
	const string &TrainTPAlnFN = opt(traintps); // "src/2025-10_reseek_tune/big_fa2/tp.mints05.maxts25.fa2";
	const string &EvalTPAlnFN = opt(evaltps);	// "src/2025-10_reseek_tune/big_fa2/tp.evalrange.fa2";
	const string &EvalFPAlnFN = opt(evalfps);	// "src/2025-10_reseek_tune/big_fa2/fp.evalrange.fa2";

	FILE *fOut = CreateStdioFile(opt(output));

	FeatureTrainer2::m_FevStr.clear();

///////////////////////////////////////////////////////////////////////////////////////
// Structure chains, must include all Train and Eval chains, may include others
///////////////////////////////////////////////////////////////////////////////////////
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
	if (!opt(evalmu))
		FeatureTrainer2::AppendAlns("traintps", TrainTPAlnFN, LabelToChainIdx, true,
		  TrainRows, TrainLabels, TrainChainIdxs, TrainsTPs_notused);

	vector<bool> EvalTPs;
	vector<string> EvalRows;
	vector<string> EvalLabels;
	vector<uint> EvalRowChainIdxs;
	vector<uint> EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec;

	vector<vector<float> > ScoreMx;
	float BestArea = FLT_MAX;
	asserta(optset_evalmu || optset_background_style);
	BACKGROUND_STYLE BS = FeatureTrainer2::StrToBS(opt(background_style));
	FILE *fSteps = 0;
	FeatureTrainer2::TrainDSSFeature(FEATURE_Mu, Chains, LabelToChainIdx,
		TrainRows, TrainLabels, TrainChainIdxs,
		EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
		EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
		ScoreMx, BS, 0);

	FeatureTrainer2::ScoreMxToTsv(fOut, ScoreMx);
	CmpMx("TrainedMx", ScoreMx);

	vector<vector<int> > IntMx(36);
	for (uint i = 0; i < 36; ++i)
		{
		IntMx[i].resize(36);
		for (uint j = 0; j < 36; ++j)
			IntMx[i][j] = (int) round(ScoreMx[i][j]*SCALE_FACTOR);
		}
	IntMxToSrc(g_fLog, TrainTPAlnFN, IntMx);

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
