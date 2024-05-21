#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"
#include "scop40bench.h"
// #include "outputfiles.h"

static double EvalFeature(const SeqDB &Input,
  const vector<PDBChain *> &Chains, const string &Name,
  vector<double> &Freqs, vector<vector<double> > &ScoreMx)
	{
	const char *FeatureName = Name.c_str();
	const FEATURE Feature = StrToFeature(Name.c_str());
	const uint FeatureIndex = uint(Feature);

	uint AlphaSize = DSS::GetAlphaSize(Feature);
	ProgressLog("Feature %s (%u)\n", FeatureName, AlphaSize);

	LogOdds LO;
	LO.Init(AlphaSize);

	double MinTM = 0.6;
	double MaxTM = 0.8;
	if (optset_mintm)
		MinTM = opt_mintm;
	if (optset_maxtm)
		MaxTM = opt_maxtm;

	const uint ChainCount = SIZE(Chains);
	map<string, uint> DomToChainIndex;
	DSS QX;
	DSS RX;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Label = Chain.m_Label;
		vector<string> Fields;
		string Dom;
		GetScopDomFromLabel(Label, Dom);
		DomToChainIndex[Dom] = ChainIndex;

		const uint QL = Chain.GetSeqLength();
		QX.Init(Chain);
		for (uint QPos = 0; QPos < QL; ++QPos)
			{
			uint Letter = QX.GetFeature(Feature, QPos);
			LO.AddBackgroundLetter(Letter);
			}
		}

	const uint SeqCount = Input.GetSeqCount();
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	uint LetterPairCount = 0;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Processing");
		const string &QLabel = Input.GetLabel(2*PairIndex);
		const string &RLabel = Input.GetLabel(2*PairIndex+1);
		vector<string> Fields;
		Split(QLabel, Fields, '/');
		asserta(SIZE(Fields) == 4);
		const string &QDom = Fields[0];
		const string &Fam = Fields[1];
		const string sTM = Fields[2];
		const string sPctId = Fields[3];
		double TM = StrToFloat(sTM);
		if (TM < MinTM || TM > MaxTM)
			continue;

		string RDom;
		GetScopDomFromLabel(RLabel, RDom);
		uint QChainIndex = DomToChainIndex[QDom];
		uint RChainIndex = DomToChainIndex[RDom];
		const PDBChain &QChain = *Chains[QChainIndex];
		const PDBChain &RChain = *Chains[RChainIndex];
		uint QL = QChain.GetSeqLength();
		uint RL = RChain.GetSeqLength();
		QX.Init(QChain);
		RX.Init(RChain);
		const string &QRow = Input.GetSeq(2*PairIndex);
		const string &RRow = Input.GetSeq(2*PairIndex+1);
		const uint ColCount = SIZE(QRow);
		asserta(SIZE(QRow) == SIZE(RRow));

		uint QPos = 0;
		uint RPos = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char q = QRow[Col];
			char r = RRow[Col];
			if (!isgap(q) && !isgap(r))
				{
				uint ValueQ = QX.GetFeature(FeatureIndex, QPos);
				uint ValueR = RX.GetFeature(FeatureIndex, RPos);
				LO.AddTruePair(ValueQ, ValueR);
				}
			if (!isgap(q))
				++QPos;
			if (!isgap(r))
				++RPos;
			}
		}
	LO.GetBackgroundFreqs(Freqs);
	double ExpectedScore = LO.GetLogOddsMx(ScoreMx);
	return ExpectedScore;
	}

void cmd_eval_feature()
	{
	vector<PDBChain *> Chains;
	ReadChains(opt_train_cal, Chains);

	vector<string> FeatureNames;
	Split(opt_features, FeatureNames, '_');

	SeqDB Input;
	Input.FromFasta(g_Arg1, true);

	const uint N = SIZE(FeatureNames);
	vector<double> ExpectedScores(N);
	vector<uint> AlphaSizes(N);
	vector<vector<double> > FreqsVec(N);
	vector<vector<vector<double> > > ScoreMxVec(N);
	for (uint i = 0; i < N; ++i)
		{
		ExpectedScores[i] = 
		  EvalFeature(Input, Chains, FeatureNames[i], FreqsVec[i], ScoreMxVec[i]);
		AlphaSizes[i] = SIZE(FreqsVec[i]);
		}

	FILE *f = CreateStdioFile(opt_output);
	asserta(f != 0);
	fprintf(f, "/********************************\n");
	fprintf(f, "Expected scores\n");
	for (uint i = 0; i < N; ++i)
		{
		const char *Name = FeatureNames[i].c_str();
		fprintf(f, "%6.4f  %2u  %s\n",
		  ExpectedScores[i], AlphaSizes[i], Name);
		}

	for (uint i = 0; i < N; ++i)
		{
		const char *Name = FeatureNames[i].c_str();
		const vector<double> &Freqs = FreqsVec[i];
		const uint AS = AlphaSizes[i];
		fprintf(f, "\n");
		fprintf(f, "Freqs %s/%u\n", Name, AS);
		double Sum = 0;
		for (uint i = 0; i < AS; ++i)
			{
			double Freq = Freqs[i];
			fprintf(f, "   [%2u]  %.6f\n", i, Freq);
			Sum += Freq;
			}
		asserta(feq(Sum, 1));
		}
	fprintf(f, "********************************/\n");

	for (uint Idx = 0; Idx < N; ++Idx)
		{
		const char *Name = FeatureNames[Idx].c_str();
		const vector<vector<double> > &ScoreMx = ScoreMxVec[Idx];
		const uint AS = AlphaSizes[Idx];

		fprintf(f, "\n");
		fprintf(f, "static double ScoreMx_%s[%u][%u] = {\n", Name, AS, AS);
		for (uint i = 0; i < AS; ++i)
			{
			fprintf(f, "  {");
			for (uint j = 0; j < AS; ++j)
				fprintf(f, " %7.4f,", ScoreMx[i][j]);
			fprintf(f, "  }, // %u\n", i);
			}
		fprintf(f, "};\n");
		}

	fprintf(f, "\n");
	fprintf(f, "float **g_ScoreMxs2[FEATURE_COUNT];\n");
	fprintf(f, "uint g_AlphaSizes2[FEATURE_COUNT];\n");
	fprintf(f, "static bool Init()\n");
	fprintf(f, "	{\n");
	for (uint Idx = 0; Idx < N; ++Idx)
		{
		const char *Name = FeatureNames[Idx].c_str();
		const vector<vector<double> > &ScoreMx = ScoreMxVec[Idx];
		const uint AS = AlphaSizes[Idx];
		fprintf(f, "	asserta(DSS::GetAlphaSize(FEATURE_%s) == %u);\n", Name, AS);
		fprintf(f, "	g_AlphaSizes2[FEATURE_%s] = %u;\n", Name, AS);
		}
	for (uint Idx = 0; Idx < N; ++Idx)
		{
		const char *Name = FeatureNames[Idx].c_str();
		const vector<vector<double> > &ScoreMx = ScoreMxVec[Idx];
		const uint AS = AlphaSizes[Idx];
		fprintf(f, "\n");
		fprintf(f, "	g_ScoreMxs2[FEATURE_%s] = myalloc(float *, %u);\n", Name, AS);
		fprintf(f, "	for (uint i = 0; i < %u; ++i)\n", AS);
		fprintf(f, "		{\n");
		fprintf(f, "		g_ScoreMxs2[FEATURE_%s][i] = myalloc(float, %u);\n", Name, AS);
		fprintf(f, "		for (uint j = 0; j < %u; ++j)\n", AS);
		fprintf(f, "			g_ScoreMxs2[FEATURE_%s][i][j] = (float) ScoreMx_%s[i][j];\n", Name, Name);
		fprintf(f, "		}\n");
		}
	fprintf(f, "	return true;\n");
	fprintf(f, "	}\n");
	fprintf(f, "static bool InitDone = Init();\n");

	CloseStdioFile(f);
	}
