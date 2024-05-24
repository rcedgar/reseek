#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"

void GetScopDomFromLabel(const string &Label, string &Dom);

static uint g_FeatureIndex = UINT_MAX;

static void ReportBins(
  const vector<double> &Values, 
  const vector<double> &Values1,
  const vector<double> &Values2,
  uint AlphaSize)
	{
	const char *FeatureName = FeatureToStr(g_FeatureIndex);
	const uint K = SIZE(Values);
	asserta(K > 0);
	vector<double> Ts;
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		{
		uint k = ((i+1)*K)/AlphaSize;
		double t = Values[k];
		Ts.push_back(t);
		}

	LogOdds LO;
	LO.Init(AlphaSize);
	const uint PairCount = SIZE(Values1);
	asserta(SIZE(Values2) == PairCount);
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		double Value1 = Values1[PairIndex];
		double Value2 = Values2[PairIndex];
		uint i1 = DSS::ValueToInt(Ts, Value1);
		uint i2 = DSS::ValueToInt(Ts, Value2);
		LO.AddBackgroundLetter(i1);
		LO.AddBackgroundLetter(i2);
		LO.AddTruePair(i1, i2);
		}

	vector<vector<double> > ScoreMx;
	double ExpectedScore = LO.GetLogOddsMx(ScoreMx);
	ProgressLog("%s: AlphaSize %u, ExpectedScore %.4g\n",
	  FeatureName, AlphaSize, ExpectedScore);

	Log("\n// %s [%2u] expected score %.4g\n",
	  FeatureName, AlphaSize, ExpectedScore);
	Log("ALPHA_SIZE(%s, %u);\n", FeatureName, AlphaSize);

	Log("BIN_T_BEGIN(%s);\n", FeatureName);
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		Log("BIN_T(%s, %u, %.4g);\n",
		  FeatureName, i, Ts[i]);
	Log("BIN_T_END(%s);\n", FeatureName);

	if (optset_alpha_size)
		{
		string MxName = string(FeatureName) + "_ScoreMx";
		Log("\n");
		LO.MxToSrc(g_fLog, MxName, ScoreMx);
		Log("\n");
		}
	}

void cmd_float_feature_bins()
	{
	asserta(optset_feature);
	g_FeatureIndex = StrToFeatureIndex(opt_feature);
	const char *FeatureName = FeatureToStr(g_FeatureIndex);

	SeqDB Input;
	Input.FromFasta(opt_float_feature_bins, true);

	vector<PDBChain *> Chains;
	ReadChains(opt_train_cal, Chains);
	const uint ChainCount = SIZE(Chains);
	map<string, uint> DomToChainIndex;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const string &Label = Chains[ChainIndex]->m_Label;
		string Dom;
		GetScopDomFromLabel(Label, Dom);
		DomToChainIndex[Dom] = ChainIndex;
		}

	const uint SeqCount = Input.GetSeqCount();
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	uint LetterPairCount = 0;
	DSS QX;
	DSS RX;
	vector<double> Values;
	vector<double> Values1;
	vector<double> Values2;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Processing");
		const string &QLabel = Input.GetLabel(2*PairIndex);
		const string &RLabel = Input.GetLabel(2*PairIndex+1);
		vector<string> Fields;
		Split(QLabel, Fields, '/');
		asserta(SIZE(Fields) == 4);
		const string &QDom = Fields[0];

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
				double ValueQ = QX.GetFloatFeature(g_FeatureIndex, QPos);
				double ValueR = RX.GetFloatFeature(g_FeatureIndex, RPos);
				if (ValueQ == DBL_MAX || ValueR == DBL_MAX)
					continue;
				Values.push_back(ValueQ);
				Values.push_back(ValueR);
				Values1.push_back(ValueQ);
				Values2.push_back(ValueR);
				}
			if (!isgap(q))
				++QPos;
			if (!isgap(r))
				++RPos;
			}
		}

	const uint K = SIZE(Values);
	asserta(K > 0);
	sort(Values.begin(), Values.end());
	ProgressLog("Value range %.3g .. %.3g\n",
	  Values[0], Values[K-1]);

	if (optset_alpha_size)
		ReportBins(Values, Values1, Values2, opt_alpha_size);
	else
		{
		ReportBins(Values, Values1, Values2, 3);
		ReportBins(Values, Values1, Values2, 4);
		ReportBins(Values, Values1, Values2, 6);
		ReportBins(Values, Values1, Values2, 8);
		ReportBins(Values, Values1, Values2, 10);
		ReportBins(Values, Values1, Values2, 12);
		ReportBins(Values, Values1, Values2, 16);
		ReportBins(Values, Values1, Values2, 24);
		ReportBins(Values, Values1, Values2, 32);
		}
	}
