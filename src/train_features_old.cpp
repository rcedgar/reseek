#if 0 // @@DELETE
#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"
#include "scop40bench.h"
#include "statsig.h"

static uint GetUngappedLength(const string &Row)
	{
	const uint ColCount = SIZE(Row);
	uint L = 0;
	for (uint i = 0; i < ColCount; ++i)
		if (!isgap(Row[i]))
			++L;
	return L;
	}

static void TruncLabel(string &Label)
	{
	size_t n = Label.find(' ');
	if (n == string::npos || n == 0)
		return;
	Label.resize(n);
	}

static void TrainFeature(const string &FeatureName, const SeqDB &Input,
  const vector<PDBChain *> &Chains, LogOdds &LO)
	{
	const char *Name = FeatureName.c_str();
	const FEATURE Feature = StrToFeature(Name);
	const uint FeatureIndex = uint(Feature);

	uint AlphaSize = DSS::GetAlphaSize(Feature);
	ProgressLog("Feature %s (%u)\n", Name, AlphaSize);

	LO.Init(AlphaSize);

	const uint ChainCount = SIZE(Chains);
	map<string, uint> DomToChainIndex;
	DSS DQ;
	DSS DR;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Label = Chain.m_Label;
		DomToChainIndex[Label] = ChainIndex;

		const uint QL = Chain.GetSeqLength();
		DQ.Init(Chain);
		for (uint QPos = 0; QPos < QL; ++QPos)
			{
			uint Letter = DQ.GetFeature(Feature, QPos);
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
		string QLabel = Input.GetLabel(2*PairIndex);
		string RLabel = Input.GetLabel(2*PairIndex+1);
		TruncLabel(QLabel);
		TruncLabel(RLabel);

		map<string, uint>::const_iterator iterq = DomToChainIndex.find(QLabel);
		map<string, uint>::const_iterator iterr = DomToChainIndex.find(RLabel);
		if (iterq == DomToChainIndex.end())
			Die("Not found >%s", QLabel.c_str());
		asserta(iterr != DomToChainIndex.end());

		uint QChainIndex = iterq->second;
		uint RChainIndex = iterr->second;

		const PDBChain &QChain = *Chains[QChainIndex];
		const PDBChain &RChain = *Chains[RChainIndex];

		uint QL = QChain.GetSeqLength();
		uint RL = RChain.GetSeqLength();

		DQ.Init(QChain);
		DR.Init(RChain);

		const string &QRow = Input.GetSeq(2*PairIndex);
		const string &RRow = Input.GetSeq(2*PairIndex+1);

		uint QL2 = GetUngappedLength(QRow);
		uint RL2 = GetUngappedLength(RRow);
		asserta(QL2 == QL);
		asserta(RL2 == RL);

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
				uint ValueQ = DQ.GetFeature(FeatureIndex, QPos);
				uint ValueR = DR.GetFeature(FeatureIndex, RPos);
				LO.AddTruePair(ValueQ, ValueR);
				}
			if (!isgap(q))
				++QPos;
			if (!isgap(r))
				++RPos;
			}
		}
	}

static void WriteLOInt8(FILE *f, const string &strName, const LogOdds &LO, int8_t MaxAbsi8)
	{
	if (f == 0)
		return;
	const char *Name = strName.c_str();

	vector<double> Freqs;
	vector<vector<double> > PairFreqs;
	vector<vector<double> > LogOddsScores;
	const uint AS = LO.m_AlphaSize;
	LO.GetBackgroundFreqs(Freqs);
	LO.GetTrueFreqMx(PairFreqs);
	double ExpectedScore = LO.GetLogOddsMx(LogOddsScores);

	vector<vector<int8_t> > LogOddsScoresi;
	LO.GetLogOddsMxInt8(LogOddsScores, LogOddsScoresi, MaxAbsi8);

	fprintf(f, "FEATURE\t%s\t%u\t%.3f\n", Name, AS, ExpectedScore);
	asserta(SIZE(LogOddsScoresi) == AS);
	for (uint i = 0; i < AS; ++i)
		{
		asserta(SIZE(LogOddsScores[i]) == AS);
		fprintf(f, "S_ij\t%u", i);
		for (uint j = 0; j < AS; ++j)
			fprintf(f, "\t%d", LogOddsScoresi[i][j]);
		fprintf(f, "\n");
		}
	}

static void WriteLO(FILE *f, const string &strName, const LogOdds &LO)
	{
	if (f == 0)
		return;
	const char *Name = strName.c_str();

	vector<double> Freqs;
	vector<vector<double> > PairFreqs;
	vector<vector<double> > LogOddsScores;
	const uint AS = LO.m_AlphaSize;
	LO.GetBackgroundFreqs(Freqs);
	LO.GetTrueFreqMx(PairFreqs);
	double ExpectedScore = LO.GetLogOddsMx(LogOddsScores);

	fprintf(f, "FEATURE\t%s\t%u\t%.3f\n", Name, AS, ExpectedScore);
	asserta(SIZE(Freqs) == AS);
	for (uint i = 0; i < AS; ++i)
		{
		fprintf(f, "f_i\t%u", i);
		fprintf(f, "\t%.4g\n", Freqs[i]);
		}

	asserta(SIZE(PairFreqs) == AS);
	for (uint i = 0; i < AS; ++i)
		{
		asserta(SIZE(PairFreqs[i]) == AS);
		fprintf(f, "f_ij\t%u", i);
		for (uint j = 0; j < AS; ++j)
			fprintf(f, "\t%.4g", PairFreqs[i][j]);
		fprintf(f, "\n");
		}

	asserta(SIZE(LogOddsScores) == AS);
	for (uint i = 0; i < AS; ++i)
		{
		asserta(SIZE(LogOddsScores[i]) == AS);
		fprintf(f, "S_ij\t%u", i);
		for (uint j = 0; j < AS; ++j)
			fprintf(f, "\t%.4g", LogOddsScores[i][j]);
		fprintf(f, "\n");
		}
	}

void cmd_train_features()
	{
	optset_fast = true;
	opt(fast) = true;
	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast);
	int8_t MaxAbsi8 = 20;
	if (optset_maxi8)
		{
		uint Max = opt(maxi8);
		MaxAbsi8 = uint8_t(Max);
		asserta(uint(MaxAbsi8) == Max);
		}

	vector<string> FeatureNames;
	if (!optset_features)
		{
		for (uint i = 0; i < Params.GetFeatureCount(); ++i)
			{
			FEATURE F = Params.m_Features[i];
			string Name = FeatureToStr(F);
			FeatureNames.push_back(Name);
			}
		}
	else
		Split(opt(features), FeatureNames, '_');

	vector<PDBChain *> Chains;
	ReadChains(opt(train_cal), Chains);

	SeqDB Input;
	Input.FromFasta(g_Arg1, true);

	const uint N = SIZE(FeatureNames);
	FILE *fOut = CreateStdioFile(opt(output));
	FILE *fOut2 = CreateStdioFile(opt(output2));
	LogOdds LO;
	for (uint i = 0; i < N; ++i)
		{
		TrainFeature(FeatureNames[i], Input, Chains, LO);
		WriteLO(fOut, FeatureNames[i], LO);

		WriteLOInt8(fOut2, FeatureNames[i], LO, MaxAbsi8);
		}
	CloseStdioFile(fOut);
	CloseStdioFile(fOut2);
	}
#endif // 0
