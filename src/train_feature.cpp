#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"
#include "scop40bench.h"

void GetFloatFeatureValues(const vector<PDBChain *> &Chains, 
						   FEATURE F, vector<float> &Values);
void ConstructBinsFromValues(const vector<float> &Values, uint AS,
							 vector<float> &BinTs);

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

void TrainFeature(const SeqDB &Input, const vector<PDBChain *> &Chains,
				  FEATURE F, uint AS, 
				  vector<float> &BinTs, LogOdds &LO)
	{
	BinTs.clear();
	LO.Init(AS);
	bool IsInt = FeatureIsInt(F);
	if (!IsInt)
		{
		vector<float> Values;
		GetFloatFeatureValues(Chains, F, Values);
		ConstructBinsFromValues(Values, AS, BinTs);
		asserta(SIZE(BinTs) + 1 == AS);
		}

	const uint ChainCount = SIZE(Chains);
	map<string, uint> DomToChainIndex;
	DSS DQ;
	DSS DR;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Background frequencies");

		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Label = Chain.m_Label;
		DomToChainIndex[Label] = ChainIndex;

		const uint QL = Chain.GetSeqLength();
		DQ.Init(Chain);
		for (uint QPos = 0; QPos < QL; ++QPos)
			{
			uint Letter = UINT_MAX;
			if (IsInt)
				Letter = DQ.GetFeature(F, QPos);
			else
				{
				float Value = DQ.GetFloatFeature(F, QPos);
				Letter = DSS::ValueToInt(BinTs, Value);
				}
			LO.AddBackgroundLetter(Letter);
			}
		}

	const uint SeqCount = Input.GetSeqCount();
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	uint LetterPairCount = 0;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Joint frequencies");
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
				uint LetterQ = UINT_MAX;
				uint LetterR = UINT_MAX;
				if (IsInt)
					{
					LetterQ = DQ.GetFeature(F, QPos);
					LetterR = DR.GetFeature(F, RPos);
					}
				else
					{
					float ValueQ = DQ.GetFloatFeature(F, QPos);
					float ValueR = DR.GetFloatFeature(F, RPos);

					LetterQ = DSS::ValueToInt(BinTs, ValueQ);
					LetterR = DSS::ValueToInt(BinTs, ValueR);
					}
				LO.AddTruePair(LetterQ, LetterR);
				}
			if (!isgap(q))
				++QPos;
			if (!isgap(r))
				++RPos;
			}
		}
	}

void WriteLOInt8(FILE *f, const string &strName, const LogOdds &LO, int8_t MaxAbsi8)
	{
	if (f == 0)
		return;
	const char *Name = strName.c_str();

	vector<float> Freqs;
	vector<vector<float> > PairFreqs;
	vector<vector<float> > LogOddsScores;
	const uint AS = LO.m_AlphaSize;
	LO.GetBackgroundFreqs(Freqs);
	LO.GetTrueFreqMx(PairFreqs);
	float ExpectedScore = LO.GetLogOddsMx(LogOddsScores);

	vector<vector<int8_t> > LogOddsScoresi;
	LO.GetLogOddsMxInt8(LogOddsScores, LogOddsScoresi, MaxAbsi8);

	fprintf(f, "FEATURE\t%s\t%u\t%.3f\n", Name, AS, ExpectedScore);
	asserta(SIZE(LogOddsScoresi) == AS);
	for (uint i = 0; i < AS; ++i)
		{
		asserta(SIZE(LogOddsScores[i]) == AS);
		fprintf(f, "S_ij_i8\t%u", i);
		for (uint j = 0; j < AS; ++j)
			fprintf(f, "\t%d", LogOddsScoresi[i][j]);
		fprintf(f, "\n");
		}
	}

void WriteLO(FILE *f, const string &strName, const LogOdds &LO)
	{
	if (f == 0)
		return;
	const char *Name = strName.c_str();

	vector<float> Freqs;
	vector<vector<float> > PairFreqs;
	vector<vector<float> > LogOddsScores;
	const uint AS = LO.m_AlphaSize;
	LO.GetBackgroundFreqs(Freqs);
	LO.GetTrueFreqMx(PairFreqs);
	float ExpectedScore = LO.GetLogOddsMx(LogOddsScores);

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

void cmd_train_feature()
	{
	asserta(optset_feature);
	const string FeatureName = string(opt(feature));
	FEATURE F = StrToFeature(FeatureName.c_str());

	vector<uint> AlphaSizes;
	bool IsInt = FeatureIsInt(F);
	if (IsInt)
		{
		uint AS = GetAlphaSize(F);
		asserta(AS != 0 && AS != UINT_MAX);
		AlphaSizes.push_back(AS);
		}
	else
		{
		if (optset_alpha_size)
			{
			uint AS = opt(alpha_size);
			asserta(AS != 0 && AS != UINT_MAX);
			AlphaSizes.push_back(AS);
			}
		else if (optset_alpha_sizes)
			{
			vector<string> Fields;
			Split(opt(alpha_sizes), Fields, ',');
			for (uint i = 0; i < SIZE(Fields); ++i)
				{
				uint AS = StrToUint(Fields[i]);
				AlphaSizes.push_back(AS);
				}
			}
		else
			Die("Must set -alpha_size[s]");
		}
	asserta(!AlphaSizes.empty());

	SeqDB Input;
	Input.FromFasta(g_Arg1, true);

	optset_fast = true;
	opt(fast) = true;
	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast, SCOP40_DBSIZE);
	int8_t MaxAbsi8 = 20;
	if (optset_maxi8)
		{
		uint Max = opt(maxi8);
		MaxAbsi8 = uint8_t(Max);
		asserta(uint(MaxAbsi8) == Max);
		}

	vector<PDBChain *> Chains;
	ReadChains(opt(train_cal), Chains);

	for (uint Idx = 0; Idx < SIZE(AlphaSizes); ++Idx)
		{
		uint AS = AlphaSizes[Idx];
		LogOdds LO;
		vector<float> BinTs;
		TrainFeature(Input, Chains, F, AS, BinTs, LO);
		WriteLO(g_fLog, FeatureName, LO);
		WriteLOInt8(g_fLog, FeatureName, LO, MaxAbsi8);
		}
	}
