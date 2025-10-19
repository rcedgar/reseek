#pragma once

#include "features.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "dss.h"
#include "peaker.h"

enum QUANTIZE_STYLE
	{
	QS_Invalid,
	QS_DiscardUndef,
	QS_UndefDistinct,
	QS_UndefOverlapMedian,
	QS_UndefNotSpecialCase,
	};

class FeatureTrainer2
	{
public:
static FEATURE m_F;
static uint m_AlphaSize;
static string m_BgMethod;

static void SetFloatFeature(FEATURE F, uint AlphaSize);
static void SetIntFeature(FEATURE F);

static void TruncLabel(string &Label);

static void TruncLabel(const string &Label,
	string &TruncatedLabel);

static void AppendAlns(
	const string &FN,
	const map<string, uint> &LabelToChainIdx,
	bool AlnsAreTPs,
	vector<string> &Rows,
	vector<string> &Labels,
	vector<uint> &ChainIdxs,
	vector<bool> &TPs);

static void AddAln(const string &Label1, const string &Row1,
	const string &Label2, const string &Row2, bool TP);

static void ReadChains(
	const string &ChainsFN,
	vector<PDBChain *> &Chains,
	map<string, uint> &LabelToChainIdx);

static void ScoreMxToSrc(
	FILE *f,
	const vector<vector<float> > &ScoreMx);

static void FreqsToSrc(
	FILE *f,
	const vector<float> &Freqs);

static void GetSortedFloatValues(
	const vector<PDBChain *> &Chains,
	DSS &D,
	bool IgnoreUndef,
	float ReplaceUndefValue,
	vector<float> &Values,
	uint &UndefCount,
	uint &ReplaceCount);

static void LogSortedValueStats(
	const string &Msg,
	const vector<float> &Values);

static void GetFloatSeqs(
	const vector<PDBChain *> &Chains,
	DSS &D,
	float ReplaceUndefValue,
	vector<vector<float> > &Seqs);

static void LogChainFloatSeqsStats(
	const string &Msg,
	const vector<vector<float> > &Seqs);

static void GetIntSeqs(
	const vector<PDBChain *> &Chains,
	DSS &D,
	uint ReplaceUndefValue,
	const vector<float> &BinTs, vector<vector<uint> > &Seqs,
	uint &UndefCount);

static void LogChainIntSeqsStats(
	const vector<vector<uint> > &Seqs);

static void GetAlignedLetterCounts(
	const vector<vector<uint> > &ChainIntSeqsNoUndefs,
	const vector<string> &Rows,
	const vector<uint> &RowChainIdxs,
	vector<uint> &Counts);

static void GetAlignedLetterPairCounts(
	const vector<vector<uint> > &ChainIntSeqsNoUndefs,
	const vector<string> &Rows,
	const vector<uint> &RowChainIdxs,
	vector<vector<uint> > &CountMx);

static void GetAllLetterCountsUniqueChains(
	const vector<vector<uint> > &ChainIntSeqsNoUndefs,
	const vector<uint> &ChainIdxs,
	vector<uint> &Counts);

static void GetIntSeqLetterCounts(
	const vector<vector<uint> > &Seqs,
	vector<uint> &Counts,
	uint &UndefCount);

static void GetFreqs(
	const vector<uint> &Counts,
	vector<float> &Freqs);

static void GetFreqMx(
	const vector<vector<uint> > &CountMx,
	vector<vector<float> > &FreqMx);

static void GetFreqDiffs(
	const vector<float> &Freqs1,
	const vector<float> &Freqs2,
	vector<float> &Diffs);

static void GetLogOddsMx(
	const vector<float> &Freqs,
	const vector<vector<float> > &FreqMx,
	vector<vector<float> > &ScoreMx);

static float GetExpectedScore(
	const vector<float> &Freqs,
	vector<vector<float> > &ScoreMx);

static float GetShannonEntropy(
	const vector<vector<float> > &FreqMx);

static float GetRelativeEntropy(
	const vector<vector<float> > &FreqMx,
	const vector<vector<float> > &ScoreMx);

static void GetGapCounts(
	const string &Row1,
	const string &Row2,
	uint &ColCount,
	uint &Opens,
	uint &Exts);

static void GetGapCountVecs(
	const vector<string> &Rows,
	vector<uint> &ColCountVec,
	vector<uint> &OpenVec,
	vector<uint> &ExtVec);

static void BinTsToTsv(
	FILE *f,
	const vector<float> &Ts);

static void Round3SigFig(
	const vector<float> &InputValues,
	vector<float> &RoundedValues);

static void GetChainIntSeqs_Int(
	const vector<PDBChain *> &Chains,
	vector<vector<uint> > &IntSeqs,
	uint &LetterCount,
	uint &UndefCount);

static void FloatSeqsToInt(
	const vector<vector<float> > &FloatSeqs,
	const vector<float> &BinTs,
	float UndefReplaceValue,
	vector<vector<uint> > &IntSeqs);

static void GetChainIntSeqs_DSS(
	const vector<PDBChain *> &Chains,
	vector<vector<uint> > &IntSeqs);

static void GetAlnSubstScores(
	const vector<vector<uint> > &ChainIntSeqsNoUndefs,
	const vector<string> &Rows,
	const vector<uint> &RowChainIdxs,
	bool UndefsAllowed,
	uint ReplaceUndefByThisLetter,
	const vector<vector<float > > &ScoreMx,
	vector<float> &SubstScores);

static void ReplaceUndefs(
	const vector<vector<uint> > &ChainIntSeqs,
	uint ReplacementLetter,
	vector<vector<uint> > &ChainIntSeqsNoUndefs);

static void Quantize(
	const vector<float> &Values,
	QUANTIZE_STYLE QS,
	vector<float> &BinTs,
	float &UndefReplaceValue);

static void Quantize_DiscardUndef(
	const vector<float> &Values,
	vector<float> &BinTs);

static void Quantize_UndefDistinct(
	const vector<float> &Values,
	vector<float> &BinTs);

static void Quantize_UndefOverlapMedian(
	const vector<float> &Values,
	vector<float> &BinTs,
	float &Median);

static void Quantize_UndefNotSpecialCase(
	const vector<float> &Values,
	vector<float> &BinTs);

static void GetAlnScores(
	const vector<float> &AlnSubstScores,
	const vector<uint> &AlnColCountVec,
	const vector<uint> &AlnOpenVec,
	const vector<uint> &AlnExtVec,
	float OpenPenalty,
	float ExtPenalty,
	float Bias,
	vector<float> &AlnScores);

static void LogAlnScoreQuarts(
	const vector<float> &AlnScores,
	const vector<bool> &TPs);

static void GetSteps(
	const vector<float> &AlnScores,
	const vector<bool> &TPs,
	vector<float> &StepScores,
	vector<float> &StepTPfs,
	vector<float> &StepFPfs);

static float CalcArea(
	const vector<float> &AlnScores,
	const vector<bool> &TPs);

static void LogFreqs(
	const vector<float> &Freqs);

static void OptimizeArea(
	const vector<float> &EvalAlnSubstScores,
	const vector<uint> &EvalAlnColCountVec,
	const vector<uint> &EvalAlnOpenVec,
	const vector<uint> &EvalAlnExtVec,
	const vector<bool> &EvalTPs,
	float &OpenPenalty,
	float &ExtPenalty,
	float &Bias,
	float &Area,
	uint Iters);

static void TrainLogOddsMx(
	const vector<uint> &Counts,
	const vector<vector<uint> > &CountMx,
	vector<vector<float> > &ScoreMx);

static void LoadEvalAlns(
	const string &EvalTPAlnFN,
	const string &EvalFPAlnFN,
	const map<string, uint> &LabelToChainIdx,
	vector<string> &EvalRows,
	vector<string> &EvalLabels,
	vector<uint> &EvalRowChainIdxs,
	vector<bool> &EvalTPs,
	vector<uint> &EvalAlnColCountVec,
	vector<uint> &EvalAlnOpenVec,
	vector<uint> &EvalAlnExtVec);

static void TrainIntFeature(
	FEATURE F,
	const vector<PDBChain *> &Chains,
	const map<string, uint> &LabelToChainIdx,
	const vector<string> &TrainRows,
	const vector<string> &TrainLabels,
	const vector<uint> &TrainChainIdxs,
	const vector<bool> &EvalTPs,
	const vector<string> &EvalRows,
	const vector<string> &EvalLabels,
	const vector<uint> &EvalRowChainIdxs,
	const vector<uint> &EvalAlnColCountVec,
	const vector<uint> &EvalAlnOpenVec,
	const vector<uint> &EvalAlnExtVec,
	bool UndefsAllowed,
	uint ReplaceUndefWithThisLetter,
	const string &BgMethod,
	vector<vector<float > > &ScoreMx,
	float &BestArea);

static void GetFloatValuesAndSeqs(
	const vector<PDBChain *> &Chains,
	bool DiscardUndefsFromValuesButNotSeqs,
	float ReplaceUndefWithThisValue,
	vector<float> &SortedValues,
	vector<vector<float> > &Seqs,
	uint &UndefCount);

static void ValuesToLetters(
	const vector<float> &Values,
	const vector<float> &BinTs,
	float UndefReplaceValue,
	vector<uint> &Letters);

static void GetLetterCounts(
	const vector<uint> &Letters,
	vector<uint> &Counts,
	uint &UndefCount);

static void LogLetterCountsFreqsAndBinTs(
	const vector<uint> &Counts,
	uint UndefCount,
	const vector<float> &BinTs);

// Special case bin threshold for undef, Letter=m_AlphaSize-1
// BinTs[AS-1]=FLT_MAX
// All values except FLT_MAX will be < this threshold
// *MUST* use < *NOT* <= in ValueToInt_xxx
static inline uint ValueToInt_UndefDistinctLetter(
	float Value, const vector<float> &Ts)
	{
	asserta(SIZE(Ts) + 1 == m_AlphaSize);
	asserta(Ts[m_AlphaSize-1] == FLT_MAX);
	for (uint i = 0; i + 1 < m_AlphaSize; ++i)
		if (Value < Ts[i])
			return i;
	return m_AlphaSize - 1;
	}

// Special case bin threshold for undef, Letter=m_AlphaSize-1
//		this achieves no special case in ValueToInt_xxx
// BinTs[AS-1]=FLT_MAX
// All values except FLT_MAX will be < this threshold
// *MUST* use < *NOT* <= in ValueToInt_xxx
static inline uint ValueToInt_UndefNotSpecialCase(
	float Value, const vector<float> &Ts)
	{
	asserta(SIZE(Ts) + 1 == m_AlphaSize);
	for (uint i = 0; i + 1 < m_AlphaSize; ++i)
		if (Value < Ts[i])
			return i;
	return m_AlphaSize - 1;
	}

static inline uint ValueToInt_UndefOverlapValue(
	float Value, const vector<float> &Ts, float UndefValue)
	{
	asserta(UndefValue < FLT_MAX);
	if (Value == FLT_MAX)
		Value = UndefValue;
	for (uint i = 0; i + 1 < m_AlphaSize; ++i)
		if (Value <= Ts[i])
			return i;
	return m_AlphaSize - 1;
	}

static inline uint ValueToInt_UndefOverlapLetter(
	float Value, const vector<float> &Ts, uint UndefLetter)
	{
	asserta(UndefLetter < m_AlphaSize);
	if (Value == FLT_MAX)
		return UndefLetter;
	for (uint i = 0; i + 1 < m_AlphaSize; ++i)
		if (Value <= Ts[i])
			return i;
	return m_AlphaSize - 1;
	}

static void EvalLogOddsMx(
	const vector<vector<uint> > &ChainIntSeqsNoUndefs,
	const vector<string> &EvalRows,
	const vector<uint> &EvalRowChainIdxs,
	const vector<bool> &EvalTPs,
	const vector<uint> &EvalAlnColCountVec,
	const vector<uint> &EvalAlnOpenVec,
	const vector<uint> &EvalAlnExtVec,
	const vector<vector<float> > &ScoreMx,
	float &BestArea);

static void GetDSSScoreMx(
	FEATURE F,
	vector<vector<float> > &ScoreMx);

static void TrainFloatFeature(
	FEATURE F,
	uint AlphaSize,
	const vector<PDBChain *> &Chains,
	const map<string, uint> &LabelToChainIdx,
	const vector<string> &TrainRows,
	const vector<string> &TrainLabels,
	const vector<uint> &TrainChainIdxs,
	const vector<bool> &EvalTPs,
	const vector<string> &EvalRows,
	const vector<string> &EvalLabels,
	const vector<uint> &EvalRowChainIdxs,
	const vector<uint> &EvalAlnColCountVec,
	const vector<uint> &EvalAlnOpenVec,
	const vector<uint> &EvalAlnExtVec,
	vector<vector<float > > &ScoreMx,
	QUANTIZE_STYLE QS,
	float &BestArea);

static void TrainDSSFeature(
	FEATURE F,
	const vector<PDBChain *> &Chains,
	const map<string, uint> &LabelToChainIdx,
	const vector<string> &TrainRows,
	const vector<string> &TrainLabels,
	const vector<uint> &TrainChainIdxs,
	const vector<bool> &EvalTPs,
	const vector<string> &EvalRows,
	const vector<string> &EvalLabels,
	const vector<uint> &EvalRowChainIdxs,
	const vector<uint> &EvalAlnColCountVec,
	const vector<uint> &EvalAlnOpenVec,
	const vector<uint> &EvalAlnExtVec,
	const string &BgMethod,
	vector<vector<float > > &ScoreMx,
	float &BestArea);
};
