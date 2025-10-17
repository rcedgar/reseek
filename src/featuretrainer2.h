#pragma once

#include "features.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "dss.h"
#include "peaker.h"

class FeatureTrainer2
	{
public:
static FEATURE m_F;
static uint m_AlphaSize;

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

static void LogChainIntSeqsStats(const string &Msg,
	const vector<vector<uint> > &Seqs);

static void GetAlignedLetterCounts(
	const vector<vector<uint> > &ChainIntSeqs,
	const vector<string> &Rows,
	const vector<uint> &RowChainIdxs,
	bool IgnoreUndef,
	vector<uint> &Counts,
	uint &LetterCount,
	uint &UndefCount);

static void GetAlignedLetterPairCounts(
	const vector<vector<uint> > &ChainIntSeqs,
	const vector<string> &Rows,
	const vector<uint> &RowChainIdxs,
	bool IgnoreUndef,
	vector<vector<uint> > &CountMx);

static void GetLetterCounts(
	const vector<vector<uint> > &ChainIntSeqs,
	const vector<uint> &ChainIdxs,
	bool IgnoreUndef,
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
	vector<vector<float> > &ScoreMx,
	const vector<float> &Freqs);

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

static void GetChainIntSeqs_Float(
	const vector<PDBChain *> &Chains,
	vector<vector<uint> > &IntSeqs,
	const vector<float> &BinTs,
	uint &LetterCount,
	uint &UndefCount);

static void GetAlnSubstScores(
	const vector<vector<uint> > &ChainIntSeqs,
	const vector<string> &Rows,
	const vector<uint> &RowChainIdxs,
	bool IgnoreUndef,
	uint ReplaceUndefByThisLetter,
	const vector<vector<float > > &ScoreMx,
	vector<float> &SubstScores);

static void Quantize(
	const vector<float> &Values,
	vector<float> &BinTs);

static void QuantizeUniques(
	const vector<float> &SortedValues,
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

static void TrainIntFeatureNoUndefs(
	FEATURE F,
	const string &ChainFN,
	const string &TrainTPAlnFN,
	const string &TrainFPAlnFN,
	const string &EvalTPAlnFN,
	const string &EvalFPAlnFN,
	vector<vector<float > > &ScoreMx);

static void TrainIntFeatureIgnoreUndefs(
	FEATURE F,
	const string &ChainFN,
	const string &TrainTPAlnFN,
	const string &TrainFPAlnFN,
	const string &EvalTPAlnFN,
	const string &EvalFPAlnFN,
	vector<vector<float > > &ScoreMx);
};
