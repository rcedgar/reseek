#include "myutils.h"
#include "featuretrainer2.h"
#include "sfasta.h"
#include "alpha.h"
#include "round3sigfig.h"
#include "quarts.h"
#include "sort.h"

/***
		 Feature  AS   Type               Int forced        Float forced    Float not forced
		NormDens  16  Float          22412 (  1.15%)     22412 (  1.15%)     22412 (  1.15%)
		 NENDist  16  Float            116 (  0.01%)       116 (  0.01%)         0 (  0.00%)
	   HelixDens  16  Float          22412 (  1.15%)     22412 (  1.15%)     22412 (  1.15%)
	  StrandDens  16  Float          22412 (  1.15%)     22412 (  1.15%)     22412 (  1.15%)
	   DstNxtHlx  16  Float         611684 ( 31.40%)    611684 ( 31.40%)         0 (  0.00%)
	   DstPrvHlx  16  Float         643390 ( 33.02%)    643390 ( 33.02%)         0 (  0.00%)
			  NX  16  Float          22412 (  1.15%)     22412 (  1.15%)     22412 (  1.15%)
		 RENDist  16  Float         291220 ( 14.95%)    291220 ( 14.95%)         0 (  0.00%)
		  PMDist  16  Float              5 (  0.00%)         5 (  0.00%)         0 (  0.00%)
***/

// If undef has its own bin:
//  DiscardUndefs=false
//  ReplaceUndefWithThisValue=FLT_MAX
void FeatureTrainer2::GetFloatValuesAndSeqs(
	const vector<PDBChain *> &Chains,
	bool DiscardUndefsFromValuesButNotSeqs,
	float ReplaceUndefWithThisValue,
	vector<float> &SortedValues,
	vector<vector<float> > &Seqs,
	uint &UndefCount)
	{
	SortedValues.clear();
	Seqs.clear();
	UndefCount = 0;
	DSS D;
	const uint ChainCount = SIZE(Chains);
	SortedValues.reserve(ChainCount*300);
	Seqs.resize(ChainCount);
	opt_force_undef = true;
	optset_force_undef = true;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *Chains[ChainIdx];
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		vector<float> &Seq = Seqs[ChainIdx];
		Seq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = D.GetFloatFeature(m_F, Pos);
			if (Value == FLT_MAX)
				{
				++UndefCount;
				if (DiscardUndefsFromValuesButNotSeqs)
					{
					Seq.push_back(FLT_MAX);
					continue;
					}
				Value = ReplaceUndefWithThisValue;
				}
			SortedValues.push_back(Value);
			Seq.push_back(Value);
			}
		}
	opt_force_undef = false;
	optset_force_undef = false;
	sort(SortedValues.begin(), SortedValues.end());
	}

// Letters will include UINT_MAX iff
//  ReplaceUndefWithThisLetter == UINT_MAX
void FeatureTrainer2::ValuesToLetters(
	const vector<float> &Values,
	const vector<float> &BinTs,
	float UndefReplaceValue,
	vector<uint> &Letters)
	{
	const uint N = SIZE(Values);
	asserta(SIZE(BinTs) + 1 == m_AlphaSize);
	Letters.clear();
	Letters.reserve(N);
	for (uint i = 0; i < N; ++i)
		{
		float Value = Values[i];
		if (Value == FLT_MAX)
			Value = UndefReplaceValue;
		uint Letter = ValueToInt_UndefNotSpecialCase(Value, BinTs);
		Letters.push_back(Letter);
		}
	}

void FeatureTrainer2::LogLetterCountsFreqsAndBinTs(
	const vector<uint> &Counts,
	uint UndefCount,
	const vector<float> &BinTs)
	{
	vector<float> Freqs;
	FeatureTrainer2::GetFreqs(Counts, Freqs);

	Log("\n");
	Log("%4.4s", "Int");
	Log("  %10.10s", "Letters");
	Log("  %6.6s", "Freq");
	Log("  %8.8s", "BinT");
	Log("\n");
	uint Total = 0;
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		uint n = Counts[Letter];
		Log("[%2u]", Letter);
		Log("  %10u", n);
		Log("  %6.4f", Freqs[Letter]);
		if (Letter + 1 < m_AlphaSize)
			Log("  %8.3g", BinTs[Letter]);
		Log("\n");
		Total += n;
		}
	Log("Undef %u / %u (%.2f%%)\n",
		UndefCount, Total, GetPct(UndefCount, Total));
	}

void FeatureTrainer2::FloatSeqsToInt(
	const vector<vector<float> > &FloatSeqs,
	const vector<float> &BinTs,
	float UndefReplaceValue,
	vector<vector<uint> > &IntSeqs)
	{
	const uint SeqCount = SIZE(FloatSeqs);
	IntSeqs.clear();
	IntSeqs.resize(SeqCount);
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		const vector<float> &FloatSeq = FloatSeqs[SeqIdx];
		vector<uint> &IntSeq = IntSeqs[SeqIdx];
		const uint L = SIZE(FloatSeq);
		IntSeq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = FloatSeq[Pos];
			if (Value == FLT_MAX)
				Value = UndefReplaceValue;
			uint Letter = ValueToInt_UndefNotSpecialCase(Value, BinTs);
			IntSeq.push_back(Letter);
			}
		}
	}

void FeatureTrainer2::TrainFloatFeature(
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
	float UndefReplaceValue,
	float &BestArea,
	FILE *fOut)
	{
	m_FevStr += "feature=" + string(FeatureToStr(F)) + ";";
	m_FevStr += "type=float;";
	Psa(m_FevStr, "AS=%u;", AlphaSize);
	if (UndefReplaceValue == FLT_MAX)
		m_FevStr += "undef_value=*;";
	else
		Psa(m_FevStr, "undef_value=%.3g;", UndefReplaceValue);
	m_FevStr += "QS=" + string(QSToStr(QS)) + ";";

	SetFloatFeature(F, AlphaSize);
	m_QS = QS;
	m_BS = BS_Float;

	vector<float> SortedValues;
	vector<vector<float> > FloatSeqs;
	uint UndefCount1 = UINT_MAX;
	GetFloatValuesAndSeqs(Chains, false, FLT_MAX,
		SortedValues, FloatSeqs, UndefCount1);

	vector<float> BinTs;
	Quantize(SortedValues, BinTs, UndefReplaceValue);

	vector<uint> Letters;
	ValuesToLetters(SortedValues, BinTs, UndefReplaceValue, Letters);

	vector<uint> LetterCounts;
	uint UndefCount2 = UINT_MAX;
	GetLetterCounts(Letters, LetterCounts, UndefCount2);

	LogLetterCountsFreqsAndBinTs(LetterCounts, UndefCount1, BinTs);

	vector<vector<uint> > ChainIntSeqsNoUndefs;
	FloatSeqsToInt(FloatSeqs, BinTs, UndefReplaceValue,
		ChainIntSeqsNoUndefs);

	vector<vector<uint> > TrainAlnLetterPairCountMx;
	GetAlignedLetterPairCounts(ChainIntSeqsNoUndefs, TrainRows,
		TrainChainIdxs, TrainAlnLetterPairCountMx);

	TrainLogOddsMx(LetterCounts, TrainAlnLetterPairCountMx, ScoreMx);

	float BestOpenPenalty, BestExtPenalty, BestBias;
	EvalLogOddsMx(ChainIntSeqsNoUndefs, EvalRows, EvalRowChainIdxs,
		EvalTPs, EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
		ScoreMx, BestOpenPenalty, BestExtPenalty, BestBias, BestArea);

	WriteSteps(fOut, ChainIntSeqsNoUndefs, EvalRows, EvalRowChainIdxs,
		EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec, EvalTPs,
		ScoreMx, BestOpenPenalty, BestExtPenalty, BestBias);
	}
