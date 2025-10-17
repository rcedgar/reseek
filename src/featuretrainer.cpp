#include "myutils.h"
#include "featuretrainer.h"
#include "sfasta.h"
#include "alpha.h"
#include "valuetointtpl.h"
#include "round3sigfig.h"
#include "quarts.h"
#include "sort.h"

static bool GetRequiredStyleOpt(const string &Style, char c)
	{
	if (Style.find(toupper(c)) != string::npos)
		return true;
	else if (Style.find(tolower(c)) != string::npos)
		return false;
	else
		Die("-style missing %c", c);
	return false;
	}

uint FeatureTrainer::GetUngappedLength(const string &Row) const
	{
	const uint ColCount = SIZE(Row);
	uint L = 0;
	for (uint i = 0; i < ColCount; ++i)
		if (!isgap(Row[i]))
			++L;
	return L;
	}

void FeatureTrainer::TruncLabel(string &Label)
	{
	size_t n = Label.find(' ');
	if (n != string::npos)
		Label.resize(n);
	n = Label.find('|');
	if (n != string::npos)
		Label.resize(n);
	n = Label.find('/');
	if (n != string::npos)
		Label.resize(n);
	}

void FeatureTrainer::ReadChains(const string &ChainsFN)
	{
	::ReadChains(ChainsFN, m_Chains);
	SetLabelToChainIndex();
	}

void FeatureTrainer::GetGapCounts(const string &Row1, const string &Row2,
	uint &Opens, uint &Exts) const
	{
	const uint ColCount = SIZE(Row1);
	asserta(SIZE(Row2) == ColCount);
	bool InGap = false;
	Opens = 0;
	Exts = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		bool gap1 = isgap(Row1[Col]);
		bool gap2 = isgap(Row2[Col]);
		asserta(!(gap1 && gap2));
		if (gap1 || gap2)
			{
			if (InGap)
				++Exts;
			else
				{
				InGap = true;
				++Opens;
				}
			}
		else
			InGap = false;
		}
	}

void FeatureTrainer::AddAln(const string &Label1, const string &Row1,
	const string &Label2, const string &Row2, bool TP)
	{
	map<string, uint>::const_iterator iter1 = m_LabelToChainIndex.find(Label1);
	map<string, uint>::const_iterator iter2 = m_LabelToChainIndex.find(Label2);
	if (iter1 == m_LabelToChainIndex.end())
		Die("Not found >%s", Label1.c_str());
	if (iter2 == m_LabelToChainIndex.end())
		Die("Not found >%s", Label2.c_str());
	uint ChainIdx1 = iter1->second;
	uint ChainIdx2 = iter2->second;
	asserta(ChainIdx1 < SIZE(m_Chains));
	asserta(ChainIdx2 < SIZE(m_Chains));
	const PDBChain &Chain1 = *m_Chains[ChainIdx1];
	const PDBChain &Chain2 = *m_Chains[ChainIdx2];
	uint L1 = SIZE(Chain1.m_Seq);
	uint L2 = SIZE(Chain2.m_Seq);
	uint UngappedL1 = GetUngappedLength(Row1);
	uint UngappedL2 = GetUngappedLength(Row2);
	asserta(UngappedL1 == L1);
	asserta(UngappedL2 == L2);
	m_AlnChainIdxs.push_back(ChainIdx1);
	m_AlnChainIdxs.push_back(ChainIdx2);
	uint Opens, Exts;
	GetGapCounts(Row1, Row2, Opens, Exts);
	m_AlnOpenCounts.push_back(Opens);
	m_AlnExtCounts.push_back(Exts);
	m_AlnTPs.push_back(TP);
	if (TP)
		++m_TPCount;
	else
		++m_FPCount;
	asserta(SIZE(m_AlnTPs) == m_TPCount + m_FPCount);
	}

void FeatureTrainer::ReadAlns(const string &FN, bool TP)
	{
	SFasta SF;
	SF.Open(FN);
	SF.m_AllowGaps = true;

	bool Row1 = true;
	string FirstLabel = "";
	string FirstRow;
	for (;;)
		{
		const char* Seq = SF.GetNextSeq();
		if (Seq == 0)
			break;
		string Label = SF.GetLabel();
		TruncLabel(Label);
		const unsigned L = SF.GetSeqLength();
		asserta(L != 0);
		string Row;
		for (unsigned i = 0; i < L; ++i)
			Row.push_back(Seq[i]);
		m_AlnLabels.push_back(Label);
		m_AlnRows.push_back(Row);
		if (Row1)
			{
			FirstLabel = Label;
			FirstRow = Row;
			}
		else
			{
			const string &OtherRow = m_AlnRows.back();
			AddAln(FirstLabel, FirstRow, Label, Row, TP);
			}
		Row1 = !Row1;
		}
	const uint N = SIZE(m_AlnRows);
	asserta(N%2 == 0);
	asserta(SIZE(m_AlnChainIdxs) == N);
	asserta(SIZE(m_AlnTPs) == N/2);
	asserta(SIZE(m_AlnOpenCounts) == N/2);
	asserta(SIZE(m_AlnExtCounts) == N/2);
	}

void FeatureTrainer::SetFeature(FEATURE F, uint AlphaSize)
	{
	if (optset_maxi8)
		{
		uint Max = opt(maxi8);
		m_MaxAbsi8 = uint8_t(Max);
		asserta(uint(m_MaxAbsi8) == Max);
		}

	m_F = F;
	m_FeatureName = FeatureToStr(F);
	m_IsInt = FeatureIsInt(F);
	LogOdds::Init(AlphaSize);
	}

void FeatureTrainer::SetUnalignedBackgroundChain(const PDBChain &Chain,
	uint UndefLetter)
	{
	m_D.Init(Chain);
	const uint L = Chain.GetSeqLength();
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		uint Letter = UINT_MAX;
		if (m_IsInt)
			Letter = m_D.GetFeature(m_F, Pos);
		else
			{
			float Value = m_D.GetFloatFeature(m_F, Pos);

			// Allow undefined letters, AddUnalignedLetter will ignore
			Letter = ValueToIntTpl<true>(
				Value, m_AlphaSize, m_BinTs, UndefLetter);
			}
		AddUnalignedLetter(Letter);
		}
	}

void FeatureTrainer::SetUnalignedBackground(bool IgnoreUndef)
	{
	const uint ChainCount = SIZE(m_Chains);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Unaligned background");
		const PDBChain &Chain = *m_Chains[ChainIdx];
		SetUnalignedBackgroundChain(Chain, IgnoreUndef);
		}
	}

void FeatureTrainer::SetChainLetterSeqs(uint UndefLetter)
	{
	if (m_IsInt)
		SetChainLetterSeqs_Int(UndefLetter);
	else
		SetChainLetterSeqs_Float(UndefLetter);
	}

void FeatureTrainer::SetChainLetterSeqs_Int(uint UndefLetter)
	{
	m_ChainLetterSeqVec.clear();
	const uint ChainCount = SIZE(m_Chains);
	m_ChainLetterSeqVec.resize(ChainCount);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Letters %s", m_FeatureName);
		vector<uint> &ChainLetterSeq = m_ChainLetterSeqVec[ChainIndex];
		const PDBChain &Chain = *m_Chains[ChainIndex];
		m_D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		ChainLetterSeq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint Letter = m_D.GetFeature(m_F, Pos);
			asserta(Letter < m_AlphaSize || Letter == UINT_MAX);
			if (Letter == UINT_MAX)
				Letter = UndefLetter;
			ChainLetterSeq.push_back(Letter);
			}
		}
	}

void FeatureTrainer::SetChainLetterSeqs_Float(uint UndefLetter)
	{
	m_ChainLetterSeqVec.clear();
	const uint ChainCount = SIZE(m_Chains);
	m_ChainLetterSeqVec.resize(ChainCount);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Letters %s", m_FeatureName);
		vector<float> &ChainFloatSeq = m_ChainFloatSeqVec[ChainIndex];
		vector<uint> &ChainLetterSeq = m_ChainLetterSeqVec[ChainIndex];

		const PDBChain &Chain = *m_Chains[ChainIndex];
		const uint L = Chain.GetSeqLength();
		ChainLetterSeq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = ChainFloatSeq[Pos];
			uint Letter = ValueToIntTpl<true>(
				Value, m_AlphaSize, m_BinTs, UINT_MAX);
			asserta(Letter < m_AlphaSize || Letter == UINT_MAX);
			if (Letter == UINT_MAX)
				ChainLetterSeq.push_back(UndefLetter);
			else
				ChainLetterSeq.push_back(Letter);
			}
		}
	}

void FeatureTrainer::SetChainFloatSeqs(float ReplaceUndefValue)
	{
	m_ChainFloatSeqVec.clear();

	const uint ChainCount = SIZE(m_Chains);
	m_ChainFloatSeqVec.resize(ChainCount);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		vector<float> &ChainFloatSeq = m_ChainFloatSeqVec[ChainIndex];
		ProgressStep(ChainIndex, ChainCount, "Values %s", m_FeatureName);
		const PDBChain &Chain = *m_Chains[ChainIndex];
		m_D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		ChainFloatSeq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = m_D.GetFloatFeature(m_F, Pos);
			if (Value == FLT_MAX)
				Value = ReplaceUndefValue;
			ChainFloatSeq.push_back(Value);
			}
		}
	}

void FeatureTrainer::SetFloatValues(bool IgnoreUndef,
	float ReplaceUndefValue)
	{
	m_SortedFloatValues.clear();

	const uint ChainCount = SIZE(m_Chains);
	uint UndefCount = 0;
	uint ReplaceCount = 0;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Values %s", m_FeatureName);
		const PDBChain &Chain = *m_Chains[ChainIndex];
		m_D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = m_D.GetFloatFeature(m_F, Pos);
			if (Value == FLT_MAX)
				{
				++UndefCount;
				if (IgnoreUndef)
					continue;
				Value = ReplaceUndefValue;
				++ReplaceCount;
				}
			m_SortedFloatValues.push_back(Value);
			}
		}
	asserta(SIZE(m_SortedFloatValues) > 100);
	sort(m_SortedFloatValues.begin(), m_SortedFloatValues.end());

	const uint N = SIZE(m_SortedFloatValues);
	float Min = m_SortedFloatValues[0];
	float Med = m_SortedFloatValues[N/2];
	float Max = m_SortedFloatValues[N-1];

	uint MinCount = 0;
	for (uint i = 0; i < N; ++i)
		{
		if (m_SortedFloatValues[i] == Min)
			++MinCount;
		else
			break;
		}
	double MinPct = GetPct(MinCount, N);

	uint MaxCount = 0;
	for (uint i = 0; i < N; ++i)
		{
		if (m_SortedFloatValues[N-i-1] == Max)
			++MaxCount;
		else
			break;
		}
	double MaxPct = GetPct(MinCount, N);

	Log("\n");
	Log("SetFloatValues(IgnoreUndef=%c, ReplaceUndefValue=", tof(IgnoreUndef));
	if (ReplaceUndefValue == FLT_MAX)
		Log("*)\n");
	else
		Log("%.3g\n", ReplaceUndefValue);
	Log("%7u  Values\n", N);
	Log("%7u  Undefined\n", UndefCount);
	Log("%7u  Replaced\n", ReplaceCount);
	Log("%7.2g  Min (%.1f%%)\n", Min, MinPct);
	Log("%7.2g  Med\n", Med);
	Log("%7.2g  Max (%.1f%%)\n", Max, MaxPct);
	}

void FeatureTrainer::SetLabelToChainIndex()
	{
	m_LabelToChainIndex.clear();
	const uint ChainCount = SIZE(m_Chains);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *m_Chains[ChainIndex];
		string Label = Chain.m_Label;
		TruncLabel(Label);
		m_LabelToChainIndex[Label] = ChainIndex;
		}
	}

void FeatureTrainer::UpdateJointCounts(uint PairIndex, bool IgnoreUndef)
	{
	asserta(PairIndex < SIZE(m_AlnTPs));
	if (!m_AlnTPs[PairIndex])
		return;

	asserta(2*PairIndex+1 < SIZE(m_AlnLabels));
	asserta(2*PairIndex+1 < SIZE(m_AlnRows));

	uint QChainIdx = m_AlnChainIdxs[2*PairIndex];
	uint RChainIdx = m_AlnChainIdxs[2*PairIndex+1];

	const string &QRow = m_AlnRows[2*PairIndex];
	const string &RRow = m_AlnRows[2*PairIndex+1];
	const uint ColCount = SIZE(QRow);
	asserta(SIZE(RRow) == ColCount);

	const vector<uint> &QLetterSeq = m_ChainLetterSeqVec[QChainIdx];
	const vector<uint> &RLetterSeq = m_ChainLetterSeqVec[RChainIdx];

	uint QPos = 0;
	uint RPos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char q = QRow[Col];
		char r = RRow[Col];
		if (isunaligned(q))
			asserta(isunaligned(r));
		if (isupper(q) && isupper(r))
			{
			uint LetterQ = QLetterSeq[QPos];
			uint LetterR = RLetterSeq[RPos];
			if (IgnoreUndef)
				{
				if (LetterQ != UINT_MAX && LetterR != UINT_MAX)
					{
					m_Lock.lock();
					++m_JointLetterPairCount;
					AddPair(LetterQ, LetterR);
					m_Lock.unlock();
					}
				}
			else
				{
				asserta(LetterQ != UINT_MAX);
				asserta(LetterR != UINT_MAX);
				m_Lock.lock();
				++m_JointLetterPairCount;
				AddPair(LetterQ, LetterR);
				m_Lock.unlock();
				}
			}
		if (!isgap(q))
			++QPos;
		if (!isgap(r))
			++RPos;
		}
	asserta(QPos == SIZE(QLetterSeq));
	asserta(RPos == SIZE(RLetterSeq));
	}

void FeatureTrainer::TrainLogOdds(bool IgnoreUndef)
	{
	ResetCountsToZero();
	m_JointLetterPairCount = 0;
	const uint SeqCount = SIZE(m_AlnLabels);
	asserta(SeqCount > 0);
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
//#pragma omp parallel for
	for (int PairIndex = 0; PairIndex < int(PairCount); ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Joint frequencies %s(%u)",
					 m_FeatureName, m_AlphaSize);
		UpdateJointCounts(PairIndex, IgnoreUndef);
		}
	if (m_UseUnalignedBackground)
		SetUnalignedBackground(IgnoreUndef);

	Log("Joint counts\n");
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		const vector<uint> &Row = m_TrueCountMx[i];
		Log("[%2u] ", i);
		for (uint j = 0; j < m_AlphaSize; ++j)
			Log("%10u", Row[j]);
		Log("\n");
		}

	GetLogOddsMx(m_ScoreMx);
	}

void FeatureTrainer::WriteSummary(FILE *f) const
	{
	if (f == 0)
		return;
	float ES = LogOdds::GetExpectedScore();
	fprintf(f, "%s(%u)", m_FeatureName, m_AlphaSize);
	fprintf(f, " ES=%.3f", ES);
	fprintf(f, "\n");
	}

void FeatureTrainer::BinTsToSrc(FILE *f) const
	{
	if (f == 0)
		return;
	if (m_IsInt)
		return;
	if (m_BinTs.empty())
		return;
	const uint N = SIZE(m_BinTs);
	fprintf(f, "\nuint DSS::ValueToInt_%s(float Value) const\n", m_FeatureName);
	fprintf(f, "	{\n");
	for (uint i = 0 ; i < N; ++i)
		fprintf(f, "	if (Value < %.4g) return %u;\n", m_BinTs[i], i);
	fprintf(f, "	return %u;\n", N);
	fprintf(f, "	}\n");
	}

void FeatureTrainer::FreqsToSrc(FILE *f) const
	{
	if (f == 0)
		return;
	vector<float> Freqs;
	GetFreqs(Freqs);
	asserta(SIZE(Freqs) == m_AlphaSize);

	fprintf(f, "\nstatic float %s_f_i[%u] = {\n",
			m_FeatureName, m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		fprintf(f, "	%.4f,\n", Freqs[i]);
	fprintf(f, "};\n");
	}

void FeatureTrainer::ScoreMxToSrc(FILE *f) const
	{
	if (f == 0)
		return;
	vector<vector<float> > ScoreMx;
	GetLogOddsMx(ScoreMx);
	asserta(SIZE(ScoreMx) == m_AlphaSize);
	fprintf(f, "\nstatic float %s_S_i[%u][%u] = {\n",
			m_FeatureName, m_AlphaSize, m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		fprintf(f, "	{");
		for (uint j = 0; j < m_AlphaSize; ++j)
			fprintf(f, " %11.4g,", ScoreMx[i][j]);
		fprintf(f, " }, // %u\n", i);
		}
	fprintf(f, "};\n");
	}

void FeatureTrainer::ToSrc(FILE *f) const
	{
	if (f == 0)
		return;
	BinTsToSrc(f);
	FreqsToSrc(f);
	ScoreMxToSrc(f);
	}

void FeatureTrainer::ToTsv(const string &FN) const
	{
	if (FN == "")
		return;

	FILE *f = CreateStdioFile(FN);
	fprintf(f, "feature\t%s\n", m_FeatureName);
	fprintf(f, "type\t%s\n", m_IsInt ? "int" : "float");
	fprintf(f, "undef\t%s\n", m_Style.c_str());
	if (!m_IsInt)
		{
		float MaxValue = GetMaxDefinedValue();
		const uint N = SIZE(m_SortedFloatValues);
		fprintf(f, "min\t%.3g\n", m_SortedFloatValues[0]);
		fprintf(f, "med\t%.3g\n", m_SortedFloatValues[N/2]);
		fprintf(f, "max\t%.3g\n", MaxValue);
		}
	uint BestUndefLetter = GetBestDefaultLetter(UINT_MAX);
	fprintf(f, "best_undef\t%u\n", BestUndefLetter);

	LogOdds::ToTsv(f);
	
	if (!m_IsInt)
		{
		for (uint i = 0; i < SIZE(m_BinTs); ++i)
			fprintf(f, "bint\t%u\t%.4g\n", i, m_BinTs[i]);
		}

	CloseStdioFile(f);
	}

void FeatureTrainer::FromTsv(const string &FN)
	{
#if 0
	FILE *f = OpenStdioFile(FN);

	string FeatureName;
	ReadStringValue(f, "feature", FeatureName);
	m_FeatureName = mystrsave(FeatureName.c_str());
	m_F = StrToFeature(FeatureName.c_str());

	string Type;
	ReadStringValue(f, "type", Type);
	if (Type == "int")
		m_IsInt = true;
	else if (Type == "float")
		m_IsInt = false;
	else
		Die("Bad feature type '%s'", Type.c_str());

	string strUB;
	ReadStringValue(f, "undef_binning", strUB);
	m_UB = StrToUB(strUB.c_str());

	if (!m_IsInt)
		{
		m_MinValue = ReadFloatValue(f, "min");
		m_MedValue = ReadFloatValue(f, "med");
		m_MaxValue = ReadFloatValue(f, "max");
		m_UndefFreq = ReadFloatValue(f, "undef");
		}

	m_BestDefaultLetter = UINT_MAX;
	if (m_UB == UB_UndefinedIsDefaultLetter)
		m_BestDefaultLetter = ReadIntValue(f, "default_letter", UINT_MAX);

	LogOdds::FromTsv(f);

	m_BinTs.clear();
	if (!m_IsInt)
		{
		uint BinThresholdCount =
			DSS::GetBinThresholdCount(m_AlphaSize, m_UB);
		for (uint i = 0; i < BinThresholdCount; ++i)
			m_BinTs.push_back(ReadFloatValue(f, "bint", i));
		}
	ScoreMxFromTsv(f);

	string Line;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(!Ok);
	CloseStdioFile(f);
#endif
	}

void FeatureTrainer::ScoreMxFromTsv(FILE* f)
	{
	m_ScoreMx.resize(m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		ReadFloatVec(f, "scoremx", i, m_ScoreMx[i]);
	}

float FeatureTrainer::GetMaxDefinedValue() const
	{
	float MaxValue = -FLT_MAX;
	for (uint i = 0; i < SIZE(m_SortedFloatValues); ++i)
		{
		float Value = m_SortedFloatValues[i];
		if (Value != FLT_MAX)
			MaxValue = max(Value, MaxValue);
		}
	return MaxValue;
	}

float FeatureTrainer::GetMidValue(uint Letter) const
	{
	asserta(Letter < m_AlphaSize || Letter == UINT_MAX);
	if (Letter == UINT_MAX)
		{
		Log("GetMidValue(UINT_MAX)=FLT_MAX\n");
		return FLT_MAX;
		}
	asserta(SIZE(m_BinTs) + 1 == m_AlphaSize);
	if (Letter == m_AlphaSize - 1)
		return m_BinTs[m_AlphaSize - 2] + 1;
	float Lo = m_BinTs[Letter];
	float Hi = m_BinTs[Letter+1];
	float Value = (Lo + Hi)/2;
	Log("GetMidValue()\n");
	Log("  Letter=%u Lo=%.3g Hi=%.3g, Value=%.3g\n",
		Letter, Lo, Hi, Value);
	return Value;
	}

void FeatureTrainer::TrainInt_IgnoreUndefs()
	{
	SetChainLetterSeqs_Int(UINT_MAX);

	// Construct scoring matrix ignoring undefineds
	TrainLogOdds(true);
	}

// Undefined value overloads a letter
void FeatureTrainer::TrainInt_UndefOverlap(bool Retrain)
	{
	Log("\n");
	Log("TrainInt_UndefOverlap(Retrain=%c)\n", tof(Retrain));
	SetChainLetterSeqs_Int(UINT_MAX);

	// Construct scoring matrix ignoring undefineds
	TrainLogOdds(true);

	if (!Retrain)
		return;

	uint BestUndefLetterBefore = GetBestUndefLetter();
	vector<float> BeforeFreqs;
	GetFreqs(BeforeFreqs);

	uint UndefCountBeforeRetrain = GetUndefChainLetterCount();
	asserta(UndefCountBeforeRetrain > 0);

	// Best letter has highest expected score vs other letters
	uint UndefLetter = GetBestUndefLetter();
	asserta(UndefLetter != UINT_MAX);

	SetChainLetterSeqs_Int(UndefLetter);
	uint UndefCountAfterRetrain = GetUndefChainLetterCount();
	uint JNBefore = m_JointLetterPairCount;
	asserta(UndefCountAfterRetrain == 0);
	TrainLogOdds(false);
	uint JNAfter = m_JointLetterPairCount;

	uint BestUndefLetterAfter = GetBestUndefLetter();
	vector<float> AfterFreqs;
	GetFreqs(AfterFreqs);
	Log("\n");
	Log("Undefs before %u\n", UndefCountBeforeRetrain);
	Log("JNs before %u, after %u\n", JNBefore, JNAfter);
	Log("Best letter before %u, after %u\n",
		BestUndefLetterBefore, BestUndefLetterAfter);
	Log("Letter   Freq1   Freq2  before/after retrain\n");
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		Log("%6u  %6.4f  %6.4f\n",
			i, BeforeFreqs[i], AfterFreqs[i]);
		}
	}

uint FeatureTrainer::GetBestUndefLetter() const
	{
	vector<float> ExpectedScores;
	GetExpectedScores(ExpectedScores);
	asserta(SIZE(ExpectedScores) == m_AlphaSize);

	uint BestLetter = UINT_MAX;
	float BestES = -1;
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		float ES = ExpectedScores[Letter];
		if (ES > BestES)
			{
			BestLetter = Letter;
			BestES = ES;
			}
		}

	Log("\n");
	Log("Expected scores\n");
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		float ES = ExpectedScores[Letter];
		Log("  [%2u]", Letter);
		Log("  %8.3g", ES);
		if (ES == BestES)
			Log("  <== best");
		Log("\n");
		}
	return BestLetter;
	}

void FeatureTrainer::TrainFloat_IgnoreUndefs()
	{
	Log("\n");
	Log("TrainFloat_IgnoreUndefs()\n");

	// Collect only defined values
	SetFloatValues(true, FLT_MAX);

	// Quantize defined values
	Quantize(m_SortedFloatValues, m_AlphaSize, m_BinTs);
	SetChainFloatSeqs(FLT_MAX);
	SetChainLetterSeqs_Float(UINT_MAX);

	// Construct scoring matrix ignoring undefineds
	TrainLogOdds(true);
	}

// Undefined value overloads a letter
void FeatureTrainer::TrainFloat_UndefOverlap(bool Retrain)
	{
	Log("\n");
	Log("TrainFloat_UndefOverlap()\n");

	// Collect only defined values
	SetFloatValues(true, FLT_MAX);

	// Quantize defined values
	Quantize(m_SortedFloatValues, m_AlphaSize, m_BinTs);
	SetChainFloatSeqs(FLT_MAX);
	SetChainLetterSeqs_Float(UINT_MAX);

	// Construct scoring matrix ignoring undefineds
	TrainLogOdds(true);
	if (!Retrain)
		return;

	vector<float> BeforeFreqs;
	GetFreqs(BeforeFreqs);
	uint BestUndefLetterBefore = GetBestUndefLetter();

	uint UndefCountBeforeRetrain = GetUndefChainLetterCount();
	asserta(UndefCountBeforeRetrain > 0);

	// Best letter has highest expected score vs other letters
	uint UndefLetter = GetBestUndefLetter();
	asserta(UndefLetter != UINT_MAX);

	// Default value is midpoint of bin for m_BestDefaultLetter
	float UndefValue = GetMidValue(UndefLetter);
	asserta(UndefValue != FLT_MAX);

	SetFloatValues(false, UndefValue);
	SetChainFloatSeqs(UndefValue);
	SetChainLetterSeqs_Float(UndefLetter);
	uint UndefCountAfterRetrain = GetUndefChainLetterCount();
	asserta(UndefCountAfterRetrain == 0);
	TrainLogOdds(false);
	uint BestUndefLetterAfter = GetBestUndefLetter();
	vector<float> AfterFreqs;
	GetFreqs(AfterFreqs);

	Log("\n");
	Log("Best letter before %u, after %u\n",
		BestUndefLetterBefore, BestUndefLetterAfter);
	Log("Letter   Freq1   Freq2  before/after retrain\n");
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		Log("%6u  %6.4f  %6.4f\n",
			i, BeforeFreqs[i], AfterFreqs[i]);
		}
	}

// Undefined value has its own letter (value AS-1),
//   leaving AS-1 other letters for defined values
void FeatureTrainer::TrainFloat_UndefDistinct()
	{
	uint UndefLetter = m_AlphaSize - 1;
	float UndefValue = FLT_MAX;

	// Collect all values
	SetFloatValues(false, FLT_MAX);

	// Quantize all values, undefineds will get their own bin
	//   without special-casing
	Quantize(m_SortedFloatValues, m_AlphaSize, m_BinTs);
	SetChainFloatSeqs(FLT_MAX);
	SetChainLetterSeqs_Float(UndefLetter);

	// Train scoring matrix including undefineds
	TrainLogOdds(false);
	}

void FeatureTrainer::SetAlnSubstScores()
	{
	m_AlnSubstScores.clear();
	const uint AlnCount = SIZE(m_AlnTPs);
	m_AlnSubstScores.reserve(AlnCount);
	asserta(AlnCount > 0);
	vector<float> TPScores;
	vector<float> FPScores;
	for (int AlnIdx = 0; AlnIdx < int(AlnCount); ++AlnIdx)
		{
		ProgressStep(AlnIdx, AlnCount, "Subst scores");
		float Score = GetAlnSubstScore(AlnIdx);
		m_AlnSubstScores.push_back(Score);
		if (m_AlnTPs[AlnIdx])
			TPScores.push_back(Score);
		else
			FPScores.push_back(Score);
		}

	Log("\n");
	Log("SetAlnSubstScores()\n");
	QuartsFloat Q;
	GetQuartsFloat(TPScores, Q);
	Log("TP scores: ");
	Q.LogMe();
	if (!FPScores.empty())
		{
		GetQuartsFloat(FPScores, Q);
		Log("FP scores: ");
		Q.LogMe();
		}
	}

float FeatureTrainer::GetAlnSubstScore(uint AlnIdx)
	{
	asserta(SIZE(m_ScoreMx) == m_AlphaSize);
	asserta(2*AlnIdx+1 < SIZE(m_AlnLabels));
	asserta(2*AlnIdx+1 < SIZE(m_AlnRows));

	uint QChainIdx = m_AlnChainIdxs[2*AlnIdx];
	uint RChainIdx = m_AlnChainIdxs[2*AlnIdx+1];

	const string &QRow = m_AlnRows[2*AlnIdx];
	const string &RRow = m_AlnRows[2*AlnIdx+1];
	const uint ColCount = SIZE(QRow);
	asserta(SIZE(RRow) == ColCount);

	const vector<uint> &QLetterSeq = m_ChainLetterSeqVec[QChainIdx];
	const vector<uint> &RLetterSeq = m_ChainLetterSeqVec[RChainIdx];

	uint QPos = 0;
	uint RPos = 0;
	float Score = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char q = QRow[Col];
		char r = RRow[Col];
		if (isunaligned(q))
			asserta(isunaligned(r));
		if (isupper(q) && isupper(r))
			{
			uint LetterQ = QLetterSeq[QPos];
			uint LetterR = RLetterSeq[RPos];
			if (LetterQ != UINT_MAX && LetterR != UINT_MAX)
				{
				asserta(SIZE(m_ScoreMx[LetterQ]) == m_AlphaSize);
				Score += m_ScoreMx[LetterQ][LetterR];
				}
			}
		if (!isgap(q))
			++QPos;
		if (!isgap(r))
			++RPos;
		}
	asserta(QPos == SIZE(QLetterSeq));
	asserta(RPos == SIZE(RLetterSeq));
	return Score;
	}

void FeatureTrainer::SetAlnScoresAndArea(float OpenPenalty, float ExtPenalty)
	{
	const uint AlnCount = SIZE(m_AlnTPs);
	asserta(SIZE(m_AlnSubstScores) == AlnCount);
	asserta(SIZE(m_AlnOpenCounts) == AlnCount);
	asserta(SIZE(m_AlnExtCounts) == AlnCount);
	m_AlnScores.clear();
	m_AlnScores.reserve(AlnCount);
	for (uint AlnIdx = 0; AlnIdx < AlnCount; ++AlnIdx)
		{
		float Score = m_AlnSubstScores[AlnIdx];
		Score -= m_AlnOpenCounts[AlnIdx]*OpenPenalty;
		Score -= m_AlnExtCounts[AlnIdx]*ExtPenalty;
		Score = round3sigfig(Score);
		m_AlnScores.push_back(Score);
		}
	SetArea();
	}

void FeatureTrainer::LogROCStepsAndArea()
	{
	m_Area = 0;
	m_ROCStepScores.clear();
	m_ROCStepTPfs.clear();
	m_ROCStepFPfs.clear();
	m_ROCOrder.clear();

	const uint AlnCount = GetAlnCount();
	asserta(AlnCount == m_TPCount + m_FPCount);
	asserta(m_TPCount > 0);
	asserta(m_FPCount > 0);

	m_ROCStepScores.reserve(AlnCount);
	m_ROCStepTPfs.reserve(AlnCount);
	m_ROCStepFPfs.reserve(AlnCount);

	m_ROCOrder.resize(AlnCount);
	QuickSortOrderDesc(m_AlnScores.data(), AlnCount, m_ROCOrder.data());

	uint NTP = 0;
	uint NFP = 0;
	float PrevTPf = 0;
	float PrevFPf = 0;
	float CurrentScore = m_AlnScores[m_ROCOrder[0]];
	for (uint k = 0; k < AlnCount; ++k)
		{
		uint i = m_ROCOrder[k];
		float Score = m_AlnScores[i];
		bool TP = m_AlnTPs[i];
		if (Score != CurrentScore)
			{
			float TPf = float(NTP)/m_TPCount;
			float FPf = float(NFP)/m_FPCount;
			m_Area += TPf*(FPf - PrevFPf);
			m_ROCStepScores.push_back(CurrentScore);
			m_ROCStepTPfs.push_back(TPf);
			m_ROCStepFPfs.push_back(FPf);
			CurrentScore = Score;
			PrevTPf = TPf;
			PrevFPf = FPf;
			}
		if (TP)
			++NTP;
		else
			++NFP;
		}
	asserta(NTP == m_TPCount);
	asserta(NFP == m_FPCount);
	asserta(NTP > 0);
	asserta(NFP > 0);

	float TPf = float(NTP)/m_TPCount;
	float FPf = float(NFP)/m_FPCount;
	m_ROCStepScores.push_back(CurrentScore);
	m_ROCStepTPfs.push_back(TPf);
	m_ROCStepFPfs.push_back(FPf);
	m_Area += TPf*(FPf - PrevFPf);

	const uint StepCount = SIZE(m_ROCStepTPfs);
	asserta(SIZE(m_ROCStepFPfs) == StepCount);
	asserta(SIZE(m_ROCStepScores) == StepCount);

	Log("\n");
	Log("SetROCSteps() %u steps\n", StepCount);
	Log("%8.8s  %10.10s  %10.10s\n", "Score", "NTP", "NFP");
	for (uint StepIdx = 0; StepIdx < StepCount; ++StepIdx)
		{
		float Score = m_ROCStepScores[StepIdx];
		float TPf = m_ROCStepTPfs[StepIdx];
		float FPf = m_ROCStepFPfs[StepIdx];
		Log("%.3g\t%.3g\t%.3g\n", Score, TPf, FPf);
		}
	Log("Area = %.3g\n", m_Area);
	}

void FeatureTrainer::SetArea()
	{
	m_Area = 0;
	const uint AlnCount = GetAlnCount();
	uint *Order = myalloc(uint, AlnCount);
	QuickSortOrderDesc(m_AlnScores.data(), AlnCount, Order);

	uint NTP = 0;
	uint NFP = 0;
	float PrevTPf = 0;
	float PrevFPf = 0;
	float CurrentScore = m_AlnScores[Order[0]];
	for (uint k = 0; k < AlnCount; ++k)
		{
		uint i = Order[k];
		float Score = m_AlnScores[i];
		bool TP = m_AlnTPs[i];
		if (Score != CurrentScore)
			{
			float TPf = float(NTP)/m_TPCount;
			float FPf = float(NFP)/m_FPCount;
			m_Area += (PrevTPf + TPf)*(FPf - PrevFPf)/2;
			CurrentScore = Score;
			PrevTPf = TPf;
			PrevFPf = FPf;
			}
		if (TP)
			++NTP;
		else
			++NFP;
		}
	myfree(Order);
	Order = 0;
	asserta(NTP == m_TPCount);
	asserta(NFP == m_FPCount);
	asserta(NTP > 0);
	asserta(NFP > 0);

	float TPf = float(NTP)/m_TPCount;
	float FPf = float(NFP)/m_FPCount;
	m_Area += (PrevTPf + TPf)*(FPf - PrevFPf)/2;
	}

static FeatureTrainer *s_FT;
double FeatureTrainer::EvalArea(const vector<double> &xv)
	{
	asserta(s_FT != 0);
	asserta(SIZE(xv) == 2);
	float OpenPenalty = float(xv[0]);
	float ExtPenalty = float(xv[1]);
	s_FT->SetAlnScoresAndArea(OpenPenalty, ExtPenalty);
	return s_FT->m_Area;
	}

uint FeatureTrainer::GetUndefChainLetterCount() const
	{
	uint n = 0;
	for (uint ChainIdx = 0; ChainIdx < SIZE(m_ChainLetterSeqVec); ++ChainIdx)
		{
		const vector<uint> &Letters = m_ChainLetterSeqVec[ChainIdx];
		uint L = SIZE(Letters);
		for (uint i = 0; i < L; ++i)
			if (Letters[i] == UINT_MAX)
				++n;
		}
	return n;
	}

void FeatureTrainer::OptimizeGapPenalties()
	{
	s_FT = this;

	vector<string> SpecLines;

	// y is Area, in range 0 to 1
	SpecLines.push_back("mindy=0.001");
	SpecLines.push_back("maxdy=0.1");
	SpecLines.push_back("minh=0.0005");
	SpecLines.push_back("latin=yes");
	SpecLines.push_back("sigfig=3");
	SpecLines.push_back("var=open\tmin=-2\tmax=4\tdelta=0.2\tbins=32\tinit=1");
	SpecLines.push_back("var=ext\tmin=-1\tmax=1\tdelta=0.04\tbins=32\tinit=0.1");

	m_Peaker.Init(SpecLines, EvalArea);
	m_Peaker.Run();
	m_Area = (float) m_Peaker.m_Best_y;
	}

void FeatureTrainer::Train(const string &Style)
	{
	m_Style = Style;

	bool IgnoreUndefs = GetRequiredStyleOpt(Style, 'I');
	m_UseUnalignedBackground = GetRequiredStyleOpt(Style, 'U');
	bool UndefOverlap = GetRequiredStyleOpt(Style, 'O');
	bool RetrainOverlap = GetRequiredStyleOpt(Style, 'R');

	Log("Style=%s\n", Style.c_str());
	Log("IgnoreUndefs=%c\n", tof(IgnoreUndefs));
	Log("UndefOverlap=%c\n", tof(UndefOverlap));
	Log("RetrainOverlap=%c\n", tof(RetrainOverlap));
	Log("m_UseUnalignedBackground=%c\n", tof(m_UseUnalignedBackground));

	if (m_IsInt)
		{
		asserta(IgnoreUndefs || UndefOverlap);
		if (IgnoreUndefs)
			TrainInt_IgnoreUndefs();
		else
			TrainInt_UndefOverlap(RetrainOverlap);
		}
	else
		{
		asserta(IgnoreUndefs || UndefOverlap);
		if (IgnoreUndefs)
			TrainFloat_IgnoreUndefs();
		else
			TrainFloat_UndefOverlap(RetrainOverlap);
		}
	}