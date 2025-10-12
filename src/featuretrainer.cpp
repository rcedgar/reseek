#include "myutils.h"
#include "featuretrainer.h"
#include "alpha.h"

static uint GetUngappedLength(const string &Row)
	{
	const uint ColCount = SIZE(Row);
	uint L = 0;
	for (uint i = 0; i < ColCount; ++i)
		if (!isgap(Row[i]))
			++L;
	return L;
	}

// d1a0rp__P
// 012345678
static void TruncLabel(string &Label)
	{
	//size_t k = Label.size();
	//if (k == 9 && Label[7] == '_' && isupper(Label[8]))
	//	{
	//	Label.resize(7);
	//	return;
	//	}
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

void FeatureTrainer::SetInput(const string &ChainsFN, const string &AlnsFN)
	{
// Params
	optset_fast = true;
	opt(fast) = true;
	DSSParams::Init(DM_DefaultFast);
	m_MaxAbsi8 = 20;
	if (optset_maxi8)
		{
		uint Max = opt(maxi8);
		m_MaxAbsi8 = uint8_t(Max);
		asserta(uint(m_MaxAbsi8) == Max);
		}

// Chains
	ReadChains(ChainsFN, m_Chains);
	SetLabelToChainIndex();

// Alignments
	m_Alns.FromFasta(AlnsFN, true);
	}

void FeatureTrainer::SetFeature(FEATURE F, uint AlphaSize)
	{
	asserta(!m_Chains.empty());
	asserta(m_Alns.GetSeqCount() > 0);
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
			Letter = DSS::ValueToInt(Value, m_AlphaSize,
				m_BinTs, UndefLetter);
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

void FeatureTrainer::SetFloatValues(bool IgnoreUndef,
	float ReplaceUndefValue)
	{
	m_UndefCount = 0;
	m_SortedFloatValues.clear();

	DSSParams::Init(DM_DefaultFast);
	const uint ChainCount = SIZE(m_Chains);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		uint N = SIZE(m_SortedFloatValues);
		ProgressStep(ChainIndex, ChainCount, "Values %s", m_FeatureName);
		const PDBChain &Chain = *m_Chains[ChainIndex];
		m_D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = m_D.GetFloatFeature(m_F, Pos);
			if (Value == FLT_MAX)
				{
				++m_UndefCount;
				if (IgnoreUndef)
					continue;
				Value = ReplaceUndefValue;
				}
			m_SortedFloatValues.push_back(Value);
			}
		}
	asserta(SIZE(m_SortedFloatValues) > 100);
	sort(m_SortedFloatValues.begin(), m_SortedFloatValues.end());
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

// >d1gaia_/a.102.1.1|E=0.908|Id=21.3%|TS=0.1274|AQ=0.6715|d3p2ca_/a.102.1.8
float FeatureTrainer::GetPctIdFromLabel(const string &Label) const
	{
	vector<string> Fields;
	Split(Label, Fields, '|');
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		const string &Field = Fields[i];
		if (StartsWith(Field, "Id="))
			{
			vector<string> Fields2;
			Split(Field, Fields2, '=');
			asserta(SIZE(Fields2) == 2);
			string sPctId = Fields2[1];
			if (EndsWith(sPctId, "%"))
				sPctId = sPctId.substr(0, sPctId.size()-1);
			float PctId = (float) StrToFloat(sPctId);
			asserta(PctId > 0 && PctId <= 100);
			return PctId;
			}
		}
	Die("Id= not found in >%s", Label.c_str());
	return FLT_MAX;
	}

float FeatureTrainer::GetAQFromLabel(const string &Label) const
	{
	vector<string> Fields;
	Split(Label, Fields, '|');
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		const string &Field = Fields[i];
		if (StartsWith(Field, "AQ="))
			{
			vector<string> Fields2;
			Split(Field, Fields2, '=');
			asserta(SIZE(Fields2) == 2);
			string sAQ = Fields2[1];
			float AQ = (float) StrToFloat(sAQ);
			asserta(AQ >= 0 && AQ <= 1);
			return AQ;
			}
		}
	Die("AQ= not found in >%s", Label.c_str());
	return FLT_MAX;
	}

void FeatureTrainer::UpdateJointCounts(uint PairIndex, bool IgnoreUndef)
	{
	string QLabel = m_Alns.GetLabel(2*PairIndex);
	string RLabel = m_Alns.GetLabel(2*PairIndex+1);
	TruncLabel(QLabel);
	TruncLabel(RLabel);

	map<string, uint>::const_iterator iterq = m_LabelToChainIndex.find(QLabel);
	map<string, uint>::const_iterator iterr = m_LabelToChainIndex.find(RLabel);
	if (iterq == m_LabelToChainIndex.end())
		Die("Not found >%s", QLabel.c_str());
	asserta(iterr != m_LabelToChainIndex.end());

	uint QChainIndex = iterq->second;
	uint RChainIndex = iterr->second;

	const PDBChain &QChain = *m_Chains[QChainIndex];
	const PDBChain &RChain = *m_Chains[RChainIndex];

	string QLabel2 = QChain.m_Label;
	string RLabel2 = RChain.m_Label;
	TruncLabel(QLabel2);
	TruncLabel(RLabel2);
	asserta(QLabel2 == QLabel);
	asserta(RLabel2 == RLabel);

	uint QL = QChain.GetSeqLength();
	uint RL = RChain.GetSeqLength();

	DSS DQ;
	DSS DR;
	DQ.Init(QChain);
	DR.Init(RChain);

	const string &QRow = m_Alns.GetSeq(2*PairIndex);
	const string &RRow = m_Alns.GetSeq(2*PairIndex+1);

	uint QL2 = GetUngappedLength(QRow);
	uint RL2 = GetUngappedLength(RRow);
	if (QL2 != QL || RL2 != RL)
		{
		Log("\n");
		Log("Q>%s\n", QChain.m_Label.c_str());
		Log("R>%s\n", RChain.m_Label.c_str());
		Log("Q: %s\n", QRow.c_str());
		Log("T: %s\n", RRow.c_str());
		Die("QL_row %u, QL %u, RL_row %u, RL %u",
			QL2, QL, RL2, RL);
		}
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
			if (m_IsInt)
				{
				LetterQ = DQ.GetFeature(m_F, QPos);
				LetterR = DR.GetFeature(m_F, RPos);
				
				if (IgnoreUndef)
					{
					asserta(LetterQ < m_AlphaSize || LetterQ == UINT_MAX);
					asserta(LetterR < m_AlphaSize || LetterR == UINT_MAX);
					if (LetterQ != UINT_MAX && LetterR != UINT_MAX)
						{
						m_Lock.lock();
						AddPair(LetterQ, LetterR);
						m_Lock.unlock();
						}
					}
				else
					{
					asserta(LetterQ < m_AlphaSize-1 || LetterQ == UINT_MAX);
					asserta(LetterR < m_AlphaSize-1 || LetterR == UINT_MAX);
					if (LetterQ == UINT_MAX)
						LetterQ = m_AlphaSize - 1;
					if (LetterR == UINT_MAX)
						LetterR = m_AlphaSize - 1;
					m_Lock.lock();
					AddPair(LetterQ, LetterR);
					m_Lock.unlock();
					}
				}
			else
				{
				float ValueQ = DQ.GetFloatFeature(m_F, QPos);
				float ValueR = DR.GetFloatFeature(m_F, RPos);

				if (IgnoreUndef)
					{
					bool IgnoreQ = (ValueQ == FLT_MAX);
					bool IgnoreR = (ValueR == FLT_MAX);
					if (!IgnoreQ && !IgnoreR)
						{
						LetterQ = DSS::ValueToInt(ValueQ, m_AlphaSize,
													m_BinTs, m_BestDefaultLetter);

						LetterR = DSS::ValueToInt(ValueR, m_AlphaSize,
													m_BinTs, m_BestDefaultLetter);
						m_Lock.lock();
						AddPair(LetterQ, LetterR);
						m_Lock.unlock();
						}
					}
				else
					{
					if (ValueQ == FLT_MAX)
						ValueQ = m_BestDefaultValue;						
					if (ValueR == FLT_MAX)
						ValueR = m_BestDefaultValue;						
					LetterQ = DSS::ValueToInt(ValueQ, m_AlphaSize,
												m_BinTs, m_BestDefaultLetter);

					LetterR = DSS::ValueToInt(ValueR, m_AlphaSize,
												m_BinTs, m_BestDefaultLetter);
					asserta(LetterQ < m_AlphaSize);
					asserta(LetterR < m_AlphaSize);
					m_Lock.lock();
					AddPair(LetterQ, LetterR);
					m_Lock.unlock();
					}
				}
			}
		if (!isgap(q))
			++QPos;
		if (!isgap(r))
			++RPos;
		}
	}

void FeatureTrainer::TrainLogOdds(bool IgnoreUndef)
	{
	const uint SeqCount = m_Alns.GetSeqCount();
	asserta(SeqCount > 0);
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	m_Counter = 0;
//#pragma omp parallel for
	for (int PairIndex = 0; PairIndex < int(PairCount); ++PairIndex)
		{
		m_Lock.lock();
		uint Count = m_Counter++;
		m_Lock.unlock();
		ProgressStep(Count, PairCount, "Joint frequencies %s(%u)",
					 m_FeatureName, m_AlphaSize);
		UpdateJointCounts(PairIndex, IgnoreUndef);
		}
	if (m_UseUnalignedBackground)
		SetUnalignedBackground(IgnoreUndef);
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
	fprintf(f, "undef\t%s\n", m_UndefStyle.c_str());
	if (!m_IsInt)
		{
		float MaxValue = GetMaxDefinedValue();
		const uint N = SIZE(m_SortedFloatValues);
		fprintf(f, "min\t%.3g\n", m_SortedFloatValues[0]);
		fprintf(f, "med\t%.3g\n", m_SortedFloatValues[N/2]);
		fprintf(f, "max\t%.3g\n", MaxValue);
		fprintf(f, "undef_count\t%u\n", m_UndefCount);
		fprintf(f, "undef_freq\t%.3g\n", double(m_UndefCount)/N);
		}
	if (m_BestDefaultLetter == UINT_MAX)
		fprintf(f, "default_letter\t*\n");
	else
		fprintf(f, "default_letter\t%u\n", m_BestDefaultLetter);

	if (m_BestDefaultValue == FLT_MAX)
		fprintf(f, "default_value\t*\n");
	else
		fprintf(f, "default_value\t%.3g\n", m_BestDefaultValue);

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

float FeatureTrainer::GetDefaultValue() const
	{
	if (m_BestDefaultLetter == UINT_MAX)
		return FLT_MAX;
	asserta(m_BestDefaultLetter < m_AlphaSize);
	asserta(SIZE(m_BinTs) + 1 == m_AlphaSize);
	if (m_BestDefaultLetter == m_AlphaSize - 1)
		return m_BinTs[m_AlphaSize - 2] + 1;
	float Lo = m_BinTs[m_BestDefaultLetter];
	float Hi = m_BinTs[m_BestDefaultLetter+1];
	return (Lo + Hi)/2;
	}

// Undefined value overloads a letter
void FeatureTrainer::TrainInt_UndefOverlap()
	{
	m_UndefStyle = "overlap";

	// Construct scoring matrix ignoring undefineds
	TrainLogOdds(true);

	// Best letter has highest expected score vs other letters
	m_BestDefaultLetter = LogOdds::GetBestDefaultLetter(UINT_MAX);
	m_BestDefaultValue = FLT_MAX;
	}

// Undefined value overloads a letter
void FeatureTrainer::TrainFloat_UndefOverlap()
	{
	m_UndefStyle = "overlap";
	// Collect only defined values
	SetFloatValues(true, FLT_MAX);

	// Quantize defined values
	Quantize(m_SortedFloatValues, m_AlphaSize, m_BinTs);

	// Construct scoring matrix ignoring undefineds
	TrainLogOdds(true);

	// Best letter has highest expected score vs other letters
	m_BestDefaultLetter = LogOdds::GetBestDefaultLetter(UINT_MAX);

	// Default value is midpoint of bin for m_BestDefaultLetter
	m_BestDefaultValue = GetDefaultValue();
	}

// Undefined value has its own letter (value AS-1),
//   leaving AS-1 other letters for defined values
void FeatureTrainer::TrainFloat_UndefDistinct()
	{
	m_UndefStyle = "distinct";

	// Collect all values
	SetFloatValues(false, FLT_MAX);

	// Quantize all values, undefineds will get their own bin
	//   without special-casing
	Quantize(m_SortedFloatValues, m_AlphaSize, m_BinTs);
	m_BestDefaultLetter = m_AlphaSize - 1;
	m_BestDefaultValue = FLT_MAX;

	// Train scoring matrix including undefineds
	TrainLogOdds(false);
	}

// Undefined value has its own letter (value AS-1),
//   leaving AS-1 other letters for defined values
void FeatureTrainer::TrainInt_UndefDistinct()
	{
	m_UndefStyle = "distinct";

	m_AlphaSize += 1;
	// Train scoring matrix including undefineds
	TrainLogOdds(false);

	m_BestDefaultLetter = m_AlphaSize - 1;
	m_BestDefaultValue = FLT_MAX;
	}

void FeatureTrainer::Quantize(const vector<float> &SortedValues,
	uint AlphaSize, vector<float> &BinTs)
	{
	BinTs.clear();

	const uint K = SIZE(SortedValues);
	asserta(K > 0);
	float MinValue = SortedValues[0];
	float MaxValue = SortedValues[K-1];
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		{
		uint k = ((i+1)*K)/AlphaSize;
		float t = SortedValues[k];
		if (i > 0)
			{
			if (feq(t, BinTs[i-1]))
				Log("Quantize tie for bin thresholds %u,%u at %.3g\n",
						i-1, i, t);
			QuantizeUniques(SortedValues, AlphaSize, BinTs);
			asserta(SIZE(BinTs) + 1 == AlphaSize);
			return;
			}
		BinTs.push_back(t);
		}
	asserta(SIZE(BinTs) + 1 == AlphaSize);
	}

void FeatureTrainer::QuantizeUniques(const vector<float> &SortedValues,
	uint AlphaSize, vector<float> &BinTs)
	{
	BinTs.clear();
	vector<float> UniqueFloatValues;
	vector<uint> UniqueFloatCounts;
	
	const uint N = SIZE(SortedValues);
	asserta(N > 100);
	float UniqueValue = SortedValues[0];
	uint Count = 1;
	for (uint i = 1; i < N; ++i)
		{
		float Value = SortedValues[i];
		if (Value == UniqueValue)
			++Count;
		else
			{
			asserta(Value > UniqueValue);
			UniqueFloatValues.push_back(UniqueValue);
			UniqueFloatCounts.push_back(Count);
			UniqueValue = Value;
			Count = 1;
			}
		}
	UniqueFloatValues.push_back(UniqueValue);
	UniqueFloatCounts.push_back(Count);
	uint UniqueValueCount = SIZE(UniqueFloatValues);
	asserta(SIZE(UniqueFloatCounts) == UniqueValueCount);
	asserta(UniqueValueCount >= AlphaSize);
	uint TargetCountPerBin = uint(double(N)/AlphaSize + 0.5);
	vector<float> TmpValues;
	for (uint i = 0; i < UniqueValueCount; ++i)
		{
		float Value = UniqueFloatValues[i];
		uint n = min(UniqueFloatCounts[i], TargetCountPerBin/2);
		for (uint j = 0; j < n; ++j)
			TmpValues.push_back(Value);
		}

	const uint K = SIZE(TmpValues);
	asserta(K > 0);
	float MinValue = TmpValues[0];
	float MaxValue = TmpValues[K-1];
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		{
		uint k = ((i+1)*K)/AlphaSize;
		float t = TmpValues[k];
		if (i > 0)
			{
			if (feq(t, BinTs[i-1]))
				Die("QuantizeUniques tie for bin thresholds %u,%u at %.3g\n",
						i-1, i, t);
			}
		BinTs.push_back(t);
		}
	asserta(SIZE(BinTs) + 1 == AlphaSize);
	}
