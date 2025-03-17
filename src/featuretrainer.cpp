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

static void TruncLabel(string &Label)
	{
	size_t n = Label.find(' ');
	if (n == string::npos || n == 0)
		return;
	Label.resize(n);
	}

void FeatureTrainer::SetInput(const string &ChainsFN, const string &AlnsFN)
	{
// Params
	optset_fast = true;
	opt(fast) = true;
	m_Params.SetDSSParams(DM_DefaultFast, SCOP40_DBSIZE);
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

void FeatureTrainer::SetFeature(FEATURE F)
	{
	asserta(!m_Chains.empty());
	asserta(m_Alns.GetSeqCount() > 0);
	m_F = F;
	m_FeatureName = FeatureToStr(F);
	m_IsInt = FeatureIsInt(F);
	if (!m_IsInt)
		SetFloatValues();
	}

void FeatureTrainer::SetAlphaSize(uint AS)
	{
	m_AlphaSize = AS;
	Init(AS);

	if (m_IsInt)
		return;

	asserta(!m_FloatValues.empty());
	DSS::Condense(m_FloatValues, m_AlphaSize, m_Wildcard,
					m_MinValue, m_MedValue, m_MaxValue,
					m_UndefFreq, m_BinTs);
	}

void FeatureTrainer::SetFloatValues()
	{
	m_MinValue = FLT_MAX;
	m_MedValue = FLT_MAX;
	m_MaxValue = -FLT_MAX;

	uint UndefCount = 0;
	m_FloatValues.clear();

	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast, SCOP40_DBSIZE);

	DSS D;
	D.SetParams(Params);

	const uint ChainCount = SIZE(m_Chains);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		uint N = SIZE(m_FloatValues);
		double Pct = GetPct(UndefCount, N);
		ProgressStep(ChainIndex, ChainCount, "Values %s (%.2f%% undefined)",
					 m_FeatureName, Pct);
		const PDBChain &Chain = *m_Chains[ChainIndex];
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = D.GetFloatFeature(m_F, Pos);
			if (Value == FLT_MAX)
				++UndefCount;
			else
				{
				m_MinValue = min(m_MinValue, Value);
				m_MaxValue = max(m_MaxValue, Value);
				}
			m_FloatValues.push_back(Value);
			}
		}
	m_UndefFreq = float(UndefCount)/float(SIZE(m_FloatValues));
	}

void FeatureTrainer::SetLabelToChainIndex()
	{
	m_LabelToChainIndex.clear();
	const uint ChainCount = SIZE(m_Chains);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *m_Chains[ChainIndex];
		const string &Label = Chain.m_Label;
		m_LabelToChainIndex[Label] = ChainIndex;
		}
	}

void FeatureTrainer::SetJointFreqsPair(uint PairIndex)
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

	uint QL = QChain.GetSeqLength();
	uint RL = RChain.GetSeqLength();

	m_DQ.Init(QChain);
	m_DR.Init(RChain);

	const string &QRow = m_Alns.GetSeq(2*PairIndex);
	const string &RRow = m_Alns.GetSeq(2*PairIndex+1);

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
			if (m_IsInt)
				{
				LetterQ = m_DQ.GetFeature(m_F, QPos);
				LetterR = m_DR.GetFeature(m_F, RPos);
				}
			else
				{
				float ValueQ = m_DQ.GetFloatFeature(m_F, QPos);
				float ValueR = m_DR.GetFloatFeature(m_F, RPos);

				LetterQ = DSS::ValueToInt(m_BinTs, ValueQ);
				LetterR = DSS::ValueToInt(m_BinTs, ValueR);
				}
			AddPair(LetterQ, LetterR);
			}
		if (!isgap(q))
			++QPos;
		if (!isgap(r))
			++RPos;
		}
	}

void FeatureTrainer::Train(bool Wildcard)
	{
	m_Wildcard = Wildcard;
	const uint SeqCount = m_Alns.GetSeqCount();
	asserta(SeqCount > 0);
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	uint LetterPairCount = 0;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Joint frequencies %s(%u)",
					 m_FeatureName, m_AlphaSize);
		SetJointFreqsPair(PairIndex);
		}
	m_BestDefaultLetter = GetBestDefaultLetter(UINT_MAX);
	}

void FeatureTrainer::WriteSummary(FILE *f) const
	{
	if (f == 0)
		return;
	float ES = GetExpectedScore();
	fprintf(f, "%s(%u) wc=%c", m_FeatureName, m_AlphaSize, yon(m_Wildcard));
	if (m_IsInt)
		fprintf(f, " integer");
	else
		{
		fprintf(f, ", condensed min %.3g, med %.3g, max %.3g, undef %.4f",
				m_MinValue,
				m_MedValue,
				m_MaxValue,
				m_UndefFreq);
		}
	fprintf(f, " ES=%.3f", ES);
	fprintf(f, "\n");
	}

void FeatureTrainer::ToTsv(const string &FN) const
	{
	if (FN == "")
		return;

	FILE *f = CreateStdioFile(FN);
	fprintf(f, "feature\t%s\n", m_FeatureName);
	fprintf(f, "type\t%s\n", m_IsInt ? "int" : "float");
	fprintf(f, "wildcard\t%s\n", m_Wildcard ? "yes" : "no");
	if (!m_IsInt)
		{
		fprintf(f, "min\t%.3g\n", m_MinValue);
		fprintf(f, "med\t%.3g\n", m_MedValue);
		fprintf(f, "max\t%.3g\n", m_MaxValue);
		fprintf(f, "undef\t%.4f\n", m_UndefFreq);
		}
	if (!m_Wildcard)
		fprintf(f, "default_letter\t%u\n", m_BestDefaultLetter);

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

	string strWildcard;
	ReadStringValue(f, "wildcard", strWildcard);
	if (strWildcard == "yes")
		m_Wildcard = true;
	else if (strWildcard == "no")
		m_Wildcard = false;
	else
		Die("Bad wildcard value '%s'", strWildcard.c_str());

	if (!m_IsInt)
		{
		m_MinValue = ReadFloatValue(f, "min");
		m_MedValue = ReadFloatValue(f, "med");
		m_MaxValue = ReadFloatValue(f, "max");
		m_UndefFreq = ReadFloatValue(f, "undef");
		}

	m_BestDefaultLetter = UINT_MAX;
	if (!m_Wildcard)
		m_BestDefaultLetter = ReadIntValue(f, "default_letter", UINT_MAX);

	LogOdds::FromTsv(f);

	m_BinTs.clear();
	if (!m_IsInt)
		{
		uint BinCount = (m_Wildcard ? m_AlphaSize - 1 : m_AlphaSize);
		for (uint i = 0; i + 1 < BinCount; ++i)
			m_BinTs.push_back(ReadFloatValue(f, "bint", i));

		string Line;
		bool Ok = ReadLineStdioFile(f, Line);
		asserta(!Ok);
		}

	CloseStdioFile(f);
	}
