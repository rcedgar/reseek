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

void FeatureTrainer::Init(const string &ChainsFN, const string &AlnsFN)
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

	m_AS = AS;
	m_LO.Init(AS);

	if (m_IsInt)
		return;

	asserta(!m_FloatValues.empty());
	DSS::Condense(m_FloatValues, m_AS, 
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

void FeatureTrainer::SetBackgroundFreqs()
	{
	const uint ChainCount = SIZE(m_Chains);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Background frequencies %s(%u)",
					 m_FeatureName, m_AS);

		const PDBChain &Chain = *m_Chains[ChainIndex];
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
				Letter = DSS::ValueToInt(m_BinTs, Value);
				}
			m_LO.AddBackgroundLetter(Letter);
			}
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
			m_LO.AddTruePair(LetterQ, LetterR);
			}
		if (!isgap(q))
			++QPos;
		if (!isgap(r))
			++RPos;
		}
	}

void FeatureTrainer::SetJointFreqs()
	{
	const uint SeqCount = m_Alns.GetSeqCount();
	asserta(SeqCount > 0);
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	uint LetterPairCount = 0;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Joint frequencies %s(%u)",
					 m_FeatureName, m_AS);
		SetJointFreqsPair(PairIndex);
		}
	}
void FeatureTrainer::Train()
	{
	SetBackgroundFreqs();
	SetJointFreqs();
	}

void FeatureTrainer::WriteSummary(FILE *f) const
	{
	if (f == 0)
		return;
	float ES = m_LO.GetExpectedScore();
	fprintf(f, "%s(%u)", m_FeatureName, m_AS);
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

void FeatureTrainer::WriteTsvHdr(FILE *f, uint ASCount) const
	{
	if (f == 0)
		return;

	fprintf(f, "feature=%s", m_FeatureName);
	fprintf(f, "\tNAS=%u", ASCount);
	fprintf(f, "\ttype=%s", m_IsInt ? "int" : "float");
	if(!m_IsInt)
		{
		fprintf(f, "\tmin=%.3g", m_MinValue);
		fprintf(f, "\tmed=%.3g", m_MedValue);
		fprintf(f, "\tmax=%.3g", m_MaxValue);
		fprintf(f, "\tundef=%.4f", m_UndefFreq);
		}
	fprintf(f, "\n");
	}

void FeatureTrainer::ToTsv(FILE *f) const
	{
	m_LO.WriteFeature(f, m_FeatureName);
	}
