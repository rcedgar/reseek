#include "myutils.h"
#include "featuretrainer.h"
#include "alpha.h"

UNDEF_BINNING StrToUB(const string &s)
	{
	if (s == "never") return UB_NeverUndefined;
	if (s == "onlyzero") return UB_UndefinedIsOnlyZero;
	if (s == "zeroov") return UB_UndefinedIsZeroOverload;
	if (s == "default") return UB_UndefinedIsDefaultLetter;
	if (s == "int") return UB_IntFeatureNoBinning;
	if (s == "ignore") return UB_IgnoreUndefined;
	if (s == "include") return UB_IncludeUndefined;
	Die("StrToUB(%s)", s.c_str());
	return UB_Invalid;
	}

const char *UBToStr(UNDEF_BINNING UB)
	{
	if (UB == UB_NeverUndefined) return "never";
	if (UB == UB_UndefinedIsOnlyZero) return "onlyzero";
	if (UB == UB_UndefinedIsZeroOverload) return "zeroov";
	if (UB == UB_UndefinedIsDefaultLetter) return "default";
	if (UB == UB_IntFeatureNoBinning) return "int";
	if (UB == UB_IgnoreUndefined) return "ignore";
	if (UB == UB_IncludeUndefined) return "include";
	asserta(false);
	return "*ERROR*";
	}

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
	m_Params.SetDSSParams(DM_DefaultFast);
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

void FeatureTrainer::SetOptionsFromCmdLine()
	{
	m_MinAQ = 0;
	m_MaxAQ = 1;
	m_MinPctId = 0;
	m_MaxPctId = 100;
	if (optset_minaq) m_MinAQ = (float) opt(minaq);
	if (optset_maxaq) m_MaxAQ = (float) opt(maxaq);
	if (optset_minpctid) m_MinPctId = (float) opt(minpctid);
	if (optset_maxpctid) m_MaxPctId = (float) opt(maxpctid);
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

void FeatureTrainer::SetUnalignedBackgroundChain(const PDBChain &Chain)
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
			double Value = m_D.GetFloatFeature(m_F, Pos);
			if (Value == DBL_MAX && m_UB == UB_IgnoreUndefined)
				continue;
			Letter = DSS::ValueToInt((float) Value, m_UB, m_AlphaSize,
									 m_BinTs, m_BestDefaultLetter);
			}
		AddUnalignedLetter(Letter);
		}
	}

void FeatureTrainer::SetUnalignedBackground()
	{
	const uint ChainCount = SIZE(m_Chains);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Unaligned background");
		const PDBChain &Chain = *m_Chains[ChainIdx];
		SetUnalignedBackgroundChain(Chain);
		}
	}

void FeatureTrainer::SetAlphaSize(uint AS, UNDEF_BINNING UB,
								  uint DefaultLetter)
	{
	m_AlphaSize = AS;
	m_UB = UB;
	Init(AS);

	if (m_IsInt)
		return;

	if (UB == UB_UndefinedIsDefaultLetter)
		m_BestDefaultLetter = DefaultLetter;

	asserta(!m_FloatValues.empty());
	DSS::Condense(m_FloatValues, m_AlphaSize,
				  m_UB, DefaultLetter, m_BestDefaultLetter,
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
	Params.SetDSSParams(DM_DefaultFast);

	DSS D;

	const uint ChainCount = SIZE(m_Chains);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		uint N = SIZE(m_FloatValues);
		double Pct = GetPct(UndefCount, N);
		ProgressStep(ChainIndex, ChainCount, "Values %s (%.2f%% undefined) %s",
					 m_FeatureName, Pct, UBToStr(m_UB));
		const PDBChain &Chain = *m_Chains[ChainIndex];
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			double Value = D.GetFloatFeature(m_F, Pos);
			if (Value == DBL_MAX)
				{
				if (m_UB == UB_IgnoreUndefined)
					continue;
				++UndefCount;
				}
			else
				{
				m_MinValue = min(m_MinValue, (float) Value);
				m_MaxValue = max(m_MaxValue, (float) Value);
				}
			m_FloatValues.push_back((float) Value);
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

bool FeatureTrainer::IncludePair(const string &QLabel) const
	{
	if (m_MinPctId > 0 || m_MaxPctId < 100)
		{
		float PctId = GetPctIdFromLabel(QLabel);
		return PctId >= m_MinPctId && PctId <= m_MaxPctId;
		}
	if (m_MinAQ > 0 || m_MaxAQ < 1)
		{
		float AQ = GetAQFromLabel(QLabel);
		return AQ >= m_MinAQ && AQ <= m_MaxAQ;
		}
	return true;
	}

void FeatureTrainer::UpdateJointCounts(uint PairIndex)
	{
	string QLabel = m_Alns.GetLabel(2*PairIndex);
	if (!IncludePair(QLabel))
		{
		++m_ExcludedPairCount;
		return;
		}

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
				m_Lock.lock();
				AddPair(LetterQ, LetterR);
				m_Lock.unlock();
				}
			else
				{
				float ValueQ = (float) DQ.GetFloatFeature(m_F, QPos);
				float ValueR = (float) DR.GetFloatFeature(m_F, RPos);
				bool IgnoreQ = (m_UB == UB_IgnoreUndefined && ValueQ == FLT_MAX);
				bool IgnoreR = (m_UB == UB_IgnoreUndefined && ValueR == FLT_MAX);

				if (!IgnoreQ && !IgnoreR)
					{
					LetterQ = DSS::ValueToInt(ValueQ, m_UB, m_AlphaSize,
											  m_BinTs, m_BestDefaultLetter);
					LetterR = DSS::ValueToInt(ValueR, m_UB, m_AlphaSize,
											  m_BinTs, m_BestDefaultLetter);
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

void FeatureTrainer::AddValue(float Value)
	{
	m_FloatValues.push_back(Value);
	}

void FeatureTrainer::Train()
	{
	const uint SeqCount = m_Alns.GetSeqCount();
	asserta(SeqCount > 0);
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	m_ExcludedPairCount = 0;
	m_Counter = 0;
//#pragma omp parallel for
	for (int PairIndex = 0; PairIndex < int(PairCount); ++PairIndex)
		{
		m_Lock.lock();
		uint Count = m_Counter++;
		m_Lock.unlock();
		ProgressStep(Count, PairCount, "Joint frequencies %s(%u) %s",
					 m_FeatureName, m_AlphaSize, UBToStr(m_UB));
		UpdateJointCounts(PairIndex);
		}
	if (m_UseUnalignedBackground)
		SetUnalignedBackground();
	if (m_UB != UB_UndefinedIsDefaultLetter)
		m_BestDefaultLetter = GetBestDefaultLetter(UINT_MAX);
	ProgressLog("%u excluded alns\n", m_ExcludedPairCount);
	}

void FeatureTrainer::WriteSummary(FILE *f) const
	{
	if (f == 0)
		return;
	float ES = GetExpectedScore();
	fprintf(f, "%s(%u) undef=%s", m_FeatureName, m_AlphaSize, UBToStr(m_UB));
	if (m_UB == UB_UndefinedIsDefaultLetter)
		fprintf(f, "(%u)", m_BestDefaultLetter);
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

void FeatureTrainer::BinTsToSrc(FILE *f) const
	{
	if (f == 0)
		return;
	if (m_IsInt)
		return;
	if (m_BinTs.empty())
		return;
	const uint N = SIZE(m_BinTs);
	fprintf(f, "\nuint DSS::ValueToInt_%s(double Value) const\n", m_FeatureName);
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

	fprintf(f, "\nstatic double %s_f_i[%u] = {\n",
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
	fprintf(f, "\nstatic double %s_S_i[%u][%u] = {\n",
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
	fprintf(f, "undef_binning\t%s\n", UBToStr(m_UB));
	if (!m_IsInt)
		{
		fprintf(f, "min\t%.3g\n", m_MinValue);
		fprintf(f, "med\t%.3g\n", m_MedValue);
		fprintf(f, "max\t%.3g\n", m_MaxValue);
		fprintf(f, "undef\t%.4f\n", m_UndefFreq);
		}
	if (m_UB == UB_UndefinedIsDefaultLetter)
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

		string Line;
		bool Ok = ReadLineStdioFile(f, Line);
		asserta(!Ok);
		}

	CloseStdioFile(f);
	}
