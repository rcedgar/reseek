#include "myutils.h"
#include "featuretrainer2.h"
#include "sfasta.h"
#include "alpha.h"
#include "valuetointtpl.h"
#include "round3sigfig.h"
#include "quarts.h"
#include "sort.h"

static float ScoreMx_Conf[16][16] = {
	{  1.2393, -6.3137, -7.4247, -5.3244, -3.4982, -5.4640, -4.4456, -0.8324, -2.1363, -2.8107, -5.2803,  0.0006, -2.0412, -4.2875, -4.2306, -3.6567},
	{ -6.3137,  1.5405,  0.8580,  0.3837, -1.4791,  0.4851, -1.8924, -5.2809, -4.1997, -3.1426,  0.5381, -3.4161, -5.4338, -1.1387, -2.5156, -4.3710},
	{ -7.4247,  0.8580,  1.9506, -1.2560, -2.8325, -0.2713, -3.5757, -7.4744, -6.2666, -5.0881, -1.1766, -5.6923, -6.7618, -3.0443, -4.7181, -7.1453},
	{ -5.3244,  0.3837, -1.2560,  2.0783,  0.6138,  0.6716, -1.0381, -2.2034, -1.6893, -1.6287, -0.0255, -0.9745, -2.2674,  0.7563,  0.6708, -2.3625},
	{ -3.4982, -1.4791, -2.8325,  0.6138,  2.3139,  0.6514, -1.2239, -3.2332,  1.1391, -2.5452, -0.9770, -0.0726, -0.7900, -0.4030,  0.3387, -1.9807},
	{ -5.4640,  0.4851, -0.2713,  0.6716,  0.6514,  2.3105, -1.5730, -4.3206, -1.8865, -3.1361, -0.5157, -1.7197, -3.4370, -0.9122, -0.9873, -3.4545},
	{ -4.4456, -1.8924, -3.5757, -1.0381, -1.2239, -1.5730,  2.5428, -2.1597,  0.4085,  0.3033,  1.3113, -2.7798, -3.0668,  0.0763, -0.4322, -0.9295},
	{ -0.8324, -5.2809, -7.4744, -2.2034, -3.2332, -4.3206, -2.1597,  2.7886, -2.0834, -0.7201, -2.9083, -0.8924, -1.6873, -0.8072, -0.1275, -0.1927},
	{ -2.1363, -4.1997, -6.2666, -1.6893,  1.1391, -1.8865,  0.4085, -2.0834,  2.4276, -0.4129, -1.2313, -0.1312,  0.7039, -1.2204,  0.2023, -0.0368},
	{ -2.8107, -3.1426, -5.0881, -1.6287, -2.5452, -3.1361,  0.3033, -0.7201, -0.4129,  2.6616, -1.3069, -2.2672, -2.9099,  1.1417, -1.0529,  1.6350},
	{ -5.2803,  0.5381, -1.1766, -0.0255, -0.9770, -0.5157,  1.3113, -2.9083, -1.2313, -1.3069,  2.3440, -2.8522, -4.2340,  0.0421, -1.1754, -2.2877},
	{  0.0006, -3.4161, -5.6923, -0.9745, -0.0726, -1.7197, -2.7798, -0.8924, -0.1312, -2.2672, -2.8522,  2.4597,  1.0263, -0.3818,  0.2272, -1.7486},
	{ -2.0412, -5.4338, -6.7618, -2.2674, -0.7900, -3.4370, -3.0668, -1.6873,  0.7039, -2.9099, -4.2340,  1.0263,  2.8322, -1.5363,  0.8756,  0.2578},
	{ -4.2875, -1.1387, -3.0443,  0.7563, -0.4030, -0.9122,  0.0763, -0.8072, -1.2204,  1.1417,  0.0421, -0.3818, -1.5363,  2.6540,  0.3376,  0.3106},
	{ -4.2306, -2.5156, -4.7181,  0.6708,  0.3387, -0.9873, -0.4322, -0.1275,  0.2023, -1.0529, -1.1754,  0.2272,  0.8756,  0.3376,  2.6673,  0.0630},
	{ -3.6567, -4.3710, -7.1453, -2.3625, -1.9807, -3.4545, -0.9295, -0.1927, -0.0368,  1.6350, -2.2877, -1.7486,  0.2578,  0.3106,  0.0630,  3.0810},
};

FEATURE FeatureTrainer2::m_F = FEATURE(-1);
uint FeatureTrainer2::m_AlphaSize = UINT_MAX;

void FeatureTrainer2::TruncLabel(string &Label)
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

void FeatureTrainer2::TruncLabel(const string &Label,
	string &TruncatedLabel)
	{
	TruncatedLabel = Label;
	size_t n = TruncatedLabel.find(' ');
	if (n != string::npos)
		TruncatedLabel.resize(n);
	n = TruncatedLabel.find('|');
	if (n != string::npos)
		TruncatedLabel.resize(n);
	n = TruncatedLabel.find('/');
	if (n != string::npos)
		TruncatedLabel.resize(n);
	}

void FeatureTrainer2::ReadChains(const string &FN,
	vector<PDBChain *> &Chains,
	map<string, uint> &LabelToChainIdx)
	{
	Chains.clear();
	LabelToChainIdx.clear();
	::ReadChains(FN, Chains);
	const uint N = SIZE(Chains);
	for (uint ChainIdx = 0; ChainIdx < N; ++ChainIdx)
		{
		string Label;
		TruncLabel(Chains[ChainIdx]->m_Label, Label);
		map<string, uint>::const_iterator iter =
			LabelToChainIdx.find(Label);
		if (iter != LabelToChainIdx.end())
			Die("Dupe chain label >%s", Label.c_str());
		LabelToChainIdx[Label] = ChainIdx;
		}
	}

void FeatureTrainer2::GetGapCounts(
	const string &Row1,
	const string &Row2,
	uint &ColCount,
	uint &Opens,
	uint &Exts)
	{
	const uint RowLength = SIZE(Row1);
	asserta(SIZE(Row2) == RowLength);
	bool InGap = false;
	ColCount = 0;
	Opens = 0;
	Exts = 0;
	for (uint Col = 0; Col < RowLength; ++Col)
		{
		char c1 = Row1[Col];
		char c2 = Row2[Col];
		if (islower(c1) || islower(c2))
			continue;
		bool gap1 = isgap(c1);
		bool gap2 = isgap(c2);
		asserta(!(gap1 && gap2));
		++ColCount;
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

void FeatureTrainer2::GetGapCountVecs(
	const vector<string> &Rows,
	vector<uint> &ColCountVec,
	vector<uint> &OpenVec,
	vector<uint> &ExtVec)
	{
	const uint RowCount = SIZE(Rows);
	asserta(RowCount%2 == 0);
	const uint AlnCount = RowCount/2;
	ColCountVec.clear();
	OpenVec.clear();
	ExtVec.clear();
	ColCountVec.reserve(AlnCount);
	OpenVec.reserve(AlnCount);
	ExtVec.reserve(AlnCount);
	for (uint AlnIdx = 0; AlnIdx < AlnCount; ++AlnIdx)
		{
		const string &Row1 = Rows[2*AlnIdx];
		const string &Row2 = Rows[2*AlnIdx+1];
		uint ColCount, Opens, Exts;
		GetGapCounts(Row1, Row2, ColCount, Opens, Exts);
		ColCountVec.push_back(ColCount);
		OpenVec.push_back(Opens);
		ExtVec.push_back(Exts);
		}
	asserta(SIZE(ColCountVec) == AlnCount);
	asserta(SIZE(OpenVec) == AlnCount);
	asserta(SIZE(ExtVec) == AlnCount);
	}

void FeatureTrainer2::AppendAlns(
	const string &FN,
	const map<string, uint> &LabelToChainIdx,
	bool AlnsAreTPs,
	vector<string> &Rows,
	vector<string> &Labels,
	vector<uint> &ChainIdxs,
	vector<bool> &TPs)
	{
	uint N = SIZE(Rows);
	asserta(N%2 == 0);
	asserta(SIZE(Labels) == N);
	asserta(SIZE(ChainIdxs) == N);
	asserta(SIZE(TPs) == N/2);

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
		map<string, uint>::const_iterator iter = LabelToChainIdx.find(Label);
		if (iter == LabelToChainIdx.end())
			Die("Label not found >%s", Label.c_str());
		uint ChainIdx = iter->second;
		ChainIdxs.push_back(ChainIdx);
		const unsigned L = SF.GetSeqLength();
		asserta(L != 0);
		string Row;
		Row.reserve(L);
		for (uint i = 0; i < L; ++i)
			Row.push_back(Seq[i]);
		if (!Row1)
			{
			asserta(SIZE(Row) == SIZE(Rows.back()));
			TPs.push_back(AlnsAreTPs);
			}
		TruncLabel(Label, Label);
		Labels.push_back(Label);
		Rows.push_back(Row);
		Row1 = !Row1;
		}

	N = SIZE(Rows);
	asserta(N%2 == 0);
	asserta(SIZE(Labels) == N);
	asserta(SIZE(ChainIdxs) == N);
	asserta(SIZE(TPs) == N/2);
	}

void FeatureTrainer2::SetFloatFeature(FEATURE F, uint AlphaSize)
	{
	asserta(!FeatureIsInt(F));
	m_F = F;
	m_AlphaSize = AlphaSize;
	}

void FeatureTrainer2::SetIntFeature(FEATURE F)
	{
	asserta(FeatureIsInt(F));
	m_F = F;
	m_AlphaSize = DSS::GetAlphaSize(F);
	}

void FeatureTrainer2::GetFloatSeqs(
	const vector<PDBChain *> &Chains,
	DSS &D,
	float ReplaceUndefValue,
	vector<vector<float> > &Seqs)
	{
	Seqs.clear();
	const uint ChainCount = SIZE(Chains);
	Seqs.resize(ChainCount);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		vector<float> &ChainFloatSeq = Seqs[ChainIndex];
		const PDBChain &Chain = *Chains[ChainIndex];
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		ChainFloatSeq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = D.GetFloatFeature(m_F, Pos);
			if (Value == FLT_MAX)
				Value = ReplaceUndefValue;
			ChainFloatSeq.push_back(Value);
			}
		}
	}

void FeatureTrainer2::GetIntSeqs(
	const vector<PDBChain *> &Chains,
	DSS &D,
	uint ReplaceUndefValue,
	const vector<float> &BinTs, vector<vector<uint> > &Seqs,
	uint &UndefCount)
	{
	Seqs.clear();
	UndefCount = 0;
	asserta(SIZE(BinTs) == m_AlphaSize-1);
	const uint ChainCount = SIZE(Chains);
	Seqs.resize(ChainCount);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		vector<uint> &Seq = Seqs[ChainIndex];
		const PDBChain &Chain = *Chains[ChainIndex];
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		Seq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = D.GetFloatFeature(m_F, Pos);
			if (Value == FLT_MAX)
				{
				++UndefCount;
				Seq.push_back(ReplaceUndefValue);
				}
			else
				{
				uint Letter = ValueToIntTpl<false>
					(Value, m_AlphaSize, BinTs, UINT_MAX);
				Seq.push_back(Letter);
				}
			}
		}
	}

void FeatureTrainer2::GetSortedFloatValues(
	const vector<PDBChain *> &Chains,
	DSS &D,
	bool IgnoreUndef,
	float ReplaceUndefValue,
	vector<float> &Values,
	uint &UndefCount,
	uint &ReplaceCount)
	{
	UndefCount = 0;
	ReplaceCount = 0;
	const uint ChainCount = SIZE(Chains);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = D.GetFloatFeature(m_F, Pos);
			if (Value == FLT_MAX)
				{
				++UndefCount;
				if (IgnoreUndef)
					continue;
				Value = ReplaceUndefValue;
				++ReplaceCount;
				}
			Values.push_back(Value);
			}
		}
	sort(Values.begin(), Values.end());
	}

void FeatureTrainer2::LogSortedValueStats(const string &Msg,
	const vector<float> &Values)
	{
	const uint N = SIZE(Values);
	asserta(N > 0);
	uint M = 0;
	float Max = 0;
	for (uint i = 0; i < N; ++i)
		{
		if (Values[i] == FLT_MAX)
			++M;
		else
			Max = max(Values[i], Max);
		}

	Log("LogSortedValueStats(%s) N=%u, FLT_MAX=%u", Msg.c_str(), N, M);
	Log(", min %.3g, med %.3g, max %.3g", Values[0], Values[N/2], Max);
	Log("\n");
	}

void FeatureTrainer2::LogChainFloatSeqsStats(const string &Msg,
	const vector<vector<float> > &Seqs)
	{
	const uint SeqCount = SIZE(Seqs);
	asserta(SeqCount > 0);
	uint M = 0;
	uint N = 0;
	float Min = 0;
	float Max = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const vector<float> &Seq = Seqs[i];
		const uint L = SIZE(Seq);
		N += L;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			if (Seq[i] == FLT_MAX)
				++M;
			else if (i == 0)
				{
				Max = Seq[i];
				Min = Seq[i];
				}
			else
				{
				Max = max(Seq[i], Max);
				Min = min(Seq[i], Min);
				}
			}
		}
	Log("LogChainFloatSeqStats(%s) N=%u, FLT_MAX=%u", Msg.c_str(), N, M);
	Log(", min %.3g, max %.3g\n", Min, Max);
	}

void FeatureTrainer2::LogChainIntSeqsStats(const string &Msg,
	const vector<vector<uint> > &Seqs)
	{
	const uint SeqCount = SIZE(Seqs);
	asserta(SeqCount > 0);
	uint N = 0;
	uint M = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const vector<uint> &Seq = Seqs[i];
		const uint L = SIZE(Seq);
		N += L;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			if (Seq[i] == FLT_MAX)
				++M;
			}
		}
	Log("LogChainIntSeqStats(%s) N=%u, UINT_MAX=%u", Msg.c_str(), N, M);
	}

void FeatureTrainer2::BinTsToTsv(
	FILE *f,
	const vector<float> &Ts)
	{
	if (f == 0)
		return;
	asserta(!FeatureIsInt(m_F));
	asserta(SIZE(Ts) + 1 == m_AlphaSize);
	for (uint i = 0 ; i + 1 < m_AlphaSize; ++i)
		fprintf(f, "%u\t%.4g\n", i, Ts[i]);
	}

void FeatureTrainer2::FreqsToSrc(
	FILE *f,
	const vector<float> &Freqs)
	{
	if (f == 0)
		return;
	asserta(SIZE(Freqs) == m_AlphaSize);

	fprintf(f, "\nstatic float %s_f_i[%u] = {\n",
			FeatureToStr(m_F), m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		fprintf(f, "	%.4f,\n", Freqs[i]);
	fprintf(f, "};\n");
	}

void FeatureTrainer2::ScoreMxToSrc(
	FILE *f,
	const vector<vector<float> > &ScoreMx)
	{
	if (f == 0)
		return;
	asserta(SIZE(ScoreMx) == m_AlphaSize);
	fprintf(f, "\nstatic float %s_S_ij[%u][%u] = {\n",
			FeatureToStr(m_F), m_AlphaSize, m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		fprintf(f, "	{");
		for (uint j = 0; j < m_AlphaSize; ++j)
			fprintf(f, " %11.4g,", ScoreMx[i][j]);
		fprintf(f, " }, // %u\n", i);
		}
	fprintf(f, "};\n");
	}

void FeatureTrainer2::GetAlignedLetterCounts(
	const vector<vector<uint> > &ChainIntSeqs,
	const vector<string> &Rows,
	const vector<uint> &RowChainIdxs,
	bool IgnoreUndef,
	vector<uint> &Counts)
	{
	Counts.clear();
	Counts.resize(m_AlphaSize);

	const uint RowCount = SIZE(Rows);
	asserta(SIZE(RowChainIdxs) == RowCount);
	asserta(RowCount > 0);
	asserta(RowCount%2 == 0);
	const uint AlnCount = RowCount/2;
	for (uint AlnIdx = 0; AlnIdx < AlnCount; ++AlnIdx)
		{
		uint RowIdx1 = AlnIdx*2;
		uint RowIdx2 = RowIdx1 + 1;

		const string &Row1 = Rows[RowIdx1];
		const string &Row2 = Rows[RowIdx2];

		const uint ChainIdx1 = RowChainIdxs[RowIdx1];
		const uint ChainIdx2 = RowChainIdxs[RowIdx2];

		const vector<uint> &IntSeq1 = ChainIntSeqs[ChainIdx1];
		const vector<uint> &IntSeq2 = ChainIntSeqs[ChainIdx2];
		uint L1 = SIZE(IntSeq1);
		uint L2 = SIZE(IntSeq2);

		const uint ColCount = SIZE(Row1);
		asserta(SIZE(Row2) == ColCount);

		uint Pos1 = 0;
		uint Pos2 = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char c1 = Row1[Col];
			char c2 = Row2[Col];

			if (isupper(c1) && isupper(c2))
				{
				uint Letter1 = IntSeq1[Pos1];
				uint Letter2 = IntSeq2[Pos2];
				if (Letter1 == UINT_MAX || Letter2 == UINT_MAX)
					asserta(IgnoreUndef);
				else
					{
					asserta(Letter1 < m_AlphaSize && Letter2 < m_AlphaSize);
					Counts[Letter1] += 1;
					Counts[Letter2] += 1;
					}
				}

			if (!isgap(c1))
				++Pos1;
			if (!isgap(c2))
				++Pos2;
			}
		asserta(Pos1 == L1);
		asserta(Pos2 == L2);
		}
	}

void FeatureTrainer2::GetLetterCounts(
	const vector<vector<uint> > &ChainIntSeqs,
	const vector<uint> &ChainIdxs,
	bool IgnoreUndef,
	vector<uint> &Counts)
	{
	Counts.clear();
	Counts.resize(m_AlphaSize);

	set<uint> ChainIdxSet;
	for (uint k = 0; k < SIZE(ChainIdxs); ++k)
		ChainIdxSet.insert(ChainIdxs[k]);

	for (set<uint>::const_iterator iter = ChainIdxSet.begin();
		iter != ChainIdxSet.end(); ++iter)
		{
		uint ChainIdx = *iter;
		const vector<uint> &IntSeq = ChainIntSeqs[ChainIdx];
		uint L = SIZE(IntSeq);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint Letter = IntSeq[Pos];
			if (Letter == UINT_MAX)
				asserta(IgnoreUndef);
			else
				{
				asserta(Letter < m_AlphaSize);
				Counts[Letter] += 1;
				}
			}
		}
	}

void FeatureTrainer2::GetAlignedLetterPairCounts(
	const vector<vector<uint> > &ChainIntSeqs,
	const vector<string> &Rows,
	const vector<uint> &RowChainIdxs,
	bool IgnoreUndef,
	vector<vector<uint> > &CountMx)
	{
	CountMx.clear();
	CountMx.resize(m_AlphaSize);
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		CountMx[Letter].resize(m_AlphaSize);

	const uint RowCount = SIZE(Rows);
	asserta(SIZE(RowChainIdxs) == RowCount);
	asserta(RowCount > 0);
	asserta(RowCount%2 == 0);
	const uint AlnCount = RowCount/2;
	for (uint AlnIdx = 0; AlnIdx < AlnCount; ++AlnIdx)
		{
		uint RowIdx1 = AlnIdx*2;
		uint RowIdx2 = RowIdx1 + 1;

		const string &Row1 = Rows[RowIdx1];
		const string &Row2 = Rows[RowIdx2];

		const uint ChainIdx1 = RowChainIdxs[RowIdx1];
		const uint ChainIdx2 = RowChainIdxs[RowIdx2];

		const vector<uint> &IntSeq1 = ChainIntSeqs[ChainIdx1];
		const vector<uint> &IntSeq2 = ChainIntSeqs[ChainIdx2];
		uint L1 = SIZE(IntSeq1);
		uint L2 = SIZE(IntSeq2);

		const uint ColCount = SIZE(Row1);
		asserta(SIZE(Row2) == ColCount);

		uint Pos1 = 0;
		uint Pos2 = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char c1 = Row1[Col];
			char c2 = Row2[Col];

			if (isupper(c1) && isupper(c2))
				{
				uint Letter1 = IntSeq1[Pos1];
				uint Letter2 = IntSeq2[Pos2];
				if (Letter1 == UINT_MAX || Letter2 == UINT_MAX)
					asserta(IgnoreUndef);
				else
					{
					asserta(Letter1 < m_AlphaSize && Letter2 < m_AlphaSize);

					// Diagonal does not need special case here
					// Off-diagonals are double-counted by assuming symmetry
					// On-diagonals are double-counted by not checking for diagonal
					// => all counts are doubled, frequencies are correct
					CountMx[Letter1][Letter2] += 1;
					CountMx[Letter2][Letter1] += 1;
					}
				}

			if (!isgap(c1))
				++Pos1;
			if (!isgap(c2))
				++Pos2;
			}
		asserta(Pos1 == L1);
		asserta(Pos2 == L2);
		//Log("Joint counts\n");//@@
		//for (uint i = 0; i < m_AlphaSize; ++i)
		//	{
		//	const vector<uint> &Row = CountMx[i];
		//	Log("[%2u] ", i);
		//	for (uint j = 0; j < m_AlphaSize; ++j)
		//		Log("%10u", Row[j]);
		//	Log("\n");
		//	}
		//Die("");//@@
		}

	Log("Joint counts\n");
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		const vector<uint> &Row = CountMx[i];
		Log("[%2u] ", i);
		for (uint j = 0; j < m_AlphaSize; ++j)
			Log("%10u", Row[j]);
		Log("\n");
		}
	}

void FeatureTrainer2::GetIntSeqLetterCounts(
	const vector<vector<uint> > &Seqs,
	vector<uint> &Counts,
	uint &UndefCount)
	{
	const uint SeqCount = SIZE(Seqs);
	asserta(SeqCount > 0);
	Counts.clear();
	Counts.resize(m_AlphaSize);
	UndefCount = 0;
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		const vector<uint> &Seq = Seqs[SeqIdx];
		const uint L = SIZE(Seq);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint Letter = Seq[Pos];
			if (Letter == UINT_MAX)
				++UndefCount;
			else if (Letter < m_AlphaSize)
				Counts[Letter] += 1;
			else
				asserta(false);
			}
		}
	}

float FeatureTrainer2::GetRelativeEntropy(
	const vector<vector<float> > &FreqMx,
	const vector<vector<float> > &ScoreMx)
	{
	asserta(SIZE(FreqMx) == m_AlphaSize);
	asserta(SIZE(ScoreMx) == m_AlphaSize);
	float Hrel = 0;
	float Sumf = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			float f = FreqMx[Letter1][Letter2];
			Hrel += f*ScoreMx[Letter1][Letter2];
			Sumf += f;
			}
		}
	asserta(feq(Sumf, 1));
	return Hrel;
	}

float FeatureTrainer2::GetShannonEntropy(
	const vector<vector<float> > &FreqMx)
	{
	asserta(SIZE(FreqMx) == m_AlphaSize);
	float H = 0;
	float Sumf = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			float f = FreqMx[Letter1][Letter2];
			if (f > 0)
				H -= f*log(f);
			Sumf += f;
			}
		}
	asserta(feq(Sumf, 1));
	return H;
	}

float FeatureTrainer2::GetExpectedScore(
	vector<vector<float> > &ScoreMx,
	const vector<float> &Freqs)
	{
	asserta(SIZE(ScoreMx) == m_AlphaSize);
	asserta(SIZE(Freqs) == m_AlphaSize);
	float ES = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		float f1 = Freqs[Letter1];
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			float f2 = Freqs[Letter1];
			ES += f1*f2*ScoreMx[Letter1][Letter2];
			}
		}
	return ES;
	}

void FeatureTrainer2::GetLogOddsMx(
	const vector<float> &Freqs,
	const vector<vector<float> > &FreqMx,
	vector<vector<float> > &ScoreMx)
	{
	ScoreMx.clear();
	ScoreMx.resize(m_AlphaSize);
	
	float SumLetterFreqs = 0;
	float SumPairFreqs = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		ScoreMx[Letter1].resize(m_AlphaSize);
		float f1 = Freqs[Letter1];
		SumLetterFreqs += f1;
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			float f2 = Freqs[Letter2];
			float PairFreq = FreqMx[Letter1][Letter2];
			float ExpectedFreq = f1*f2;
			float Score = FLT_MAX;
			if (PairFreq == 0 || ExpectedFreq == 0)
				Score = 0;
			else
				{
				float Ratio = PairFreq/ExpectedFreq;
				Score = log(Ratio);
				}
			ScoreMx[Letter1][Letter2] = Score;
			SumPairFreqs += PairFreq;
			}
		}
	asserta(feq(SumLetterFreqs, 1));
	asserta(feq(SumPairFreqs, 1));
	}

void FeatureTrainer2::GetFreqs(
	const vector<uint> &Counts,
	vector<float> &Freqs)
	{
	asserta(SIZE(Counts) == m_AlphaSize);
	uint Sum = 0;
	for (uint i = 0; i < m_AlphaSize; ++i)
		Sum += Counts[i];
	const float fSum = float(Sum);

	Freqs.clear();
	float SumFreqs = 0;
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		float Freq = Counts[i]/fSum;
		Freqs.push_back(Freq);
		SumFreqs += Freq;
		}
	asserta(feq(SumFreqs, 1));
	}

void FeatureTrainer2::GetFreqMx(
	const vector<vector<uint> > &CountMx,
	vector<vector<float> > &FreqMx)
	{
	asserta(SIZE(CountMx) == m_AlphaSize);
	FreqMx.clear();
	FreqMx.resize(m_AlphaSize);
	uint N = 0;
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		FreqMx[i].resize(m_AlphaSize);
		const vector<uint> &Row = CountMx[i];
		asserta(SIZE(Row) == m_AlphaSize);
		for (uint j = 0; j < m_AlphaSize; ++j)
			{
			uint n = CountMx[i][j];
			N += n;
			}
		}
	asserta(N > 0);
	const float fN = float(N);
	float SumFreq = 0;
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		const vector<uint> &Row = CountMx[i];
		asserta(SIZE(Row) == m_AlphaSize);
		for (uint j = 0; j < m_AlphaSize; ++j)
			{
			uint n = CountMx[i][j];
			float Freq = n/fN;
			FreqMx[i][j] = Freq;
			SumFreq += Freq;
			}
		}
	asserta(feq(SumFreq, 1));
	}

void FeatureTrainer2::GetAlnSubstScores(
	const vector<vector<uint> > &ChainIntSeqs,
	const vector<string> &Rows,
	const vector<uint> &RowChainIdxs,
	bool IgnoreUndef,
	uint ReplaceUndefByThisLetter,
	const vector<vector<float> > &ScoreMx,
	vector<float > &SubstScores)
	{
	SubstScores.clear();
	asserta(SIZE(ScoreMx) == m_AlphaSize);
	asserta(SIZE(ScoreMx[0]) == m_AlphaSize);
	const uint RowCount = SIZE(Rows);
	asserta(SIZE(RowChainIdxs) == RowCount);
	asserta(RowCount > 0);
	asserta(RowCount%2 == 0);
	const uint AlnCount = RowCount/2;
	SubstScores.reserve(AlnCount);
	for (uint AlnIdx = 0; AlnIdx < AlnCount; ++AlnIdx)
		{
		uint RowIdx1 = AlnIdx*2;
		uint RowIdx2 = RowIdx1 + 1;

		const string &Row1 = Rows[RowIdx1];
		const string &Row2 = Rows[RowIdx2];

		const uint ChainIdx1 = RowChainIdxs[RowIdx1];
		const uint ChainIdx2 = RowChainIdxs[RowIdx2];

		const vector<uint> &IntSeq1 = ChainIntSeqs[ChainIdx1];
		const vector<uint> &IntSeq2 = ChainIntSeqs[ChainIdx2];
		uint L1 = SIZE(IntSeq1);
		uint L2 = SIZE(IntSeq2);

		const uint ColCount = SIZE(Row1);
		asserta(SIZE(Row2) == ColCount);

		uint Pos1 = 0;
		uint Pos2 = 0;
		float SubstScore = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char c1 = Row1[Col];
			char c2 = Row2[Col];

			if (isupper(c1) && isupper(c2))
				{
				uint Letter1 = IntSeq1[Pos1];
				uint Letter2 = IntSeq2[Pos2];
				if (Letter1 == UINT_MAX || Letter2 == UINT_MAX)
					asserta(IgnoreUndef);
				else
					{
					if (Letter1 == UINT_MAX)
						Letter1 = ReplaceUndefByThisLetter;
					if (Letter2 == UINT_MAX)
						Letter2 = ReplaceUndefByThisLetter;
					asserta(Letter1 < m_AlphaSize && Letter2 < m_AlphaSize);
					SubstScore += ScoreMx[Letter1][Letter2];
					}
				}

			if (!isgap(c1))
				++Pos1;
			if (!isgap(c2))
				++Pos2;
			}
		asserta(Pos1 == L1);
		asserta(Pos2 == L2);
		SubstScores.push_back(SubstScore);
		if (AlnIdx < 5)//@@
			Log("GetAlnSubstScore(%u) = %.3g\n", AlnIdx, SubstScore);//@@
		}
	asserta(SIZE(SubstScores) == AlnCount);
	}

void FeatureTrainer2::GetSteps(
	const vector<float> &AlnScores,
	const vector<bool> &TPs,
	vector<float> &StepScores,
	vector<float> &StepTPfs,
	vector<float> &StepFPfs)
	{
	const uint AlnCount = SIZE(AlnScores);
	asserta(SIZE(TPs) == AlnCount);
	asserta(AlnCount > 100);

	StepScores.clear();
	StepTPfs.clear();
	StepFPfs.clear();

	StepScores.reserve(AlnCount);
	StepTPfs.reserve(AlnCount);
	StepFPfs.reserve(AlnCount);

	uint TPCount = 0;
	uint FPCount = 0;
	for (uint k = 0; k < AlnCount; ++k)
		if (TPs[k])
			++TPCount;
		else
			++FPCount;
	asserta(TPCount > 100);
	asserta(FPCount > 100);

	uint *Order = myalloc(uint, AlnCount);
	QuickSortOrderDesc(AlnScores.data(), AlnCount, Order);

	uint NTP = 0;
	uint NFP = 0;
	float CurrentScore = AlnScores[Order[0]];
	for (uint k = 0; k < AlnCount; ++k)
		{
		uint i = Order[k];
		float Score = AlnScores[i];
		bool TP = TPs[i];
		if (Score != CurrentScore)
			{
			asserta(Score < CurrentScore);
			float TPf = float(NTP)/TPCount;
			float FPf = float(NFP)/FPCount;
			StepScores.push_back(CurrentScore);
			StepTPfs.push_back(TPf);
			StepFPfs.push_back(FPf);
			CurrentScore = Score;
			}
		if (TP)
			++NTP;
		else
			++NFP;
		}
	myfree(Order);
	asserta(NTP == TPCount);
	asserta(NFP == FPCount);
	asserta(NTP > 0);
	asserta(NFP > 0);
	}

float FeatureTrainer2::CalcArea(
	const vector<float> &AlnScores,
	const vector<bool> &TPs)
	{
	const uint AlnCount = SIZE(AlnScores);
	asserta(SIZE(TPs) == AlnCount);
	asserta(AlnCount > 100);

	uint TPCount = 0;
	uint FPCount = 0;
	for (uint k = 0; k < AlnCount; ++k)
		if (TPs[k])
			++TPCount;
		else
			++FPCount;
	asserta(TPCount > 100);
	asserta(FPCount > 100);

	uint *Order = myalloc(uint, AlnCount);
	QuickSortOrderDesc(AlnScores.data(), AlnCount, Order);

	uint NTP = 0;
	uint NFP = 0;
	float PrevTPf = 0;
	float PrevFPf = 0;
	float CurrentScore = AlnScores[Order[0]];
	float Area = 0;
	for (uint k = 0; k < AlnCount; ++k)
		{
		uint i = Order[k];
		float Score = AlnScores[i];
		bool TP = TPs[i];
		if (Score != CurrentScore)
			{
			float TPf = float(NTP)/TPCount;
			float FPf = float(NFP)/FPCount;
			Area += TPf*(FPf - PrevFPf)/2;
			}
		if (TP)
			++NTP;
		else
			++NFP;
		}
	asserta(NTP == TPCount);
	asserta(NFP == FPCount);
	asserta(NTP > 0);
	asserta(NFP > 0);

	float TPf = float(NTP)/TPCount;
	float FPf = float(NFP)/FPCount;
	Area += TPf*(FPf - PrevFPf)/2;
	return Area;
	}

void FeatureTrainer2::GetChainIntSeqs_Int(
	const vector<PDBChain *> &Chains,
	vector<vector<uint> > &IntSeqs,
	uint &LetterCount,
	uint &UndefCount)
	{
	asserta(FeatureIsInt(m_F));
	const uint ChainCount = SIZE(Chains);
	IntSeqs.clear();
	IntSeqs.resize(ChainCount);
	LetterCount = 0;
	UndefCount = 0;
	DSS D;
	opt_force_undef = true;
	optset_force_undef = true;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *Chains[ChainIdx];
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		LetterCount += L;
		vector<uint> &IntSeq = IntSeqs[ChainIdx];
		IntSeq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint Letter = D.GetFeature(m_F, Pos);
			if (Letter == UINT_MAX)
				++UndefCount;
			IntSeq.push_back(Letter);
			}
		}
	opt_force_undef = false;
	optset_force_undef = false;
	}

void FeatureTrainer2::GetChainIntSeqs_Float(
	const vector<PDBChain *> &Chains,
	vector<vector<uint> > &IntSeqs,
	const vector<float> &BinTs,
	uint &LetterCount,
	uint &UndefCount)
	{
	asserta(FeatureIsInt(m_F));
	const uint ChainCount = SIZE(Chains);
	IntSeqs.clear();
	IntSeqs.resize(ChainCount);
	UndefCount = 0;
	LetterCount = 0;
	DSS D;
	opt_force_undef = true;
	optset_force_undef = true;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *Chains[ChainIdx];
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		vector<uint> &IntSeq = IntSeqs[ChainIdx];
		IntSeq.reserve(L);
		LetterCount += L;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = D.GetFloatFeature(m_F, Pos);
			if (Value == FLT_MAX)
				{
				++UndefCount;
				IntSeq.push_back(UINT_MAX);
				}
			else
				{
				uint Letter = ValueToIntTpl<false>(
					Value, m_AlphaSize, BinTs, UINT_MAX);
				IntSeq.push_back(Letter);
				}
			}
		}
	opt_force_undef = false;
	optset_force_undef = false;
	}

void FeatureTrainer2::LogFreqs(
	const vector<float> &Freqs)
	{
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		Log("[%2u] %c", Letter, g_LetterToCharAmino[Letter]);
		float Freq = Freqs[Letter];
		Log("  %6.3f", Freq);
		Log("\n");
		}
	}

void FeatureTrainer2::GetFreqDiffs(
	const vector<float> &Freqs1,
	const vector<float> &Freqs2,
	vector<float> &Diffs)
	{
	Diffs.clear();
	asserta(SIZE(Freqs1) == m_AlphaSize);
	asserta(SIZE(Freqs2) == m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		Diffs.push_back(Freqs1[i] - Freqs2[i]);
	}

void FeatureTrainer2::GetAlnScores(
	const vector<float> &AlnSubstScores,
	const vector<uint> &AlnColCountVec,
	const vector<uint> &AlnOpenVec,
	const vector<uint> &AlnExtVec,
	float OpenPenalty,
	float ExtPenalty,
	float Bias,
	vector<float> &AlnScores)
	{
	const uint AlnCount = SIZE(AlnSubstScores);
	asserta(SIZE(AlnOpenVec) == AlnCount);
	asserta(SIZE(AlnExtVec) == AlnCount);

	AlnScores.clear();
	AlnScores.reserve(AlnCount);
	for (uint AlnIndex = 0; AlnIndex < AlnCount; ++AlnIndex)
		{
		uint ColCount = AlnColCountVec[AlnIndex];
		uint OpenCount = AlnOpenVec[AlnIndex];
		uint ExtCount = AlnExtVec[AlnIndex];
		float SubstScore = AlnSubstScores[AlnIndex];
		float AlnScore = SubstScore;
		AlnScore -= OpenCount*OpenPenalty;
		AlnScore -= ExtCount*ExtPenalty;
		AlnScores.push_back(AlnScore);
		}
	}

void FeatureTrainer2::Round3SigFig(
	const vector<float> &InputValues,
	vector<float> &RoundedValues)
	{
	const uint N = SIZE(InputValues);
	RoundedValues.clear();
	RoundedValues.reserve(N);
	for (uint i = 0; i < N; ++i)
		RoundedValues.push_back(round3sigfig(InputValues[i]));
	}

void FeatureTrainer2::LogAlnScoreQuarts(
	const vector<float> &AlnScores,
	const vector<bool> &TPs)
	{
	const uint AlnCount = SIZE(AlnScores);
	asserta(SIZE(TPs) == AlnCount);
	vector<float> TPScores;
	vector<float> FPScores;
	TPScores.reserve(AlnCount);
	FPScores.reserve(AlnCount);
	for (uint i = 0; i < AlnCount; ++i)
		{
		float Score = AlnScores[i];
		if (TPs[i])
			TPScores.push_back(Score);
		else
			FPScores.push_back(Score);
		}
	QuartsFloat QF;
	GetQuartsFloat(TPScores, QF);
	Log("TP scores: ");
	QF.LogMe();

	GetQuartsFloat(FPScores, QF);
	Log("FP scores: ");
	QF.LogMe();
	}

void FeatureTrainer2::TrainIntFeatureNoUndefs(
	FEATURE F,
	const string &ChainFN,
	const string &TrainTPAlnFN,
	const string &TrainFPAlnFN,
	const string &EvalTPAlnFN,
	const string &EvalFPAlnFN,
	vector<vector<float > > &ScoreMx)
	{
	asserta(false);
	}

void FeatureTrainer2::TrainIntFeatureIgnoreUndefs(
	FEATURE F,
	const string &ChainFN,
	const string &TrainTPAlnFN,
	const string &TrainFPAlnFN,
	const string &EvalTPAlnFN,
	const string &EvalFPAlnFN,
	vector<vector<float > > &ScoreMx)
	{
	const bool IgnoreUndef = true;

	vector<PDBChain *> Chains;
	map<string, uint> LabelToChainIdx;
	ReadChains(ChainFN, Chains, LabelToChainIdx);

	vector<vector<uint> > ChainIntSeqs;
	uint UndefCount = 0;
	uint LetterCount = 0;
	GetChainIntSeqs_Int(Chains, ChainIntSeqs, LetterCount, UndefCount);
	Log("%u / %u (%.2f%%) undefined letters in chains\n",
		UndefCount, LetterCount, GetPct(UndefCount, LetterCount));

	vector<bool> TrainTPs;
	vector<string> TPTrainRows;
	vector<string> TPTrainLabels;
	vector<uint> TPTrainChainIdxs;
	AppendAlns(TrainTPAlnFN, LabelToChainIdx, true,
	  TPTrainRows, TPTrainLabels, TPTrainChainIdxs, TrainTPs);

	vector<string> AllTrainRows = TPTrainRows;
	vector<string> AllTrainLabels = TPTrainLabels;
	vector<uint> AllTrainChainIdxs = TPTrainChainIdxs;

	AppendAlns(TrainFPAlnFN, LabelToChainIdx, false,
	  AllTrainRows, AllTrainLabels, AllTrainChainIdxs, TrainTPs);

	vector<uint> TrainTPLetterCounts;
	vector<float> TrainTPLetterFreqs;
	GetAlignedLetterCounts(ChainIntSeqs, TPTrainRows,
	  TPTrainChainIdxs, IgnoreUndef, TrainTPLetterCounts);
	GetLetterCounts(ChainIntSeqs, AllTrainChainIdxs,
		IgnoreUndef, TrainTPLetterCounts);
	GetFreqs(TrainTPLetterCounts, TrainTPLetterFreqs);

	vector<vector<uint> > TrainAlnLetterPairCountMx;
	vector<vector<float> > TrainAlnLetterPairFreqMx;
	GetAlignedLetterPairCounts(ChainIntSeqs, TPTrainRows,
	  TPTrainChainIdxs, IgnoreUndef, TrainAlnLetterPairCountMx);
	GetFreqMx(TrainAlnLetterPairCountMx, TrainAlnLetterPairFreqMx);

	vector<vector<float> > ScoreMxAlnBg;
	GetLogOddsMx(TrainTPLetterFreqs, TrainAlnLetterPairFreqMx, ScoreMxAlnBg);

	
	{//@@
	for (uint i = 0; i < 16; ++i)
		for (uint j = 0; j < 16; ++j)
			ScoreMxAlnBg[i][j] = ScoreMx_Conf[i][j];
	}//@@

	LogFreqs(TrainTPLetterFreqs);

	Log("\n");
	Log("// ScoreMxAlnBg\n");
	ScoreMxToSrc(g_fLog, ScoreMxAlnBg);
	//Log("\n");
	//Log("// ScoreMxAllBg\n");
	//ScoreMxToSrc(g_fLog, ScoreMxAllBg);

	float ES = GetExpectedScore(ScoreMxAlnBg, TrainTPLetterFreqs);
	float HAln = GetShannonEntropy(TrainAlnLetterPairFreqMx);
	float Hrel = GetRelativeEntropy(TrainAlnLetterPairFreqMx, ScoreMxAlnBg);
	Log("ES=%.3g\n", ES);
	Log("Shannon entropy %.3g\n", HAln);
	Log("Relative entropy %.3g\n", Hrel);

	vector<bool> EvalTPs;
	vector<string> EvalRows;
	vector<string> EvalLabels;
	vector<uint> EvalRowChainIdxs;
	AppendAlns(EvalTPAlnFN, LabelToChainIdx, true,
	  EvalRows, EvalLabels, EvalRowChainIdxs, EvalTPs);
	AppendAlns(EvalFPAlnFN, LabelToChainIdx, false,
	  EvalRows, EvalLabels, EvalRowChainIdxs, EvalTPs);

	vector<uint> EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec;
	GetGapCountVecs(EvalRows, EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec);

	vector<float> EvalAlnSubstScores;
	GetAlnSubstScores(ChainIntSeqs, EvalRows, EvalRowChainIdxs,
		IgnoreUndef, UINT_MAX, ScoreMxAlnBg, EvalAlnSubstScores);

	Log("\nQuarts subst. scores only\n");
	LogAlnScoreQuarts(EvalAlnSubstScores, EvalTPs);

	float OpenPenalty, ExtPenalty, Bias, Area;
	OptimizeArea(EvalAlnSubstScores, EvalAlnColCountVec, EvalAlnOpenVec,
		EvalAlnExtVec, EvalTPs, OpenPenalty, ExtPenalty, Bias, Area);

	Log("Area=%.3g, open %.3g, ext %.3g, bias %.3g\n",
		Area, OpenPenalty, ExtPenalty, Bias);

	vector<float> EvalAlnSubstScores3SigFig, EvalAlnScores3SigFig;
	Round3SigFig(EvalAlnSubstScores, EvalAlnSubstScores3SigFig);

	OpenPenalty = 3.89f;//@@@@@@@@@@@@@@@@@@@
	ExtPenalty = 0.242f;//@@@@@@@@@@@@@@@@@@@
	vector<float> EvalAlnScores;
	GetAlnScores(EvalAlnSubstScores, EvalAlnColCountVec,
		EvalAlnOpenVec, EvalAlnExtVec, OpenPenalty, ExtPenalty, Bias,
		EvalAlnScores);
	Round3SigFig(EvalAlnScores, EvalAlnScores3SigFig);
	{//@@
	vector<float> StepScores, StepTPfs, StepFPfs;
	GetSteps(EvalAlnScores3SigFig, EvalTPs, StepScores, StepTPfs, StepFPfs);
	FILE *f = CreateStdioFile("steps2.tsv");
	for (uint i = 0; i < SIZE(StepScores); ++i)
		fprintf(f, "%.3g\t%.3g\t%.3g\n",
			StepScores[i], StepTPfs[i], StepFPfs[i]);
	CloseStdioFile(f);
	return;
	}//@@

	float Area_SubstScores = CalcArea(EvalAlnSubstScores, EvalTPs);
	float Area_Gaps = CalcArea(EvalAlnScores, EvalTPs);

	Log("\nQuarts with optimized gaps, area=%.3g:\n", Area_Gaps);
	LogAlnScoreQuarts(EvalAlnScores, EvalTPs);


	{//@@@@@@@@@@@@@
	Log("\nQuarts 3 sig fig\n");
	LogAlnScoreQuarts(EvalAlnSubstScores3SigFig, EvalTPs);
	FILE *f = CreateStdioFile("alnscores2.tsv");
	uint n = SIZE(EvalAlnScores3SigFig);
	asserta(SIZE(EvalRowChainIdxs) == 2*n);
	for (uint i = 0; i < n; ++i)
		{
		uint ChainIdx1 = EvalRowChainIdxs[2*i];
		uint ChainIdx2 = EvalRowChainIdxs[2*i+1];
		fprintf(f, "%s\t%s\t%.6g\t%c\n",
			Chains[ChainIdx1]->m_Label.c_str(),
			Chains[ChainIdx2]->m_Label.c_str(),
			EvalAlnSubstScores3SigFig[i],
			tof(EvalTPs[i]));
		}
	CloseStdioFile(f);
	}//@@@@@@@@@@@@@

	//vector<float> StepScores, StepTPfs, StepFPfs;
	//GetSteps(EvalAlnScores3SigFig, EvalTPs, StepScores, StepTPfs, StepFPfs);
	//{
	//FILE *f = CreateStdioFile("steps.tsv");
	//for (uint i = 0; i < SIZE(StepScores); ++i)
	//	fprintf(f, "%.3g\t%.3g\t%.3g\n",
	//		StepScores[i], StepTPfs[i], StepFPfs[i]);
	//CloseStdioFile(f);
	//}

	vector<float> SubstStepScores, SubstStepTPfs, SubstStepFPfs;
	GetSteps(EvalAlnSubstScores3SigFig, EvalTPs, SubstStepScores,
		SubstStepTPfs, SubstStepFPfs);

	{
	FILE *f = CreateStdioFile("stepsub2.tsv");
	for (uint i = 0; i < SIZE(SubstStepScores); ++i)
		fprintf(f, "%.3g\t%.3g\t%.3g\n",
			SubstStepScores[i], SubstStepTPfs[i], SubstStepFPfs[i]);
	CloseStdioFile(f);
	}
	}
