#include "myutils.h"
#include "scop40bench.h"
#include "sort.h"
#include <set>

static const uint s_LogTopHits = 0;

/***
https://en.wikipedia.org/wiki/Receiver_operating_characteristic

TPR
	= sensitivity
	= recall
	= P(predicted_true | true)
	= TP / (TP + FN)
	= TP / NT

FPR
	= probability of false alarm
	= 1 - specificity
	= P(predicted_false | false)
	= FP / (FP + TN)
	= FP / NF
***/

uint SCOP40Bench::GetNTPAtEPQThreshold(const vector<uint> &NTPs,
  const vector<uint> &NFPs, float EPQThreshold) const
	{
	uint ntp = 0;
	uint QueryCount = SIZE(m_Doms);
	const uint N = SIZE(NTPs);
	for (uint Idx = 0; Idx < N; ++Idx)
		{
		float EPQ = float(NFPs[Idx])/QueryCount;
		if (Idx > 0)
			ntp = NTPs[Idx];
		if (EPQ >= EPQThreshold)
			break;
		}
	return ntp;
	}

float SCOP40Bench::GetTPRAtEPQThreshold(const vector<uint> &NTPs,
  const vector<uint> &NFPs, float EPQThreshold) const
	{
	uint ntp = GetNTPAtEPQThreshold(NTPs, NFPs, EPQThreshold);
	return float(ntp)/m_NT;
	}

// Return 1=TP, 0=FP, -1=ignore
int SCOP40Bench::IsT(uint DomIdx1, uint DomIdx2) const
	{
	if (DomIdx1 == UINT_MAX && DomIdx2 == UINT_MAX)
		return -1;
	if (DomIdx1 == UINT_MAX && DomIdx2 != UINT_MAX)
		return 0;
	if (DomIdx1 != UINT_MAX && DomIdx2 == UINT_MAX)
		return 0;

	if (DomIdx1 == DomIdx2)
		return -1;

	assert(DomIdx1 < SIZE(m_DomIdxToSFIdx));
	assert(DomIdx2 < SIZE(m_DomIdxToSFIdx));
	uint SFIdx1 = m_DomIdxToSFIdx[DomIdx1];
	uint SFIdx2 = m_DomIdxToSFIdx[DomIdx2];

	assert(DomIdx1 < SIZE(m_DomIdxToFoldIdx));
	assert(DomIdx2 < SIZE(m_DomIdxToFoldIdx));
	uint FoldIdx1 = m_DomIdxToFoldIdx[DomIdx1];
	uint FoldIdx2 = m_DomIdxToFoldIdx[DomIdx2];

	if (m_Level == "sf")
		{
		if (SFIdx1 == SFIdx2)
			return 1;
		else
			return 0;
		}
	if (m_Level == "fold")
		{
		if (FoldIdx1 == FoldIdx2)
			return 1;
		else
			return 0;
		}
	if (m_Level == "ignore")
		{
		if (FoldIdx1 == FoldIdx2)
			{
			if (SFIdx1 == SFIdx2)
				return 1;
			else
				return -1;
			}
		else
			return 0;
		}
	Die("IsT(), SCOP40Bench::m_Level='%s'", m_Level.c_str());
	return INT_MAX;
	}

void SCOP40Bench::SetTFs()
	{
	uint HitCount = GetHitCount();
	m_TFs.clear();
	m_TFs.reserve(HitCount);
	SetNXs();
	uint DomCount = SIZE(m_Doms);
	uint SFCount = SIZE(m_SFs);
	for (uint i = 0; i < HitCount; ++i)
		{
		uint DomIdx1 = m_DomIdx1s[i];
		uint DomIdx2 = m_DomIdx2s[i];
		int T = IsT(DomIdx1, DomIdx2);
		m_TFs.push_back(T);
		}
	}

void SCOP40Bench::SetScoreOrder()
	{
	const uint HitCount = GetHitCount();
	m_ScoreOrder.resize(HitCount);
	if (m_SBS == SBS_Evalue)
		QuickSortOrder(m_Scores.data(), HitCount, m_ScoreOrder.data());
	else
		QuickSortOrderDesc(m_Scores.data(), HitCount, m_ScoreOrder.data());
	}

// Project onto common X axis (Sensitivity=TPR) 
//  with N+1 ticks
void SCOP40Bench::WriteCVE(FILE *f) const
	{
	if (f == 0)
		return;
	const uint N = SIZE(m_CVEEPQVec);
	asserta(SIZE(m_CVESensVec) == N);
	asserta(SIZE(m_CVEEPQVec) == N);
	fprintf(f, "TPR\tEPQ\tScore/E\n");
	for (uint i= 0; i < N; ++i)
		{
		fprintf(f, "%.3f", m_CVESensVec[i]);
		fprintf(f, "\t%.3g", m_CVEEPQVec[i]);
		fprintf(f, "\t%.3g", m_CVEScoreVec[i]);
		fprintf(f, "\n");
		}
	}

static float Interp2(float TargetSens,
	float Sens1, float Sens2,
	float EPQ1, float EPQ2)
	{
	asserta(TargetSens != FLT_MAX);
	asserta(Sens1 != FLT_MAX);
	asserta(Sens2 != FLT_MAX);
	asserta(EPQ1 != FLT_MAX);
	asserta(EPQ2 != FLT_MAX);

	if (EPQ1 == EPQ2)
		return EPQ1;

	float r = (Sens2 - Sens1)/(EPQ2 - EPQ1);
	float EPQ = EPQ1 + (TargetSens - Sens1)*r;
	return EPQ;
	}

static float Interp(float TargetSens,
	float PrevSens, float Sens, float NextSens,
	float PrevEPQ, float EPQ, float NextEPQ)
	{
	if (Sens == TargetSens)
		return EPQ;

	float InterpEPQ = FLT_MAX;
	if (PrevSens != FLT_MAX && PrevEPQ != FLT_MAX)
		InterpEPQ = Interp2(TargetSens, PrevSens, Sens, PrevEPQ, EPQ);
	else if (NextSens != FLT_MAX && NextEPQ != FLT_MAX)
		InterpEPQ = Interp2(TargetSens, Sens, NextSens, EPQ, NextEPQ);
	else
		asserta(false);

	//Log("Interp(TargetSens=%.3g, PrevSens=%.3g, Sens=%.3g, NextSens=%.3g PrevEPQ=%.3g, EPQ=%.3g, NextEPQ=%.3g)=%.3g\n",
	//  TargetSens, PrevSens, Sens, NextSens, PrevEPQ, EPQ, NextEPQ, InterpEPQ);
	return InterpEPQ;
	}

void SCOP40Bench::SetArea()
	{
	m_Area0 = FLT_MAX;
	m_Area3 = FLT_MAX;
	const uint N = SIZE(m_CVESensVec);
	asserta(SIZE(m_CVEEPQVec) == N);
	asserta(SIZE(m_CVEScoreVec) == N);
	vector<float> Log10EPQVec;
	for (uint i = 0; i < N; ++i)
		{
		float Sens = m_CVESensVec[i];
		float EPQ = m_CVEEPQVec[i];
		asserta(i == 0 || Sens > m_CVESensVec[i-1]);
		// may not be strictly increasing
		//asserta(i == 0 || EPQ >= m_CVEEPQVec[i-1]);
		asserta(EPQ >= 0.01f);
		asserta(EPQ <= 10.0f);
		float Log10EPQ = log10f(EPQ);
		Log10EPQVec.push_back(Log10EPQ);
		}
	float SensLo = m_CVESensVec[0];
	m_Area0 = 0;
	for (uint i = 0; i + 1 < N; ++i)
		{
		float x1 = m_CVESensVec[i];
		float x2 = m_CVESensVec[i+1];
		float LogEPQ1 = Log10EPQVec[i];
		float LogEPQ2 = Log10EPQVec[i+1];
		asserta(x2 > x1);
		// may not be strictly increasing
		//asserta(LogEPQ2 >= LogEPQ1);
		m_Area0 += (x2 + x1)*fabs(LogEPQ2 - LogEPQ1)/2;
		}

	float SensEPQ0_1 = float(m_nt_epq0_1)/m_NT;
	float SensEPQ1 = float(m_nt_epq1)/m_NT;
	float SensEPQ10 = float(m_nt_epq10)/m_NT;
	m_Area3 = m_Area0 + (SensEPQ0_1 + SensEPQ1 + SensEPQ10)/3;
	}

void SCOP40Bench::SetCVE()
	{
	const uint N = 100;
	const float MinEPQ = 0.01f;
	const float MaxEPQ = 10.0f;

	m_CVESensVec.clear();
	m_CVEEPQVec.clear();
	m_CVEScoreVec.clear();

// First pass to find range IdxLo .. IdxHi
//   such that MinEpq <= EPQ <= MaxEPQ
	const uint StepCount = SIZE(m_ROCStepScores);
	asserta(SIZE(m_ROCStepSenss) == StepCount);
	asserta(SIZE(m_ROCStepEPQs) == StepCount);
	uint IdxLo = UINT_MAX;
	uint IdxHi = UINT_MAX;
	for (uint Idx = 0; Idx < StepCount; ++Idx)
		{
		float Sens = m_ROCStepSenss[Idx];
		float EPQ = m_ROCStepEPQs[Idx];
		if (IdxLo == UINT_MAX && EPQ >= MinEPQ)
			IdxLo = Idx;

		// Need this in case EPQ does not exceed MaxEPQ
		if (EPQ >= MinEPQ && EPQ < MaxEPQ)
			IdxHi = Idx;

		if (Idx > 0 && IdxHi == UINT_MAX && EPQ >= MaxEPQ)
			{
			if (EPQ > MaxEPQ)
				IdxHi = Idx - 1;
			else
				IdxHi = Idx;
			break;
			}
		}
	asserta(IdxLo != UINT_MAX);
	asserta(IdxHi != UINT_MAX);
	asserta(IdxLo < IdxHi);
	asserta(m_ROCStepEPQs[IdxLo] >= MinEPQ);
	asserta(m_ROCStepEPQs[IdxHi] <= MaxEPQ);
	const float SensLo = m_ROCStepSenss[IdxLo];
	const float SensHi = m_ROCStepSenss[IdxHi];
	asserta(SensHi > SensLo);
	float dSens = (SensHi - SensLo)/(N - 1);
	asserta(dSens > 0);

// Second pass to interpolate
	float CurrSens = SensLo;
	for (uint Idx = IdxLo; Idx <= IdxHi; ++Idx)
		{
		float Sens = m_ROCStepSenss[Idx];
		float EPQ = m_ROCStepEPQs[Idx];
		while (Sens >= CurrSens)
			{
			float PrevSens = (Idx == 0 ? FLT_MAX : m_ROCStepSenss[Idx-1]);
			float PrevEPQ = (Idx == 0 ? FLT_MAX : m_ROCStepEPQs[Idx-1]);

			float NextSens = (Idx +1 >= StepCount ? FLT_MAX : m_ROCStepSenss[Idx+1]);
			float NextEPQ = (Idx +1 >= StepCount ? FLT_MAX : m_ROCStepEPQs[Idx+1]);

			float InterpEPQ = Interp(CurrSens,
				PrevSens, Sens, NextSens,
				PrevEPQ, EPQ, NextEPQ);

			m_CVESensVec.push_back(CurrSens);
			m_CVEEPQVec.push_back(InterpEPQ);
			m_CVEScoreVec.push_back(m_ROCStepScores[Idx]);
			CurrSens += dSens;
			}
		}
	while (SIZE(m_CVESensVec) < N)
		{
		float PrevSens = m_ROCStepSenss[IdxHi-1];
		float PrevEPQ = m_ROCStepEPQs[IdxHi-1];
		float Sens = m_ROCStepSenss[IdxHi];
		float EPQ = m_ROCStepEPQs[IdxHi];
		float InterpEPQ = Interp(CurrSens,
			PrevSens, Sens, FLT_MAX,
			PrevEPQ, EPQ, FLT_MAX);
		m_CVESensVec.push_back(CurrSens);
		m_CVEEPQVec.push_back(InterpEPQ);
		m_CVEScoreVec.push_back(m_ROCStepScores[IdxHi]);
		CurrSens += dSens;
		}
	if (!feq(CurrSens, SensHi))
		Warning("SCOP40Bench::SetCVE() CurrSens=%.3g, SensHi=%.3g",
			CurrSens, SensHi);
	asserta(SIZE(m_CVESensVec) == N);
	asserta(SIZE(m_CVEEPQVec) == N);
	asserta(SIZE(m_CVEScoreVec) == N);
	}

void SCOP40Bench::SetROCSteps()
	{
	m_ROCStepScores.clear();
	m_ROCStepNTPs.clear();
	m_ROCStepNFPs.clear();
	m_ROCStepSenss.clear();
	m_ROCStepEPQs.clear();

	SetNXs();
	const uint HitCount = GetHitCount();
	if (HitCount == 0)
		return;

	asserta(SIZE(m_TFs) == HitCount);
	SetScoreOrder();
	const vector<uint> &Order = m_ScoreOrder;
	asserta(SIZE(Order) == HitCount);

	float CurrentScore = m_Scores[Order[0]];
	uint NTP = 0;
	uint NFP = 0;
	for (uint k = 0; k < HitCount; ++k)
		{
		uint i = Order[k];
		uint Dom1 = m_DomIdx1s[i];
		uint Dom2 = m_DomIdx2s[i];
		float Score = m_Scores[i];
		int T = m_TFs[i];
		if (k < s_LogTopHits)
			{
			Log("k=%u", k);
			Log(" doms %s(%u), %s(%u)", 
			  m_Doms[Dom1].c_str(), Dom1,
			  m_Doms[Dom2].c_str(), Dom2);
			Log(" score %.3g", Score);
			Log(" T=%d\n", T);
			}

		if (Dom1 == Dom2)
			continue;
		if (Score != CurrentScore)
			{
			m_ROCStepScores.push_back(CurrentScore);
			m_ROCStepNTPs.push_back(NTP);
			m_ROCStepNFPs.push_back(NFP);
			CurrentScore = Score;
			}
		if (T == 1)
			++NTP;
		else if (T == 0)
			++NFP;
		}
	m_ROCStepScores.push_back(CurrentScore);
	m_ROCStepNTPs.push_back(NTP);
	m_ROCStepNFPs.push_back(NFP);

	const uint StepCount = SIZE(m_ROCStepNTPs);
	asserta(SIZE(m_ROCStepNFPs) == StepCount);
	asserta(SIZE(m_ROCStepScores) == StepCount);

	const uint DBSize = SIZE(m_Doms);
	for (uint StepIdx = 0; StepIdx < StepCount; ++StepIdx)
		{
		uint ntp = m_ROCStepNTPs[StepIdx];
		uint nfp = m_ROCStepNFPs[StepIdx];

		float Sens = float(ntp)/m_NT;
		float EPQ = float(nfp)/DBSize;

		m_ROCStepSenss.push_back(Sens);
		m_ROCStepEPQs.push_back(EPQ);
		}
	}

void SCOP40Bench::SetNXs()
	{
	m_NT = 0;
	m_NF = 0;
	m_NI = 0;
	const uint DomCount = GetDomCount();
	const uint FoldCount = GetFoldCount();
	vector<vector<uint> > FoldToDoms(FoldCount);
	asserta(SIZE(m_DomIdxToFoldIdx) == DomCount);
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		uint FoldIdx = m_DomIdxToFoldIdx[DomIdx];
		asserta(FoldIdx < FoldCount);
		FoldToDoms[FoldIdx].push_back(DomIdx);
		}

	uint NonSelfPairCount = DomCount*DomCount - DomCount;

	asserta(SIZE(m_DomIdxToSFIdx) == DomCount);
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		uint SFIdx = m_DomIdxToSFIdx[DomIdx];
		uint FoldIdx = m_DomIdxToFoldIdx[DomIdx];
		const vector<uint> &FoldDoms = FoldToDoms[FoldIdx];
		bool Found = false;
		for (uint i = 0; i < SIZE(FoldDoms); ++i)
			{
			uint DomIdx2 = FoldDoms[i];
			if (DomIdx2 == DomIdx)
				{
				Found = true;
				continue;
				}
			uint SFIdx2 = m_DomIdxToSFIdx[DomIdx2];
			uint FoldIdx2 = m_DomIdxToFoldIdx[DomIdx2];
			if (m_Level == "sf")
				{
				if (SFIdx2 == SFIdx)
					++m_NT;
				}
			else if (m_Level == "fold")
				{
				if (FoldIdx2 == FoldIdx)
					++m_NT;
				}
			else if (m_Level == "ignore")
				{
				if (SFIdx2 == SFIdx)
					++m_NT;
				else
					++m_NI;
				}
			else
				asserta(false);
			}
		asserta(Found);
		}
	m_NF = NonSelfPairCount - m_NT - m_NI;
	}

void SCOP40Bench::ReadBit(const string &FileName)
	{
	uint DomCount = UINT_MAX;
	uint HitCount = UINT_MAX;
	FILE *f = OpenStdioFile(FileName);
	ReadStdioFile(f, &DomCount, sizeof(DomCount));
	ReadStdioFile(f, &HitCount, sizeof(HitCount));
	Progress("%s hits, %u doms %s\n",
	  IntToStr(HitCount), DomCount, FileName.c_str());
	uint32 FileSize = GetStdioFileSize32(f);
	m_DomIdx1s.resize(HitCount, UINT_MAX);
	m_DomIdx2s.resize(HitCount, UINT_MAX);
	m_Scores.resize(HitCount, FLT_MAX);
	ReadStdioFile(f, m_DomIdx1s.data(), HitCount*sizeof(uint));
	ReadStdioFile(f, m_DomIdx2s.data(), HitCount*sizeof(uint));
	ReadStdioFile(f, m_Scores.data(), HitCount*sizeof(float));
	CloseStdioFile(f);
	}

void SCOP40Bench::LoadHitsFromTsv(const string &FileName)
	{
	asserta(!m_Doms.empty());
	asserta(!m_SFs.empty());
	asserta(!m_DomToIdx.empty());
	asserta(!m_SFToIdx.empty());
	asserta(!m_DomIdxToSFIdx.empty());

	uint ScoreFieldNr = (optset_scorefieldnr ? opt(scorefieldnr)-1 : 2);

	char SplitChar = '\t';
	if (EndsWith(FileName, "dalialn"))
		SplitChar = ' ';
	if (EndsWith(FileName, "tmaln"))
		SplitChar = ' ';

	asserta(!FileName.empty());
	m_DomIdx1s.clear();
	m_DomIdx2s.clear();
	m_Scores.clear();

	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, SplitChar);
		asserta(SIZE(Fields) >= ScoreFieldNr+1);
		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];
		string Dom1;
		string Dom2;
		float Score = (float) StrToFloat(Fields[ScoreFieldNr]);
		if (Label1.find('/') != string::npos)
			{
			string Cls1;
			string Cls2;
			string Fold1;
			string Fold2;
			string SF1;
			string SF2;
			string Fmy1;
			string Fmy2;
			ParseScopLabel(Label1, Dom1, Cls1, Fold1, SF1, Fmy1);
			ParseScopLabel(Label2, Dom2, Cls2, Fold2, SF2, Fmy2);
			}
		else
			{
			Dom1 = Label1;
			Dom2 = Label2;
			}
		map<string, uint>::const_iterator iter1 =
			m_DomToIdx.find(Dom1);
		map<string, uint>::const_iterator iter2 =
			m_DomToIdx.find(Dom2);
		if (iter1 == m_DomToIdx.end() || iter2 == m_DomToIdx.end())
			continue;
		uint DomIdx1 = iter1->second;
		uint DomIdx2 = iter2->second;
		m_DomIdx1s.push_back(DomIdx1);
		m_DomIdx2s.push_back(DomIdx2);
		m_Scores.push_back(Score);
		}
	CloseStdioFile(f);
	}

void cmd_scop40bit2tsv()
	{
	asserta(optset_output);
	asserta(optset_lookup || optset_input);
	FILE *fOut = CreateStdioFile(opt(output));
	SCOP40Bench SB;
	SB.ReadBit(g_Arg1);
	if (optset_lookup)
		SB.ReadLookup(opt(lookup));
	else
		{
		SB.LoadDB(opt(input));
		SB.BuildDomSFIndexesFromDBChainLabels();
		}
	SB.m_Level = "sf";
	const uint HitCount = SB.GetHitCount();
	ProgressLog("%u hits\n", HitCount);
	for (uint i = 0; i < HitCount; ++i)
		{
		uint DomIdx1 = SB.m_DomIdx1s[i];
		uint DomIdx2 = SB.m_DomIdx2s[i];
		const char *Dom1 = SB.m_Doms[DomIdx1].c_str();
		const char *Dom2 = SB.m_Doms[DomIdx2].c_str();
		float Score = SB.m_Scores[i];
		if (optset_tsv_topn)
			{
			if (i >= opt(tsv_topn))
				break;
			Log("%s(%u)", Dom1, DomIdx1);
			Log("  %s(%u)", Dom2, DomIdx2);
			Log("  %.3g", Score);
			Log("\n");
			}
		fprintf(fOut, "%s", Dom1);
		fprintf(fOut, "\t%s", Dom2);
		fprintf(fOut, "\t%.6g", Score);
		fprintf(fOut, "\n");
		}
	CloseStdioFile(fOut);
	}

float SCOP40Bench::GetEPQAtEvalueThreshold(const vector<float> &Evalues,
  const vector<uint> &NFPs, float Evalue) const
	{
	float QueryCount = (float) SIZE(m_Doms);
	const uint N = SIZE(Evalues);
	asserta(SIZE(NFPs) == N);
	for (uint i = 0; i < N; ++i)
		{
		if (Evalues[i] >= Evalue)
			{
			uint NFP = NFPs[i];
			float EPQ = NFP/QueryCount;
			return EPQ;
			}
		}
	float EPQ = m_NF/QueryCount;
	return EPQ;
	}

float SCOP40Bench::GetEvalueAtEPQThreshold(const vector<float> &Evalues,
  const vector<uint> &NFPs, float EPQ) const
	{
	float DBSize = (float) SIZE(m_Doms);
	const uint N = SIZE(Evalues);
	asserta(SIZE(NFPs) == N);
	for (uint i = 0; i < N; ++i)
		{
		uint NFP = NFPs[i];
		float ThisEPQ = NFP/DBSize;
		if (ThisEPQ >= EPQ)
			return Evalues[i];
		}
	return 9e9f;
	}

void cmd_scop40tsv2bit()
	{
	SCOP40Bench SB;
	asserta(optset_input);
	SB.LoadDB(opt(input));
	SB.LoadHitsFromTsv(g_Arg1);
	SB.WriteBit(opt(output));
	uint HitCount = SB.GetHitCount();
	ProgressLog("%u hits\n", HitCount);
	}

void cmd_scop40bench_tsv()
	{
	asserta(optset_lookup);
	SCOP40Bench SB;
	SB.ReadLookup(opt(lookup));
	SB.ReadHits(g_Arg1);
	SB.WriteOutput();
	}

void cmd_scop40bit_roc()
	{
	asserta(optset_lookup || optset_input);
	SCOP40Bench SB;
	SB.ReadBit(g_Arg1);
	if (optset_lookup)
		SB.ReadLookup(opt(lookup));
	else
		{
		SB.LoadDB(opt(input));
		SB.BuildDomSFIndexesFromDBChainLabels();
		}
	SB.WriteOutput();
	}
