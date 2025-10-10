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

/***

  EPQ
   |       |
   |       |
E1 ++++++++ 
   +++++++/ 
   ++++++/  
E0 +++++    
   |        
    -------------------------- TPR
	    T0 T1

dArea = (E1 + E0)*(T1 + T0)/4
***/

float SCOP40Bench::GetArea(const vector<float> &TPRs,
	const vector<float> &Log10EPQs)
	{
	const uint N = SIZE(TPRs);
	asserta(SIZE(Log10EPQs) == N);
	float Area = 0;
	for (uint i = 1; i < N; ++i)
		{
		float Ti_1 = TPRs[i-1];
		float Ti = TPRs[i];

		float Ei_1 = Log10EPQs[i-1];
		float Ei = Log10EPQs[i];

		float dA = (Ti + Ti_1)*(Ei - Ei_1)/2;
		Area += dA;
		}
	return Area;
	}

void SCOP40Bench::GetSmoothCurve(const vector<float> &TPRs,
						const vector<float> &Es,
						float dE,
						vector<float> &SmoothTPRs,
						vector<float> &SmoothEs) const
	{
	}

void SCOP40Bench::GetCurve(const vector<float> &Scores,
	const vector<uint> &NTPs,
	const vector<uint> &NFPs,
	float MinEPQ, float MaxEPQ,
	vector<float> &CurveScores,
	vector<float> &CurveTPRs,
	vector<float> &CurveEPQs,
	vector<float> &CurveLog10EPQs) const
	{
	CurveScores.clear();
	CurveTPRs.clear();
	CurveEPQs.clear();
	CurveLog10EPQs.clear();

	const uint QueryCount = SIZE(m_Doms);
	const uint N = SIZE(Scores);
	const float StartScore = (m_ScoresAreEvalues ? 0 : FLT_MAX);
	float LastScore = StartScore;
	float LastTPR = 0;
	float LastEPQ = 0;
	float Sum = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint NTP = NTPs[i];
		uint NFP = NFPs[i];

		float Score = Scores[i];
		if (m_ScoresAreEvalues)
			asserta(i == 0 || Score > LastScore);
		else
			asserta(Score < LastScore);
		float TPR = NTP/float(m_NT);
		float EPQ = NFP/float(QueryCount);
		if (TPR == LastTPR || EPQ == LastEPQ || EPQ < MinEPQ)
			{
			LastScore = Score;
			LastTPR = TPR;
			LastEPQ = EPQ;
			continue;
			}
		float Log10EPQ = log10f(EPQ);
		if (EPQ >= MinEPQ && LastEPQ < MinEPQ)
			{
			if (i > 0)
				{
				CurveScores.push_back(LastScore);
				CurveTPRs.push_back(LastTPR);
				CurveEPQs.push_back(LastEPQ);
				if (LastEPQ > 0)
					CurveLog10EPQs.push_back(log10f(LastEPQ));
				else
					CurveLog10EPQs.push_back(0);
				}
			}

		if (EPQ >= MinEPQ && LastEPQ <= MaxEPQ)
			{
			CurveScores.push_back(Score);
			CurveTPRs.push_back(TPR);
			CurveEPQs.push_back(EPQ);
			CurveLog10EPQs.push_back(Log10EPQ);
			if (LastEPQ >= MaxEPQ)
				break;
			}

		LastScore = Score;
		LastTPR = TPR;
		LastEPQ = EPQ;
		}
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
		ProgressStep(i, HitCount, "Set TFs");
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
	if (m_ScoresAreEvalues)
		QuickSortOrder(m_Scores.data(), HitCount, m_ScoreOrder.data());
	else
		QuickSortOrderDesc(m_Scores.data(), HitCount, m_ScoreOrder.data());
	}

void SCOP40Bench::SetTSOrder()
	{
	const uint HitCount = GetHitCount();
	asserta(SIZE(m_TSs) == HitCount);
	m_TSOrder.resize(HitCount);
	QuickSortOrderDesc(m_TSs.data(), HitCount, m_TSOrder.data());
	}

void SCOP40Bench::ROCStepsToTsv(const string &FileName,
  const vector<float> &Scores, 
  const vector<uint> &NTPs, const vector<uint> &NFPs) const
	{
	if (FileName.empty())
		return;
	const uint N = SIZE(Scores);
	asserta(SIZE(NTPs) == N);
	asserta(SIZE(NFPs) == N);
	float DBSize = (float) SIZE(m_Doms);

	FILE *f = CreateStdioFile(FileName);
	fprintf(f, "Score\tNTP\tNFP\n");
	for (uint i = 0; i < N; ++i)
		fprintf(f, "%.8e\t%u\t%u\n", Scores[i], NTPs[i], NFPs[i]);
	CloseStdioFile(f);
	}


void SCOP40Bench::SmoothROCStepsToTsv(const string &FileName,
  const vector<float> &Scores, 
  const vector<uint> &NTPs, const vector<uint> &NFPs,
  const vector<float> &TPRs, const vector<float> &FPRs) const
	{
	if (FileName.empty())
		return;
	const uint N = SIZE(Scores);
	asserta(SIZE(TPRs) == N);
	asserta(SIZE(FPRs) == N);
	asserta(SIZE(NTPs) == N);
	asserta(SIZE(NFPs) == N);
	float DBSize = (float) SIZE(m_Doms);

	FILE *f = CreateStdioFile(FileName);
	fprintf(f, "Score\tNTP\tNFP\tTPR\tFPR\tTPQ\tEPQ\n");
	for (uint i = 0; i < N; ++i)
		{
		float TPQ = float((NTPs[i])/DBSize);
		float EPQ = float(NFPs[i]/DBSize);
		fprintf(f, "%.4g\t%u\t%u\t%.4g\t%.4g\t%.4g\t%.4g\n",
		  Scores[i], NTPs[i], NFPs[i], TPRs[i], FPRs[i], TPQ, EPQ);
		}
	CloseStdioFile(f);
	}

// Project onto common X axis (Sensitivity=TPR) 
//  with N+1 ticks
void SCOP40Bench::WriteCVE(FILE *f, uint N)
	{
	if (f == 0)
		return;
	vector<float> EPQs(N+1, -1);
	vector<float> BinScores(N+1, FLT_MAX);
	vector<float> BinErrs(N+1, 99);
	float SensStep = 1.0f/N;

	vector<float> Scores;
	vector<uint> NTPs;
	vector<uint> NFPs;
	GetROCSteps(Scores, NTPs, NFPs);
	uint DBSize = SIZE(m_Doms);
	const uint NS = SIZE(Scores);
	for (uint i = 0; i < NS; ++i)
		{
		float Score = Scores[i];
		uint ntp = NTPs[i];
		uint nfp = NFPs[i];
		asserta(ntp <= m_NT);
		uint nfn = m_NT - ntp;
		asserta(ntp + nfn == m_NT);
		float Sens = float(ntp)/(ntp + nfn);
		float Sens2 = float(ntp)/m_NT;
		asserta(feq(Sens, Sens2));
	// per query = divided by number of queries = dbsize
		float EPQ = float(nfp)/DBSize;
		if (Sens > 1)
			Die("Sens=%.5g m_NT=%u m_NF=%u ntp=%u nfn=%u Score=%.3g",
			  Sens, m_NT, m_NF, ntp, nfn, Score);
		uint Bin = uint(Sens/SensStep);
		float Err = fabsf(Sens - Bin*SensStep);
		if (Err < BinErrs[Bin])
			{
			asserta(Bin <= N);
			EPQs[Bin] = EPQ;
			BinScores[Bin] = Score;
			}
		}

	float LastEPQ = 0;
	for (uint Bin = 0; Bin < N; ++Bin)
		{
		float EPQ = EPQs[Bin];
		if (EPQ < 0)
			EPQs[Bin] = LastEPQ;
		else
			LastEPQ = EPQ;
		}

	fprintf(f, "=TPR\tEPQ\tScore/E\n");
	for (uint Bin = 0; Bin <= N; ++Bin)
		{
		float TPR = Bin*SensStep;
		float Score = BinScores[Bin];
		if (Score == FLT_MAX)
			break;
		fprintf(f, "%.3f", TPR);
		fprintf(f, "\t%.3g", EPQs[Bin]);
		fprintf(f, "\t%.3g", Score);
		fprintf(f, "\n");
		}
	}

bool SCOP40Bench::SmoothROCSteps(const vector<float> &Scores,
	const vector<uint> &NTPs, const vector<uint> &NFPs,
	uint N, float MaxFPR, vector<float> &SmoothScores,
	vector<uint> &SmoothNTPs, vector<uint> &SmoothNFPs,
	vector<float> &SmoothTPRs, vector<float> &SmoothFPRs) const
	{
	SmoothScores.clear();
	SmoothTPRs.clear();
	SmoothFPRs.clear();
	SmoothNTPs.clear();
	SmoothNFPs.clear();

	const uint NS = SIZE(Scores);
	if (NS < 100)
		{
		Warning("SmoothROCSteps, %u scores", NS);
		return false;
		}
	uint n = NS - 1;
	for (uint i = 0; i < NS; ++i)
		{
		float FPR = NFPs[i]/float(m_NF);
		if (FPR >= MaxFPR)
			{
			n = i;
			break;
			}
		}
	if (n == 0)
		return false;
	asserta(SIZE(NTPs) >= n);
	asserta(SIZE(NFPs) >= n);
	if (n < 2*N)
		return false;
	Progress("SmoothSteps\n");
	const uint HitCount = GetHitCount();
	for (uint Bin = 0; Bin < N; ++Bin)
		{
		uint Idx = UINT_MAX;
		if (Bin == 0)
			Idx = 0;
		else if (Bin + 1 == N)
			Idx = n - 1;
		else
			{
			Idx = (Bin*n)/N;
			asserta(Idx > 0 && Idx < n - 1);
			}

		uint NTP = NTPs[Idx];
		uint NFP = NFPs[Idx];

		SmoothScores.push_back(Scores[Idx]);
		SmoothTPRs.push_back(NTP/float(m_NT));
		SmoothFPRs.push_back(NFP/float(m_NF));
		SmoothNTPs.push_back(NTP);
		SmoothNFPs.push_back(NFP);
		}
	return true;
	}

void SCOP40Bench::GetROCSteps(vector<float> &Scores,
  vector<uint> &NTPs, vector<uint> &NFPs, bool UseTS)
	{
	SetNXs();
	Scores.clear();
	NTPs.clear();
	NFPs.clear();
	const uint HitCount = GetHitCount();
	if (HitCount == 0)
		return;

	asserta(SIZE(m_TFs) == HitCount);
	Progress("Sort scores UseTS=%c m_ScoresAreEvalues=%c\n",
	  tof(UseTS), tof(m_ScoresAreEvalues));
	if (UseTS)
		SetTSOrder();
	else
		SetScoreOrder();
	const vector<uint> &Order = (UseTS ? m_TSOrder : m_ScoreOrder);
	asserta(SIZE(Order) == HitCount);

	float CurrentScore = (UseTS ? m_TSs[Order[0]] : m_Scores[Order[0]]);
	uint NTP = 0;
	uint NFP = 0;
	Progress("ROC Steps\n");
	for (uint k = 0; k < HitCount; ++k)
		{
		uint i = Order[k];
		uint Dom1 = m_DomIdx1s[i];
		uint Dom2 = m_DomIdx2s[i];
		float Score = (UseTS ? m_TSs[i] : m_Scores[i]);
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
			Scores.push_back(CurrentScore);
			NTPs.push_back(NTP);
			NFPs.push_back(NFP);
			CurrentScore = Score;
			}
		if (T == 1)
			++NTP;
		else if (T == 0)
			++NFP;
		}
	Scores.push_back(CurrentScore);
	NTPs.push_back(NTP);
	NFPs.push_back(NFP);
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

void SCOP40Bench::ROCToTsv(const string &FileName, float MaxFPR)
	{
	if (FileName.empty())
		return;
	vector<float> Scores;
	vector<uint> NTPs;
	vector<uint> NFPs;
	GetROCSteps(Scores, NTPs, NFPs);

	vector<float> SmoothTPRs;
	vector<float> SmoothFPRs;
	vector<float> SmoothScores;
	vector<float> TPRThresholds;
	vector<uint> SmoothNTPs;
	vector<uint> SmoothNFPs;
	SmoothROCSteps(Scores, NTPs, NFPs, 100, MaxFPR,
	  SmoothScores, SmoothNTPs, SmoothNFPs, SmoothTPRs, SmoothFPRs);
	SmoothROCStepsToTsv(FileName, SmoothScores, NTPs, NFPs,
	  SmoothTPRs, SmoothFPRs);
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
	uint Sens = SB.GetSens1stFP();
	//SB.LogFirstFewDoms();
	//SB.LogFirstFewHits();
	const uint HitCount = SB.GetHitCount();
	ProgressLog("%u hits, Sens1FP %u\n", HitCount, Sens);
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
	uint Sens1FP = SB.GetSens1stFP();
	ProgressLog("%u hits, Sens1FP %u\n", HitCount, Sens1FP);
	}

void cmd_scop40bench_tsv()
	{
	asserta(optset_lookup);
	SCOP40Bench SB;
	if (opt(scores_are_not_evalues))
		SB.m_ScoresAreEvalues = false;
	else
		SB.m_ScoresAreEvalues = true;
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

void cmd_test()
	{
	vector<float> TPRs;
	vector<float> Log10EPQs;

	uint N = 10;
	//float LoTPR = 0.1f;
	//float HiTPR = 0.4f;
	//float LoE = -2.0f;
	//float HiE = 1.0f;
	float LoTPR = 1;
	float HiTPR = 2;
	float LoE = 1;
	float HiE = 2;
	float CorrectA = (HiTPR + LoTPR)*(HiE - LoE)/2;
	for (uint i = 0; i < N; ++i)
		{
		float TPR = LoTPR + i*(HiTPR - LoTPR)/(N-1);
		float E = LoE + i*(HiE - LoE)/(N-1);
		TPRs.push_back(TPR);
		Log10EPQs.push_back(E);
		}
	float A = SCOP40Bench::GetArea(TPRs, Log10EPQs);
	ProgressLog("A = %.4g, correct = %.4g\n", A, CorrectA);
	Log("");
	}
