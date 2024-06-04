#include "myutils.h"
#include "scop40bench.h"
#include "sort.h"
#include <set>

static FILE *g_ftsv;

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

void SCOP40Bench::ScoresToEvalues()
	{
// slope -7.313169328241681 intercept 6.098028716536934
	const float Slope = -7.3f;
	const float Intercept = 6.1f;
	const uint HitCount = GetHitCount();
	uint DBSize = 9354;
	for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
		{
		float Score = m_Scores[HitIdx];
		float logNF = Slope*Score + Intercept;
		float NF = powf(10, logNF);
		float EPQ = NF/DBSize;
		m_Scores[HitIdx] = EPQ;
		}
	}

uint SCOP40Bench::GetNTPAtEPQThreshold(const vector<uint> &NTPs,
  const vector<uint> &NFPs, float EPQThreshold) const
	{
	uint ntp = 0;
	const uint N = SIZE(NTPs);
	for (uint Idx = 0; Idx < N; ++Idx)
		{
		float EPQ = float(NFPs[Idx])/m_ChainCount;
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

void SCOP40Bench::ScoreDist(const string &FileName)
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	uint ChainCount = GetChainCount();
	uint PairCount = ChainCount*ChainCount;
	ProgressLog("%u chains, %u pairs, NT=%u NF=%u NF+NT=%u\n",
	  ChainCount, PairCount, m_NT, m_NF, m_NT + m_NF);

	vector<float> Scores;
	vector<uint> NTPs;
	vector<uint> NFPs;
	GetROCSteps(Scores, NTPs, NFPs);

	vector<float> SmoothTPRs;
	vector<float> SmoothFPRs;
	vector<float> SmoothScores;
	vector<uint> SmoothNTPs;
	vector<uint> SmoothNFPs;
	const uint BINS = 100;
	SmoothROCSteps(Scores, NTPs, NFPs, BINS, 1,
	  SmoothScores, SmoothNTPs, SmoothNFPs, SmoothTPRs, SmoothFPRs);

	for (uint Bin = 0; Bin < BINS; ++Bin)
		{
		float Score = Scores[Bin];
		uint NTP = SmoothNTPs[Bin];
		uint NFP = SmoothNFPs[Bin];
		uint AlnCount = NTP + NFP;
		double Sens = double(NTP)/m_NT;
		double FilterEff = double(m_NF - NFP)/m_NF;
		ProgressLog("Score %.1f NT %u NF %u Sens %.3g FilterEff %.3g\n",
		  Score, NTP, NFP, Sens, FilterEff);
		}

	CloseStdioFile(f);
	}

int SCOP40Bench::IsT(uint DomIdx1, uint DomIdx2) const
	{
	assert(DomIdx1 < SIZE(m_DomIdxToFamIdx));
	assert(DomIdx2 < SIZE(m_DomIdxToFamIdx));
	uint FamIdx1 = m_DomIdxToFamIdx[DomIdx1];
	uint FamIdx2 = m_DomIdxToFamIdx[DomIdx2];

	assert(DomIdx1 < SIZE(m_DomIdxToFoldIdx));
	assert(DomIdx2 < SIZE(m_DomIdxToFoldIdx));
	uint FoldIdx1 = m_DomIdxToFoldIdx[DomIdx1];
	uint FoldIdx2 = m_DomIdxToFoldIdx[DomIdx2];

	if (m_Mode == "family")
		{
		if (FamIdx1 == FamIdx2)
			return 1;
		else
			return 0;
		}
	if (m_Mode == "fold")
		{
		if (FoldIdx1 == FoldIdx2)
			return 1;
		else
			return 0;
		}
	if (m_Mode == "ignore")
		{
		if (FoldIdx1 == FoldIdx2)
			{
			if (FamIdx1 == FamIdx2)
				return 1;
			else
				return -1;
			}
		else
			return 0;
		}
	asserta(false);
	return INT_MAX;
	}

void SCOP40Bench::SetTFs()
	{
	uint HitCount = GetHitCount();
	m_TFs.clear();
	m_TFs.reserve(HitCount);
	SetNXs();
	uint DomCount = SIZE(m_Doms);
	uint FamCount = SIZE(m_Fams);
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

void SCOP40Bench::ROCStepsToTsv(const string &FileName,
  vector<float> &Scores, 
  vector<uint> &NTPs, vector<uint> &NFPs,
  vector<float> &TPRs, vector<float> &FPRs) const
	{
	if (FileName.empty())
		return;
	const uint N = SIZE(Scores);
	asserta(SIZE(TPRs) == N);
	asserta(SIZE(FPRs) == N);
	asserta(SIZE(NTPs) == N);
	asserta(SIZE(NFPs) == N);
	float DBSize = (float) GetChainCount();

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
void SCOP40Bench::EPQToTsvX(FILE *f, uint N)
	{
	if (f == 0)
		return;
	vector<float> EPQs(N+1, -1);
	vector<float> BinScores(N+1, -1);
	vector<float> BinErrs(N+1, 99);
	float SensStep = 1.0f/N;

	vector<float> Scores;
	vector<uint> NTPs;
	vector<uint> NFPs;
	GetROCSteps(Scores, NTPs, NFPs);
	uint DBSize = GetChainCount();
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
	fprintf(f, "Bin\tScore\tSens=TPR\tEPQ\n");
	for (uint Bin = 0; Bin <= N; ++Bin)
		{
		fprintf(f, "%u", Bin);
		fprintf(f, "\t%.3g", BinScores[Bin]);
		fprintf(f, "\t%.3f", Bin*SensStep);
		fprintf(f, "\t%.3g", EPQs[Bin]);
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
	ProgressLog("SmoothSteps\n");
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
  vector<uint> &NTPs, vector<uint> &NFPs)
	{
	SetNXs();
	Scores.clear();
	NTPs.clear();
	NFPs.clear();
	const uint HitCount = GetHitCount();
	if (HitCount == 0)
		return;

	asserta(SIZE(m_TFs) == HitCount);
	ProgressLogPrefix("Sort scores\n");
	SetScoreOrder();
	const vector<uint> &Order = m_ScoreOrder;
	asserta(SIZE(Order) == HitCount);

	float CurrentScore = m_Scores[Order[0]];
	uint NTP = 0;
	uint NFP = 0;
	ProgressLogPrefix("ROC Steps\n");
	for (uint k = 0; k < HitCount; ++k)
		{
		uint i = Order[k];
		if (m_DomIdx1s[i] == m_DomIdx2s[i])
			continue;
		float Score = m_Scores[i];
		if (Score != CurrentScore)
			{
			Scores.push_back(CurrentScore);
			NTPs.push_back(NTP);
			NFPs.push_back(NFP);
			CurrentScore = Score;
			}
		int T = m_TFs[i];
		if (T == 1)
			++NTP;
		else if (T == 0)
			++NFP;
		}
	Scores.push_back(CurrentScore);
	NTPs.push_back(NTP);
	NFPs.push_back(NFP);
	}

uint SCOP40Bench::GetFamSize(uint FamIdx) const
	{
	uint Size = 0;
	uint DomCount = GetDomCount();
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		uint FamIdx2 = m_DomIdxToFamIdx[DomIdx];
		if (FamIdx2 == FamIdx)
			++Size;
		}
	return Size;
	}

void SCOP40Bench::GetFamSizes(vector<uint> &FamSizes) const
	{
	FamSizes.clear();
	uint FamCount = GetFamCount();
	uint DomCount = GetDomCount();
	FamSizes.resize(FamCount, 0);
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		uint FamIdx = m_DomIdxToFamIdx[DomIdx];
		assert(FamIdx < FamCount);
		++FamSizes[FamIdx];
		}
	}

void SCOP40Bench::GetFoldSizes(vector<uint> &FoldSizes) const
	{
	FoldSizes.clear();
	uint FoldCount = GetFoldCount();
	uint DomCount = GetDomCount();
	FoldSizes.resize(FoldCount, 0);
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		uint FoldIdx = m_DomIdxToFoldIdx[DomIdx];
		assert(FoldIdx < FoldCount);
		++FoldSizes[FoldIdx];
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
	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		uint FoldIdx = m_DomIdxToFoldIdx[DomIdx];
		asserta(FoldIdx < FoldCount);
		FoldToDoms[FoldIdx].push_back(DomIdx);
		}

	uint NonSelfPairCount = DomCount*DomCount - DomCount;

	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		uint FamIdx = m_DomIdxToFamIdx[DomIdx];
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
			uint FamIdx2 = m_DomIdxToFamIdx[DomIdx2];
			uint FoldIdx2 = m_DomIdxToFoldIdx[DomIdx2];
			if (m_Mode == "family")
				{
				if (FamIdx2 == FamIdx)
					++m_NT;
				}
			else if (m_Mode == "fold")
				{
				if (FoldIdx2 == FoldIdx)
					++m_NT;
				}
			else if (m_Mode == "ignore")
				{
				if (FamIdx2 == FamIdx)
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
	m_DomIdxToFamIdx.resize(DomCount, UINT_MAX);
	m_DomIdx1s.resize(HitCount, UINT_MAX);
	m_DomIdx2s.resize(HitCount, UINT_MAX);
	m_Scores.resize(HitCount, FLT_MAX);
	ReadStdioFile(f, m_DomIdxToFamIdx.data(), DomCount*sizeof(uint));
	ReadStdioFile(f, m_DomIdx1s.data(), HitCount*sizeof(uint));
	ReadStdioFile(f, m_DomIdx2s.data(), HitCount*sizeof(uint));
	ReadStdioFile(f, m_Scores.data(), HitCount*sizeof(float));
	CloseStdioFile(f);
	}

void SCOP40Bench::LoadHitsFromTsv(const string &FileName)
	{
	asserta(!m_Doms.empty());
	asserta(!m_Fams.empty());
	asserta(!m_DomToIdx.empty());
	asserta(!m_FamToIdx.empty());
	asserta(!m_DomIdxToFamIdx.empty());

	uint ScoreFieldNr = (optset_scorefieldnr ? opt_scorefieldnr-1 : 2);

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
			string Fam1;
			string Fam2;
			GetDomFamFromDomSlashFam(Label1, Dom1, Fam1);
			GetDomFamFromDomSlashFam(Label2, Dom2, Fam2);
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

void SCOP40Bench::LoadLabels(const string &FileName)
	{
	asserta(!FileName.empty());

	m_Doms.clear();
	m_DomToIdx.clear();
	m_Fams.clear();
	m_Folds.clear();
	m_FamToIdx.clear();
	m_FoldToIdx.clear();
	m_DomIdxToFamIdx.clear();
	m_DomIdxToFoldIdx.clear();

	vector<string> Labels;
	ReadLinesFromFile(FileName, Labels);
	m_DomToIdx.clear();
	for (uint Idx = 0; Idx < SIZE(Labels); ++Idx)
		{
		const string &Label = Labels[Idx];

		string Dom;
		string Fam;
		string Fold;
		GetDomFamFoldFromDomSlashFam(Label, Dom, Fam, Fold);

		uint FamIdx = UINT_MAX;
		if (m_FamToIdx.find(Fam) == m_FamToIdx.end())
			{
			FamIdx = SIZE(m_Fams);
			m_Fams.push_back(Fam);
			m_FamToIdx[Fam] = FamIdx;
			}
		else
			FamIdx = m_FamToIdx[Fam];

		uint FoldIdx = UINT_MAX;
		if (m_FoldToIdx.find(Fold) == m_FoldToIdx.end())
			{
			FoldIdx = SIZE(m_Folds);
			m_Folds.push_back(Fold);
			m_FoldToIdx[Fold] = FoldIdx;
			}
		else
			FoldIdx = m_FoldToIdx[Fold];
 
		uint DomIdx = UINT_MAX;
		if (m_DomToIdx.find(Dom) == m_DomToIdx.end())
			{
			DomIdx = SIZE(m_Doms);
			m_Doms.push_back(Dom + "/" + Fam);
			m_DomToIdx[Dom] = DomIdx;
			m_DomIdxToFamIdx.push_back(FamIdx);
			}
		else
			DomIdx = m_DomToIdx[Dom];

		m_DomIdxs.push_back(DomIdx);
		}
	}

void SCOP40Bench::LogFams() const
	{
	const uint FamCount = GetFamCount();
	vector<uint> Sizes;
	GetFamSizes(Sizes);
	asserta(SIZE(Sizes) == FamCount);
	vector<uint> Order(FamCount);
	QuickSortOrderDesc(Sizes.data(), FamCount, Order.data());
	Log("%u fams\n", FamCount);
	for (uint k = 0; k < FamCount; ++k)
		{
		uint FamIdx = Order[k];
		uint Size = Sizes[FamIdx];
		Log("[%5u]  %s\n", Size, m_Fams[FamIdx].c_str());
		}
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
	ROCStepsToTsv(FileName, SmoothScores, NTPs, NFPs,
	  SmoothTPRs, SmoothFPRs);
	}

void SCOP40Bench::EvalEval()
	{
	const uint HitCount = GetHitCount();
// Evalue bins 1e2, 1e1, 1e0, 1e-1, 1e-2 ... 1e-10...1e-M
	const int M = 14;
	const int MM = 1;
	const int MMM = M + MM + 1;
	const float MinE = powf(10.0f, -(M+1));
	const float MaxE = 999;
	vector<uint> NFs(MMM);
	float MinEFound = m_Scores[0];
	float MaxEFound = m_Scores[0];
	for (uint i = 0; i < HitCount; ++i)
		{
		float E = m_Scores[i];
		MinEFound = min(E, MinEFound);
		MaxEFound = max(E, MaxEFound);
		if (E <= MinE || E > MaxE)
			continue;
		float log10E = log10f(E);
		int ilog10E = int(log10E+0.5);
		if (ilog10E < -M || ilog10E > MM)
			continue;
		bool T = m_TFs[i];
		if (!T)
			{
			asserta(ilog10E+M >= 0 && ilog10E+M < M+3);
			NFs[ilog10E+M] += 1;
			}
		}

	uint NF = 0;
	uint DBSize = GetChainCount();
	if (g_ftsv != 0)
		fprintf(g_ftsv, "N\tNF\tEPQ\tlog10(EPQ)\tlog10(E)\n");
	for (int i = 0; i < MMM; ++i)
		{
		int ilog10E = -M + i;
		asserta(ilog10E+M >= 0 && ilog10E+M-1 < M+3);
		uint n = NFs[i];
		NF += n;
		float EPQ = float(NF)/DBSize;
		float log10EPQ = 0;
		if (EPQ > 1e-20)
			log10EPQ = log10f(EPQ);
		ProgressLog("N %10u  NF %10u  EPQ %10.3f  log10EPQ %10.3g  log10E %3d\n",
		  n, NF, EPQ, log10EPQ, ilog10E);
		if (g_ftsv != 0)
			fprintf(g_ftsv, "%u\t%u\t%.3g\t%.3g\t%d\n",
			  n, NF, EPQ, log10EPQ, ilog10E);
		}
	ProgressLog("Evalues %.3g .. %.3g\n", MinEFound, MaxEFound);
	}

void cmd_scop40bit2tsv()
	{
	asserta(optset_output);
	asserta(optset_input);
	asserta(g_ftsv != 0);
	SCOP40Bench SB;
	SB.ReadBit(g_Arg1);
	vector<uint> SavedDomIdxToFamIdx = SB.m_DomIdxToFamIdx;
	SB.m_DomIdxToFamIdx.clear();
	SB.ReadChains(opt_input);
	asserta(SB.m_DomIdxToFamIdx == SavedDomIdxToFamIdx);
	uint Sens = SB.GetSens1stFP();
	const uint HitCount = SB.GetHitCount();
	ProgressLog("%u hits, Sens1FP %u\n", HitCount, Sens);
	for (uint i = 0; i < HitCount; ++i)
		{
		uint DomIdx1 = SB.m_DomIdx1s[i];
		uint DomIdx2 = SB.m_DomIdx2s[i];
		fprintf(g_ftsv, "%s", SB.m_Doms[DomIdx1].c_str());
		fprintf(g_ftsv, "\t%s", SB.m_Doms[DomIdx2].c_str());
		fprintf(g_ftsv, "\t%.6g", SB.m_Scores[i]);
		fprintf(g_ftsv, "\n");
		}
	}

void cmd_scop40bit_roce()
	{
	asserta(optset_input);
	SCOP40Bench SB;
	SB.ReadBit(g_Arg1);
	vector<uint> SavedDomIdxToFamIdx = SB.m_DomIdxToFamIdx;
	SB.m_DomIdxToFamIdx.clear();
	SB.ReadChains(opt_input);
	asserta(SB.m_DomIdxToFamIdx == SavedDomIdxToFamIdx);
	SB.SetTFs();

	SB.ScoresToEvalues();
	SB.m_ScoresAreEvalues = true;
	opt_scores_are_evalues = true;
	optset_scores_are_evalues = true;

	SB.EvalEval();
	SB.WriteBit(opt_output);
	}

void cmd_scop40bit_roc()
	{
	asserta(optset_input);
	SCOP40Bench SB;
	asserta(optset_benchmode);
	SB.m_Mode = string(opt_benchmode);
	SB.ReadBit(g_Arg1);
	vector<uint> SavedDomIdxToFamIdx = SB.m_DomIdxToFamIdx;
	SB.m_DomIdxToFamIdx.clear();
	SB.ReadChains(opt_input);
	SB.BuildDomFamIndexesFromQueryChainLabels();
	SB.SetFamSizes();
	//asserta(SB.m_DomIdxToFamIdx == SavedDomIdxToFamIdx);
	float MaxFPR = 0.01f;
	if (optset_maxfpr)
		MaxFPR = (float) opt_maxfpr;
	SB.SetStats(MaxFPR);
	SB.WriteOutputFiles();
	double FoundFract1 = double(SB.m_DomsWithHomologAndTP1Count)/
	  SB.m_DomsWithHomologCount;

	uint nt_firstfp = SB.m_nt_firstfp;
	uint nt_epq1 = SB.m_nt_epq1;
	float SensFirstFP = float(nt_firstfp)/SB.m_NT;
	float SensEPQ1 = float(nt_epq1)/SB.m_NT;
	ProgressLog("SEPQ1=%.4f", SensEPQ1);
	ProgressLog(" S1FP=%.4f", SensFirstFP);
	ProgressLog(" FF1=%.4f", FoundFract1);
	ProgressLog("\n");
	}

float SCOP40Bench::GetEPQAtEvalueThreshold(const vector<float> &Evalues,
  const vector<uint> &NFPs, float Evalue) const
	{
	float DBSize = (float) GetChainCount();
	const uint N = SIZE(Evalues);
	asserta(SIZE(NFPs) == N);
	for (uint i = 0; i < N; ++i)
		{
		if (Evalues[i] >= Evalue)
			{
			uint NFP = NFPs[i];
			float EPQ = NFP/DBSize;
			return EPQ;
			}
		}
	float EPQ = m_NF/DBSize;
	return EPQ;
	}

float SCOP40Bench::GetEvalueAtEPQThreshold(const vector<float> &Evalues,
  const vector<uint> &NFPs, float EPQ) const
	{
	float DBSize = (float) GetChainCount();
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
	//SB.LoadLabels(opt_labels);
	SB.ReadChains(opt_input);
	SB.LoadHitsFromTsv(g_Arg1);
	SB.WriteBit(opt_output);
	uint HitCount = SB.GetHitCount();
	uint Sens1FP = SB.GetSens1stFP();
	ProgressLog("%u hits, Sens1FP %u\n", HitCount, Sens1FP);
	}

void cmd_scop40bit_thresh()
	{
	asserta(optset_thresh);
	asserta(optset_input);
	SCOP40Bench SB;
	SB.ReadBit(g_Arg1);
	vector<uint> SavedDomIdxToFamIdx = SB.m_DomIdxToFamIdx;
	SB.m_DomIdxToFamIdx.clear();
	SB.ReadChains(opt_input);
	asserta(SB.m_DomIdxToFamIdx == SavedDomIdxToFamIdx);
	SB.SetTFs();
	const uint HitCount = SB.GetHitCount();
	uint NT = 0;
	uint NF = 0;
	for (uint i = 0; i < HitCount; ++i)
		{
		float Score = SB.m_Scores[i];
		if (SB.m_ScoresAreEvalues)
			{
			if (Score > opt_thresh)
				continue;
			}
		else
			{
			if (Score < opt_thresh)
				continue;
			}
		uint DomIdx1 = SB.m_DomIdx1s[i];
		uint DomIdx2 = SB.m_DomIdx2s[i];
		bool T = SB.m_TFs[i];
		if (T)
			++NT;
		else
			++NF;

		if (g_ftsv)
			{
			uint FamIdx1 = SB.m_DomIdxToFamIdx[DomIdx1];
			uint FamIdx2 = SB.m_DomIdxToFamIdx[DomIdx2];
			const string &Fam1 = SB.m_Fams[FamIdx1];
			const string &Fam2 = SB.m_Fams[FamIdx2];
			fprintf(g_ftsv, "%s/%s", SB.m_Doms[DomIdx1].c_str(), Fam1.c_str());
			fprintf(g_ftsv, "\t%s/%s", SB.m_Doms[DomIdx2].c_str(), Fam2.c_str());
			fprintf(g_ftsv, "\t%.6g", Score);
			fprintf(g_ftsv, "\n");
			}
		}
	ProgressLog("thresh %.4g NT %u NF %u\n", opt_thresh, NT, NF);
	}

void cmd_scop40bit_evaleval()
	{
	opt_scores_are_evalues = true;
	optset_scores_are_evalues = true;
	asserta(optset_input);
	SCOP40Bench SB;
	SB.ReadBit(g_Arg1);
	vector<uint> SavedDomIdxToFamIdx = SB.m_DomIdxToFamIdx;
	SB.m_DomIdxToFamIdx.clear();
	SB.ReadChains(opt_input);
	asserta(SB.m_DomIdxToFamIdx == SavedDomIdxToFamIdx);
	SB.SetTFs();
	SB.EvalEval();
	}

void cmd_scop40bit_scoredist()
	{
	asserta(optset_input);
	SCOP40Bench SB;
	SB.ReadBit(g_Arg1);
	vector<uint> SavedDomIdxToFamIdx = SB.m_DomIdxToFamIdx;
	SB.m_DomIdxToFamIdx.clear();
	SB.ReadChains(opt_input);
	asserta(SB.m_DomIdxToFamIdx == SavedDomIdxToFamIdx);
	SB.SetTFs();
	SB.SetScoreOrder();
	SB.ScoreDist(opt_scoredist);
	}
