#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"
#include "scop40bench.h"
#include "sort.h"
#include "quarts.h"
#include <set>

static int MinDiagSpacing = 5;
static int MaxDiagSpacing = 64;
static float MinHSPScore = 0.05f;

uint GetPatternOnes(const string &Str);
void GetKmers(const vector<uint> &Letters,
  const string &PatternStr, uint AlphaSize1, uint AlphaSize,
  vector<uint> &Kmers);
void GetLetters(DSS &D, FEATURE Feature, uint AlphaSize1, 
  vector<uint> &Letters);

static string Symbols =
 "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890@#";

static char LetterToChar(uint Letter)
	{
	if (Letter == UINT_MAX)
		return '~';
	asserta(Letter < SIZE(Symbols));
	return Symbols[Letter];
	}

void GetColPosVecs(const string &Row,
  vector<uint> &ColToPos, vector<uint> &PosToCol)
	{
	ColToPos.clear();
	PosToCol.clear();
	const uint ColCount = SIZE(Row);
	uint Pos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		if (isgap(Row[Col]))
			ColToPos.push_back(UINT_MAX);
		else
			{
			PosToCol.push_back(Col);
			ColToPos.push_back(Pos);
			++Pos;
			}
		}
	}

static void LogPatternAln(
  const vector<uint> &LettersQ, const vector<uint> &LettersR,
  const string &QRow, const string &RRow,
  const vector<uint> &ColToPosQ, vector<uint> &ColToPosR,
  const vector<uint> &KmersQ, const vector<uint> &KmersR)
	{
	const uint ColCount = SIZE(QRow);
	asserta(SIZE(RRow) == ColCount);
	string FRowQ;
	string FRowR;

	uint PosQ = 0;
	uint PosR = 0;
	string AnnotRow(ColCount, ' ');
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char q = QRow[Col];
		char r = RRow[Col];
		if (!isgap(q) && !isgap(r))
			{
			if (PosQ < SIZE(KmersQ) && PosR < SIZE(KmersR))
				{
				if (KmersQ[PosQ] == KmersR[PosR])
					AnnotRow[Col] = '|';
				else
					AnnotRow[Col] = '.';
				}
			}

		if (isgap(q))
			FRowQ.push_back('.');
		else
			{
			asserta(PosQ < SIZE(LettersQ));
			uint Letter = LettersQ[PosQ];
			FRowQ.push_back(LetterToChar(Letter));
			++PosQ;
			}
		if (isgap(r))
			FRowR.push_back('.');
		else
			{
			asserta(PosR < SIZE(LettersR));
			uint Letter = LettersR[PosR];
			FRowR.push_back(LetterToChar(Letter));
			++PosR;
			}
		}

	Log("\n");
	Log("Qaa %s\n", QRow.c_str());
	Log(" QF %s\n", FRowQ.c_str());
	Log("    %s\n", AnnotRow.c_str());
	Log(" RF %s\n", FRowR.c_str());
	Log("Raa %s\n", RRow.c_str());
	}

static void GetAlignedKmers(
  const vector<uint> &ColToPosQ, const vector<uint> &ColToPosR,
  const vector<uint> &KmersQ, const vector<uint> &KmersR,
  vector<uint> &AlignedKmersQ, vector<uint> &AlignedKmersR)
	{
	const uint ColCount = SIZE(ColToPosQ);
	asserta(SIZE(ColToPosR) == ColCount);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint PosQ = ColToPosQ[Col];
		uint PosR = ColToPosR[Col];
		if (PosQ < SIZE(KmersQ) && PosR < SIZE(KmersR))
			{
			AlignedKmersQ.push_back(KmersQ[PosQ]);
			AlignedKmersR.push_back(KmersR[PosR]);
			}
		}
	}

static void MakeKmerToCoords(const vector<uint> &Kmers,
  map<uint, vector<uint> > &KmerToCoords)
	{
	KmerToCoords.clear();
	const uint n = SIZE(Kmers);
	for (uint Pos = 0; Pos < n; ++Pos)
		{
		uint Kmer = Kmers[Pos];
		if (KmerToCoords.find(Kmer) != KmerToCoords.end())
			KmerToCoords[Kmer].push_back(Pos);
		else
			{
			vector<uint> v;
			v.push_back(Pos);
			KmerToCoords[Kmer] = v;
			}
		}
	}

static uint GetSeedCount(DSSAligner &DA,
  const vector<vector<byte> > &ProfileQ,
  const vector<vector<byte> > &ProfileR,
  const vector<uint> &KmersQ,
  const vector<uint> &KmersR,
  const vector<float> &BackgroundKmerFreqs,
  float MaxFreq, uint PatternLength, float MinScore,
  vector<uint> &SeedPosQs, vector<uint> &SeedPosRs)
	{
	uint Count = 0;
	SeedPosQs.clear();
	SeedPosRs.clear();
	map<uint, vector<uint> > KmerToCoordsQ;
	map<uint, vector<uint> > KmerToCoordsR;
	MakeKmerToCoords(KmersQ, KmerToCoordsQ);
	MakeKmerToCoords(KmersR, KmerToCoordsR);
	for (map<uint, vector<uint> >::const_iterator iterq = KmerToCoordsQ.begin();
	  iterq != KmerToCoordsQ.end(); ++iterq)
		{
		uint Kmer = iterq->first;
		if (Kmer == UINT_MAX)
			continue;
		if (BackgroundKmerFreqs[Kmer] > MaxFreq)
			continue;
		map<uint, vector<uint> >::const_iterator iterr = KmerToCoordsR.find(Kmer);
		if (iterr == KmerToCoordsR.end())
			continue;
		asserta(iterr->first == Kmer);
		const vector<uint> &CoordsQ = iterq->second;
		const vector<uint> &CoordsR = iterr->second;
		const uint nq = SIZE(CoordsQ);
		const uint nr = SIZE(CoordsR);
		if (nq*nr > 8)
			continue;
		for (uint iq = 0; iq < nq; ++iq)
			{
			uint PosQ = CoordsQ[iq];
			for (uint ir = 0; ir < nr; ++ir)
				{
				uint PosR = CoordsR[ir];
				float Score = DA.GetScoreSegPair(
				  ProfileQ, ProfileR, PosQ, PosR, PatternLength);
				if (Score >= MinScore)
					{
					SeedPosQs.push_back(PosQ);
					SeedPosRs.push_back(PosR);
					Count += 1;
					}
				}
			}
		}
	return Count;
	}

static bool IsDiagPair(uint PosQi, uint PosQj, uint PosRi, uint PosRj)
	{
	int spacing = abs(int(PosQi) - int(PosQj));
	if (spacing < MinDiagSpacing || spacing > MaxDiagSpacing)
		return false;
	int di = int(PosQi) - int(PosRi);
	int dj = int(PosQj) - int(PosRj);
	if (di == dj)
		return true;
	return false;
	}

static bool GetDiagPairs(
  const vector<uint> &SeedPosQs, const vector<uint> &SeedPosRs,
  vector<uint> &DiagPosQs, vector<uint> &DiagPosRs,
  vector<uint> &Lengths)
	{
	DiagPosQs.clear();
	DiagPosRs.clear();
	Lengths.clear();
	const uint n = SIZE(SeedPosQs);
	asserta(SIZE(SeedPosRs) == n);
	for (uint i = 0; i < n; ++i)
		{
		uint PosQi = SeedPosQs[i];
		uint PosRi = SeedPosRs[i];
		for (uint j = i+1; j < n; ++j)
			{
			uint PosQj = SeedPosQs[j];
			uint PosRj = SeedPosRs[j];
			if (IsDiagPair(PosQi, PosQj, PosRi, PosRj))
				{
				DiagPosQs.push_back(min(PosQi, PosQj));
				DiagPosRs.push_back(min(PosRi, PosRj));
				Lengths.push_back(max(PosQi, PosQj) - min(PosQi, PosQj));
				}
			}
		}
	return SIZE(Lengths) > 0;
	}

void cmd_kmer_eval()
	{
	const float MAXFX = optset_maxfx ? (float) opt_maxfx : 4.0f;
	const float MinSeedScore = optset_minseedscore ? (float) opt_minseedscore : 0.01f;
	const string &PatternStr = optset_pattern ? string(opt_pattern) : "10001";

	vector<FEATURE> ComboFeatures;
	asserta(optset_features);
	vector<string> Fields;
	Split(opt_features, Fields, '_');
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		FEATURE F = StrToFeature(Fields[i].c_str());
		ComboFeatures.push_back(F);
		}
	DSSParams::SetComboFeatures(ComboFeatures);

	uint PatternLength = SIZE(PatternStr);
	uint PatternOnes = GetPatternOnes(PatternStr);

	uint AlphaSize1 = DSS::GetAlphaSize(FEATURE_Combo);
	uint AlphaSize = AlphaSize1;
	for (uint i = 1; i < PatternOnes; ++i)
		AlphaSize *= AlphaSize1;

	DSSParams Params;
	Params.SetDefaults();
	DSS D;
	D.m_Params = &Params;

	DSSAligner DA;
	DA.m_Params = &Params;

	SeqDB Input;
	Input.FromFasta(g_Arg1, true);

	vector<PDBChain *> Chains;
	ReadChains(opt_train_cal, Chains);
	const uint ChainCount = SIZE(Chains);
	map<string, uint> DomToChainIndex;

	vector<uint> BackgroundKmerCounts(AlphaSize);
	uint TotalBackgroundKmerCount = 0;

	vector<vector<uint> > LettersVec(ChainCount);
	vector<vector<uint> > KmersVec(ChainCount);
	vector<vector<vector<byte> > > Profiles(ChainCount);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Set profiles");
		vector<uint> &Letters = LettersVec[ChainIndex];
		vector<uint> &Kmers = KmersVec[ChainIndex];
		vector<vector<byte> > &Profile = Profiles[ChainIndex];

		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Label = Chain.m_Label;
		vector<string> Fields;
		string Dom;
		SCOP40Bench::GetDomFromLabel(Label, Dom);
		DomToChainIndex[Dom] = ChainIndex;

		const uint QL = Chain.GetSeqLength();
		D.Init(Chain);
		D.GetProfile(Profile);
		GetLetters(D, FEATURE_Combo, AlphaSize1, Letters);
		GetKmers(Letters, PatternStr, AlphaSize1, AlphaSize, Kmers);
		for (uint i = 0; i < SIZE(Kmers); ++i)
			{
			uint Kmer = Kmers[i];
			if (Kmer == UINT_MAX)
				continue;
			asserta(Kmer < AlphaSize);
			BackgroundKmerCounts[Kmer] += 1;
			TotalBackgroundKmerCount += 1;
			}
		}
	vector<float> BackgroundKmerFreqs;
	for (uint Kmer = 0; Kmer < AlphaSize; ++Kmer)
		BackgroundKmerFreqs.push_back(BackgroundKmerCounts[Kmer]/
		  float(TotalBackgroundKmerCount));
	uint *SortOrder = myalloc(uint, AlphaSize);
	QuickSortOrderDesc(BackgroundKmerCounts.data(), AlphaSize, SortOrder);
	const float Mean = 1.0f/AlphaSize;
	const float MaxFreq = Mean*MAXFX;
	if (0)
		{
		for (uint i = 0; i < AlphaSize; ++i)
			{
			float Freq = BackgroundKmerFreqs[SortOrder[i]];
			float Rat = Freq/Mean;
			Log("%7u  %8.4g  %7.4f\n", i, Freq, Rat);
			}
		Log("Mean = %8.4g\n", Mean);
		}

	const uint SeqCount = Input.GetSeqCount();
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;

	uint TotalFPSeedCount = 0;
	uint PairWithFPDiagCount = 0;
	uint PairWithFPHSPCount = 0;
	uint FPPairCount = 0;
	uint Skips = 0;
	vector<float> FPDiagScores;
	while (FPPairCount < PairCount)
		{
		uint ChainIndexQ = randu32()%ChainCount;
		uint ChainIndexR = randu32()%ChainCount;

		const string &LabelQ = Chains[ChainIndexQ]->m_Label;
		const string &LabelR = Chains[ChainIndexR]->m_Label;
		string DomQ, ClsQ, FoldQ, SFQ, FmyQ;
		string DomR, ClsR, FoldR, SFR, FmyR;
		SCOP40Bench::ParseScopLabel(LabelQ, DomQ, ClsQ, FoldQ, SFQ, FmyQ);
		SCOP40Bench::ParseScopLabel(LabelR, DomR, ClsR, FoldR, SFR, FmyR);
		if (FmyQ == FmyR)
			{
			++Skips;
			continue;
			}

		const vector<uint> &KmersQ = KmersVec[ChainIndexQ];
		const vector<uint> &KmersR = KmersVec[ChainIndexR];

		const vector<vector<byte> > &ProfileQ = Profiles[ChainIndexQ];
		const vector<vector<byte> > &ProfileR = Profiles[ChainIndexR];

		vector<uint> SeedPosQs;
		vector<uint> SeedPosRs;
		uint SeedCount = GetSeedCount(DA, ProfileQ, ProfileR, KmersQ, KmersR,
		  BackgroundKmerFreqs, MaxFreq, PatternLength, MinSeedScore,
		  SeedPosQs, SeedPosRs);

		vector<uint> DiagPosQs;
		vector<uint> DiagPosRs;
		vector<uint> DiagLengths;
		if (GetDiagPairs(SeedPosQs, SeedPosRs, DiagPosQs, DiagPosRs, DiagLengths))
			++PairWithFPDiagCount;

		bool HasHSP = false;
		const uint DiagCount = SIZE(DiagLengths);
		for (uint i = 0; i < DiagCount; ++i)
			{
			uint DiagPosQ = DiagPosQs[i];
			uint DiagPosR = DiagPosRs[i];
			uint Length = DiagLengths[i];
			float Score = DA.GetScoreSegPair(ProfileQ, ProfileR, DiagPosQ, DiagPosR, Length);
			if (Score >= MinHSPScore)
				HasHSP = true;
			FPDiagScores.push_back(Score);
			}

		if (HasHSP)
			PairWithFPHSPCount += 1;
		TotalFPSeedCount += SeedCount;
		++FPPairCount;
		}

	uint TotalTPSeedCount = 0;
	uint PairWithTPDiagCount = 0;
	uint PairWithTPHSPCount = 0;
	vector<float> TPDiagScores;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Processing");

		const string &LabelQ = Input.GetLabel(2*PairIndex);
		const string &LabelR = Input.GetLabel(2*PairIndex+1);

		const string &RowQ = Input.GetSeq(2*PairIndex);
		const string &RowR = Input.GetSeq(2*PairIndex+1);

		string DomQ;
		string DomR;
		SCOP40Bench::GetDomFromLabel(LabelQ, DomQ);
		SCOP40Bench::GetDomFromLabel(LabelR, DomR);
		asserta(DomQ != DomR);

		uint ChainIndexQ = DomToChainIndex[DomQ];
		uint ChainIndexR = DomToChainIndex[DomR];

		const vector<vector<byte> > &ProfileQ = Profiles[ChainIndexQ];
		const vector<vector<byte> > &ProfileR = Profiles[ChainIndexR];

		const vector<uint> &KmersQ = KmersVec[ChainIndexQ];
		const vector<uint> &KmersR = KmersVec[ChainIndexR];

		vector<uint> SeedPosQs;
		vector<uint> SeedPosRs;
		uint SeedCount = GetSeedCount(DA, ProfileQ, ProfileR, KmersQ, KmersR,
		  BackgroundKmerFreqs, MaxFreq, PatternLength, MinSeedScore,
		  SeedPosQs, SeedPosRs);
		TotalTPSeedCount += SeedCount;

		asserta(SIZE(SeedPosQs) == SeedCount);
		asserta(SIZE(SeedPosRs) == SeedCount);

		vector<uint> DiagPosQs;
		vector<uint> DiagPosRs;
		vector<uint> DiagLengths;
		if (GetDiagPairs(SeedPosQs, SeedPosRs, DiagPosQs, DiagPosRs, DiagLengths))
			++PairWithTPDiagCount;

		const uint DiagCount = SIZE(DiagLengths);
		bool HasHSP = false;
		for (uint i = 0; i < DiagCount; ++i)
			{
			uint DiagPosQ = DiagPosQs[i];
			uint DiagPosR = DiagPosRs[i];
			uint Length = DiagLengths[i];
			float Score = DA.GetScoreSegPair(ProfileQ, ProfileR, DiagPosQ, DiagPosR, Length);
			if (Score >= MinHSPScore)
				HasHSP = true;
			TPDiagScores.push_back(Score);
			}
		if (HasHSP)
			PairWithTPHSPCount += 1;
		}

	double TPPct = GetPct(PairWithTPDiagCount, PairCount);
	double FPPct = GetPct(PairWithFPDiagCount, PairCount);
	double SeedsPerPair = (TotalTPSeedCount + TotalTPSeedCount)/float(2*PairCount);
	double Qual = TPPct - FPPct/2;
	double HSPTPPct = GetPct(PairWithTPHSPCount, PairCount);
	double HSPFPPct = GetPct(PairWithFPHSPCount, PairCount);

	ProgressLog("@");
	ProgressLog("\t%.1f", Qual);
	ProgressLog("\t%s", PatternStr.c_str());
	ProgressLog("\t%s", opt_features);
	ProgressLog("\t%.4f", MinSeedScore);
	ProgressLog("\tfx=%.1f", MAXFX);
	ProgressLog("\tminhsp=%.1f", MinHSPScore);
	ProgressLog("\tdiags=%d,%d", MinDiagSpacing, MaxDiagSpacing);
	ProgressLog("\ttpseeds=%u(%.1f%%)", TotalTPSeedCount, TPPct);
	ProgressLog("\tfpseeds=%u(%.1f%%)", TotalFPSeedCount, FPPct);
	ProgressLog("\tseeds/pair=%.1f", SeedsPerPair);
	ProgressLog("\tsz=%u", AlphaSize);
	ProgressLog("\tHSPsTP=%.1f%%", HSPTPPct);
	ProgressLog("\tHSPsFP=%.1f%%", HSPFPPct);
	ProgressLog("\n");

	QuartsFloat QFTP, QFFP;
	GetQuartsFloat(TPDiagScores, QFTP);
	GetQuartsFloat(FPDiagScores, QFFP);
	Log("TP diag scores ");
	QFTP.LogMe();
	Log("FP diag scores ");
	QFFP.LogMe();
	}
