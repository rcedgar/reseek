#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"
#include "scop40bench.h"
#include <set>

static string Symbols =
 "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890@#";

static char LetterToChar(uint Letter)
	{
	if (Letter == UINT_MAX)
		return '~';
	asserta(Letter < SIZE(Symbols));
	return Symbols[Letter];
	}

static void GetColPosVecs(const string &Row,
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

void GetLetters(DSS &D, FEATURE Feature, uint AlphaSize1,
  vector<uint> &Letters)
	{
	Letters.clear();
	const uint L = D.GetSeqLength();
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		uint Letter = D.GetFeature(Feature, Pos);
		asserta(Letter < AlphaSize1 || Letter == UINT_MAX);
		Letters.push_back(Letter);
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
			AlignedKmersR.push_back(KmersQ[PosR]);
			}
		}
	}

void GetKmers(const vector<uint> &Letters, const string &PatternStr,
  uint AlphaSize1, uint AlphaSize, vector<uint> &Kmers)
	{
	Kmers.clear();
	const uint PatternLength = SIZE(PatternStr);
	const uint L = SIZE(Letters);
	for (uint Pos = 0; Pos + PatternLength < L; ++Pos)
		{
		uint Kmer = 0;
		for (uint j = 0; j < PatternLength; ++j)
			{
			if (PatternStr[j] == '1')
				{
				asserta(Pos + j < SIZE(Letters));
				uint Letter = Letters[Pos + j];
				if (Letter == UINT_MAX)
					{
					Kmer = UINT_MAX;
					break;
					}
				Kmer = Kmer*AlphaSize1 + Letter;
				}
			}
		asserta(Kmer == UINT_MAX || Kmer < AlphaSize);
		Kmers.push_back(Kmer);
		}
	}

static void GetKmerToCoords(const vector<uint> &Kmers, 
  map<uint, vector<uint> > &KmerToCoords)
	{
	KmerToCoords.clear();
	for (uint Pos = 0; Pos < SIZE(Kmers); ++Pos)
		{
		uint Kmer = Kmers[Pos];
		if (Kmer == UINT_MAX)
			continue;
		if (KmerToCoords.find(Kmer) == KmerToCoords.end())
			{
			vector<uint> v;
			v.push_back(Pos);
			KmerToCoords[Kmer] = v;
			}
		else
			KmerToCoords[Kmer].push_back(Pos);
		}
	}

static bool HasDiagSeedPair(
  const map<uint, vector<uint> > &KmerToCoords1,
  const map<uint, vector<uint> > &KmerToCoords2)
	{
	map<int, uint> DiagToCount;
	bool Has2 = false;
	vector<uint> HitCoords1;
	vector<uint> HitCoords2;
	for (map<uint, vector<uint> >::const_iterator iter1 = KmerToCoords1.begin();
	  iter1 != KmerToCoords1.end(); ++iter1)
		{
		uint Kmer = iter1->first;
		const vector<uint> &Coords1 = iter1->second;
		if (SIZE(Coords1) > 2)
			continue;
		map<uint, vector<uint> >::const_iterator iter2 = KmerToCoords2.find(Kmer);
		if (iter2 == KmerToCoords2.end())
			continue;
		asserta(iter2->first == Kmer);
		const vector<uint> &Coords2 = iter2->second;
		if (SIZE(Coords2) > 2)
			continue;
		for (uint i1 = 0; i1 < SIZE(Coords1); ++i1)
			{
			uint Pos1 = Coords1[i1];
			for (uint i2 = 0; i2 < SIZE(Coords2); ++i2)
				{
				uint Pos2 = Coords2[i2];
				HitCoords1.push_back(Pos1);
				HitCoords2.push_back(Pos2);
				}
			}
		}
	const uint N = SIZE(HitCoords1);
	for (uint i = 0; i < N; ++i)
		{
		uint Pos1i = HitCoords1[i];
		uint Pos2i = HitCoords2[i];
		int Diagi = int(Pos1i) - int(Pos2i);
		for (uint j = i+1; j < N; ++j)
			{
			uint Pos1j = HitCoords1[j];
			uint Pos2j = HitCoords2[j];
			int Diagj = int(Pos1j) - int(Pos2j);
			if (Diagi == Diagj)
				{
				int d = abs(int(Pos1i) - int(Pos2i));
				d /= 2;
				if (d >= 3)
					return true;
				}
			}
		}
	return false;
	}

uint GetPatternOnes(const string &Str);

void cmd_kmer_eval2()
	{
	SCOP40Bench SB;
	SB.LoadDB(g_Arg1);

	const string &PatternStr = "11";

	uint PatternLength = SIZE(PatternStr);
	uint PatternOnes = GetPatternOnes(PatternStr);

	uint AlphaSize1 = DSS::GetAlphaSize(FEATURE_Mu);
	uint AlphaSize = AlphaSize1;
	for (uint i = 1; i < PatternOnes; ++i)
		AlphaSize *= AlphaSize1;
	ProgressLog("Pattern %s len %u ones %u alpha size %u\n",
	  PatternStr.c_str(), PatternLength, PatternOnes, AlphaSize);

	const uint ChainCount = SB.GetDBChainCount();
	map<string, uint> DomToChainIndex;

	Warning("non-standard combo features");
	vector<FEATURE> MuFeatures;
	MuFeatures.push_back(FEATURE_NENSS3);
	MuFeatures.push_back(FEATURE_NENDist);
	DSS D;
	DSSParams::SetMuFeatures(MuFeatures);
	
	vector<vector<uint> > LettersVec(ChainCount);
	vector<vector<uint> > KmersVec(ChainCount);
	vector<map<uint, vector<uint> > > KmerToCoordsVec(ChainCount);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *SB.m_DBChains[ChainIndex];

		const uint QL = Chain.GetSeqLength();
		D.Init(Chain);
		vector<uint> &Letters = LettersVec[ChainIndex];
		vector<uint> &Kmers = KmersVec[ChainIndex];
		GetLetters(D, FEATURE_Mu, AlphaSize1, Letters);
		GetKmers(Letters, PatternStr, AlphaSize1, AlphaSize, Kmers);
		GetKmerToCoords(Kmers, KmerToCoordsVec[ChainIndex]);
		}

	const uint PairCount = (ChainCount*(ChainCount-1))/2;
	uint PairIndex = 0;
	uint NT = 0;
	uint NF = 0;
	uint NTSP = 0;
	uint NFSP = 0;
	for (uint ChainIndex1 = 0; ChainIndex1 < ChainCount; ++ChainIndex1)
		{
		const map<uint, vector<uint> > &KmerToCoords1 =
		  KmerToCoordsVec[ChainIndex1];
		for (uint ChainIndex2 = ChainIndex1+1; ChainIndex2 < ChainCount;
		  ++ChainIndex2)
			{
			ProgressStep(PairIndex, PairCount, "Pairs %u %u %u %u",
			  NT, NF, NTSP, NFSP);
			++PairIndex;
			const map<uint, vector<uint> > &KmerToCoords2 =
			  KmerToCoordsVec[ChainIndex2];
			bool SP = HasDiagSeedPair(KmerToCoords1, KmerToCoords2);
			int T = SB.IsT(ChainIndex1, ChainIndex2);
			if (T == 1)
				++NT;
			else if (T == 0)
				++NF;
			if (SP)
				{
				if (T)
					NTSP += 1;
				else
					NFSP += 1;
				}
			}
		}
	ProgressLog("NT %u, NF %u, NTSP %u, NFSP %u\n", NT, NF, NTSP, NFSP);
	}
