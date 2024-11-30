#include "myutils.h"
#include "seqdb.h"
#include "twohitdiag.h"
#include "mudex.h"
#include "mermx.h"
#include "diaghsp.h"
#include "diag.h"
#include "alpha.h"
#include "scop40bench.h"

#define MAP_CHECK	0	// Verify by keeping count map

#if MAP_CHECK
struct DiagInfo
	{
	uint SeqIdxQ = UINT_MAX;
	uint SeqIdxT = UINT_MAX;
	uint16_t Diag = UINT16_MAX;
	bool operator<(const DiagInfo &rhs) const
		{
		if (SeqIdxQ < rhs.SeqIdxQ)
			return true;
		if (SeqIdxQ > rhs.SeqIdxQ)
			return false;
		if (SeqIdxT < rhs.SeqIdxT)
			return true;
		if (SeqIdxT > rhs.SeqIdxT)
			return false;
		return Diag < rhs.Diag;
		}
	};
#endif

void StrToMuLetters(const string &StrSeq, byte *Letters);

void cmd_twohit()
	{
	asserta(optset_hsminscore);
	asserta(optset_mindiagscore);
	asserta(optset_ref);

	const int HSMinScore = opt_hsminscore;
	const int MinDiagScore = opt_mindiagscore;
	const string &RefFN = opt_ref; // "c:/int/reseek_bench/alns/reseek_fast.tsv";

	set<pair<string, string> > RefHitSet;
	uint NTP_ref = 0;
	uint NH_ref = 0;
	{
	FILE *fRef = OpenStdioFile(RefFN);
	string Line;
	vector<string> Fields;
	Progress("Reading ref...");
	while (ReadLineStdioFile(fRef, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		const string &Q = Fields[0];
		const string &T = Fields[1];
		bool IsT = SCOP40Bench::IsTP_SF(Q, T);
		if (IsT)
			++NTP_ref;
		++NH_ref;
		RefHitSet.insert(pair<string, string>(Q, T));
		}
	Progress(" done.\n");
	CloseStdioFile(fRef);
	}

	bool Self = false;
	FILE *fOut = CreateStdioFile(opt_output);

#if MAP_CHECK
	map<DiagInfo, uint> DiagInfoToCount;
#endif

	SeqDB Input;
	Input.FromFasta(g_Arg1);
	Input.SetLabelToIndex();
	const uint SeqCount = Input.GetSeqCount();

	byte **MuSeqs = myalloc(byte *, SeqCount);
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		uint L = Input.GetSeqLength(SeqIdx);
		byte *MuSeq = myalloc(byte, L);
		StrToMuLetters(Input.GetSeq(SeqIdx), MuSeq);
		MuSeqs[SeqIdx] = MuSeq;
		}

	MuDex MD;
	MD.FromSeqDB(Input);
	const uint k = MD.m_k;
	const uint DictSize = MD.m_DictSize;

	MerMx MM;
	extern const short * const *Mu_S_ij_short;
	MM.Init(Mu_S_ij_short, 5, 36, 2);
	asserta(MM.m_AS_pow[5] == DictSize);

	DiagHSP DH;

// Buffer for high-scoring k-mers
	uint *HSKmers = myalloc(uint, DictSize);

	int *SeqIdxTToBestDiagScore = myalloc(int, SeqCount);
	int *SeqIdxTToBestDiag = myalloc(int, SeqCount);
	int *SeqIdxTToBestDiagLo = myalloc(int, SeqCount);
	int *SeqIdxTToBestDiagLen = myalloc(int, SeqCount);
	uint *SeqIdxTs = myalloc(uint, SeqCount);
	uint HitCount = 0;

	zero_array(SeqIdxTToBestDiagScore, SeqCount);
	zero_array(SeqIdxTToBestDiag, SeqCount);
	zero_array(SeqIdxTToBestDiagLo, SeqCount);
	zero_array(SeqIdxTToBestDiagLen, SeqCount);
	zero_array(SeqIdxTs, SeqCount);

	uint NH2 = 0;
	uint NT2 = 0;
	uint NT_both = 0;
	TwoHitDiag TH;
	for (uint SeqIdxQ = 0; SeqIdxQ < SeqCount; ++SeqIdxQ)
		{
		ProgressStep(SeqIdxQ, SeqCount, "Searching");
		const char *SeqQ = Input.GetSeq(SeqIdxQ).c_str();
		const string &LabelQ = Input.GetLabel(SeqIdxQ);
		uint LQ32 = Input.GetSeqLength(SeqIdxQ);
		asserta(LQ32 < UINT16_MAX);
		uint16_t LQ = uint16_t(LQ32);
#if 0
		{
		const char *LabelQ = Input.GetLabel(SeqIdxQ).c_str();
		Log("Q>%s\n", LabelQ);
		}
#endif

	// Initialize k-mer scan of Query
		uint Kmer = 0;
		for (uint SeqPosQ = 0; SeqPosQ < k-1; ++SeqPosQ)
			{
			byte Letter = g_CharToLetterMu[SeqQ[SeqPosQ]];
			assert(Letter < 36);
			Kmer = Kmer*36 + Letter;
			}

	// k-mer scan of Query
		for (uint SeqPosQ = k-1; SeqPosQ < LQ; ++SeqPosQ)
			{
			byte Letter = g_CharToLetterMu[SeqQ[SeqPosQ]];
			assert(Letter < 36);
			Kmer = Kmer*36 + Letter;
			Kmer %= DictSize;
		// Construct high-scoring neighborhood of current k-mer (Kmer)
			const uint HSKmerCount =
				MM.GetHighScoring5mers(Kmer, HSMinScore, HSKmers);
#if 0
			{
			if (HSKmerCount > 0)
				{
				string Tmp;
				Log(" %5u  %s  HSKmerCount=%u\n",
					SeqPosQ-4, MM.KmerToStr(Kmer, k, Tmp), HSKmerCount);
				}
			}
#endif
			for (uint HSKmerIdx = 0; HSKmerIdx < HSKmerCount; ++HSKmerIdx)
				{
				uint HSKmer = HSKmers[HSKmerIdx];
			// Look up HSKmer in MuDex (Mu k-mer index of db)
				uint RowSize = MD.GetRowSize(HSKmer);
				uint DataOffset = MD.GetRowStart(HSKmer);
#if 0
				{
				if (RowSize > 1)
					{
					string Tmp;
					Log("   %s row=%u", MM.KmerToStr(HSKmer, k, Tmp), RowSize);
					}
				}
#endif
				for (uint ColIdx = 0; ColIdx < RowSize; ++ColIdx)
					{
					uint32_t SeqIdxT;
					uint16_t SeqPosT;
					MD.Get(DataOffset++, SeqIdxT, SeqPosT);
					if (SeqIdxT == SeqIdxQ)
						continue;
					if (Self && SeqIdxT < SeqIdxQ)
						continue;
					uint LT32 = Input.GetSeqLength(SeqIdxT);
					asserta(LT32 < UINT16_MAX);
					uint16_t LT = uint16_t(LT32);
					diag dg(LQ, LT);
					uint16_t Diag = dg.getd(SeqPosQ, SeqPosT);
					TH.Add(SeqIdxT, Diag);
#if MAP_CHECK
					{
					DiagInfo DI;
					DI.SeqIdxQ = SeqIdxQ;
					DI.SeqIdxT = SeqIdxT;
					DI.Diag = Diag;
					map<DiagInfo, uint>::iterator iter = DiagInfoToCount.find(DI);
					if (iter == DiagInfoToCount.end())
						DiagInfoToCount[DI] = 1;
					else
						DiagInfoToCount[DI] += 1;
					}
#endif
#if 0
					{
					Log(" %u:%u", SeqIdxT, Diag);
					}
#endif
					}
#if 0
				{
				if (RowSize > 1)
					Log("\n");
				}
#endif
				}
			}
		TH.ClearDupes();
		TH.SetDupes();
#if DEBUG
		TH.Validate(SeqCount, INT16_MAX);
#endif
		vector<pair<uint32_t, uint16_t> > SeqIdxDiagPairs;
		uint DupeCount = TH.m_DupeCount;
		const byte *MuSeqQ = MuSeqs[SeqIdxQ];
		DH.SetQ(MuSeqQ, LQ);
#if 0
		Log("%u dupes\n", DupeCount);
		Log(" Idx   Diag  Score     Lo    Len\n");
		//  12345  12345  12345  12345  12345
#endif
		assert(HitCount == 0);
		for (uint i = 0; i < DupeCount; ++i)
			{
			uint32_t SeqIdxT = TH.m_DupeSeqIdxs[i];
			uint16_t Diag = TH.m_DupeDiags[i];
#if MAP_CHECK
			{
			DiagInfo DI;
			DI.SeqIdxQ = SeqIdxQ;
			DI.SeqIdxT = SeqIdxT;
			DI.Diag = Diag;
			map<DiagInfo, uint>::const_iterator iter = DiagInfoToCount.find(DI);
			asserta(iter != DiagInfoToCount.end());
			asserta(iter->second > 1);
			}
#endif
			const byte *MuSeqT = MuSeqs[SeqIdxT];
			const uint LT = Input.GetSeqLength(SeqIdxT);
			DH.SetT(MuSeqT, LT);
			int Lo, Len;
			int DiagScore = DH.Search(Diag, Lo, Len);
#if 0
			Log("%5u  %5u  %5d  %5d  %5d  >%s (%u)\n",
				SeqIdxT, Diag, DiagScore, Lo, Len,
				Input.GetLabel(SeqIdxT).c_str(), LT);
#endif
			if (DiagScore >= MinDiagScore)
				{
				int BestDiagScoreT = SeqIdxTToBestDiagScore[SeqIdxT];
				if (BestDiagScoreT == 0)
					{
					SeqIdxTs[HitCount++] = SeqIdxT;
					SeqIdxTToBestDiagScore[SeqIdxT] = DiagScore;
					SeqIdxTToBestDiag[SeqIdxT] = Diag;
					SeqIdxTToBestDiagLo[SeqIdxT] = Lo;
					SeqIdxTToBestDiagLen[SeqIdxT] = Len;
					}
				else if (DiagScore > SeqIdxTToBestDiagScore[SeqIdxT])
					{
					SeqIdxTToBestDiagScore[SeqIdxT] = DiagScore;
					SeqIdxTToBestDiag[SeqIdxT] = Diag;
					SeqIdxTToBestDiagLo[SeqIdxT] = Lo;
					SeqIdxTToBestDiagLen[SeqIdxT] = Len;
					}
				}
			}
#if 1
		{
		if (HitCount > 0)
			{
			for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
				{
				uint SeqIdxT = SeqIdxTs[HitIdx];
				const string &LabelT = Input.GetLabel(SeqIdxT);
				bool IsT = SCOP40Bench::IsTP_SF(LabelQ, LabelT);
				++NH2;
				if (IsT)
					{
					++NT2;
					if (RefHitSet.find(pair<string, string>(LabelQ, LabelT)) !=
									   RefHitSet.end())
						++NT_both;
					}
				//int BestDiagScore = SeqIdxTToBestDiagScore[SeqIdxT];
				//int BestDiag = SeqIdxTToBestDiag[SeqIdxT];
				//int BestDiagLo = SeqIdxTToBestDiagLo[SeqIdxT];
				//int BestDiagLen = SeqIdxTToBestDiagLen[SeqIdxT];
				if (fOut != 0)
					fprintf(fOut, "%s\t%s\n", LabelQ.c_str(), LabelT.c_str());
#if LOGDIAGALNS
				Log("%6d  %6d  >%s\n", BestDiagScore, BestDiag, LabelT);
				LogDiagAln(MuSeqQ, LQ, LabelQ,
						   MuSeqT, LT, LabelT,
						   BestDiag, BestDiagLo, BestDiagLen);
#endif

				}
			}
		}
#endif
		for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
			{
			uint SeqIdxT = SeqIdxTs[HitIdx];
			SeqIdxTToBestDiagScore[SeqIdxT] = 0;
			}
#if DEBUG
		{
		for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
			{
			assert(SeqIdxTToBestDiagScore[SeqIdx] == 0);
			}
		}
#endif
		HitCount = 0;
		TH.Reset();
		}
	CloseStdioFile(fOut);
	ProgressLog("NT2=%u", NT2);
	ProgressLog("\tNT_both=%u", NT_both);
	ProgressLog("\tNH=%u", NH2);
	ProgressLog("\tNT_ref=%u", NTP_ref);
	ProgressLog("\ths=%u", opt_hsminscore);
	ProgressLog("\tdg=%u", opt_mindiagscore);
	}
