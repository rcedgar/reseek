#include "myutils.h"
#include "seqdb.h"
#include "twohitdiag.h"
#include "duper.h"
#include "mudex.h"
#include "mermx.h"
#include "alpha.h"
#include "diaghsp.h"
#include "diag.h"
#include "sort.h"

extern const short * const *Mu_S_ij_short;
void LogDiagAln(const byte *MuSeqQ, uint LQ, const char *LabelQ,
				const byte *MuSeqT, uint LT, const char *LabelT,
				int Diag, int Lo, int Len);

void StrToMuLetters(const string &StrSeq, byte *Letters);

static void DoDiag(const byte *MuSeqQ, uint LQ, const char *LabelQ,
				const byte *MuSeqT, uint LT, const char *LabelT,
				int Diag, int Lo, int Len)
	{
	LogDiagAln(MuSeqQ, LQ, LabelQ, 
				MuSeqT, LT, LabelT,
				Diag, Lo, Len);

	diag dg(LQ, LT);
	int mini = dg.getmini(Diag);
	int minj = dg.getminj(Diag);
	int diaglen = dg.getlen(Diag);
	vector<int> KmerScores;
	for (int diagpos = 0; diagpos + 4 < diaglen; ++diagpos)
		{
		int KmerScore = 0;
		for (int kmerpos = 0; kmerpos < 5; ++kmerpos)
			{
			int posq = mini+diagpos+kmerpos;
			int post = minj+diagpos+kmerpos;
			asserta(posq < int(LQ));
			asserta(post < int(LT));
			byte iq = MuSeqQ[posq];
			byte it = MuSeqT[post];
			KmerScore += Mu_S_ij_short[iq][it];
			}
		KmerScores.push_back(KmerScore);
		}
	const uint N = SIZE(KmerScores);
	asserta(N >= 2);
	uint *Order = myalloc(uint, N);
	QuickSortOrderDesc<int>(KmerScores.data(), N, Order);
	for (uint i = 0; i < 2; ++i)
		{
		int Pos = (int) Order[i];
		int Score = KmerScores[Pos];
		Log("  Pos=%u Score=%d", Pos, Score);
		if (Pos >= Lo && Pos < Lo + Len)
			Log(" in HSP");
		else
			Log(" NOT in HSP");
		Log("\n");
		}
	}

// Systematic search for diagonals found by 2-kmer seeds
void cmd_twohit2()
	{
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

	DiagHSP DH;
	for (uint SeqIdxQ = 0; SeqIdxQ < SeqCount; ++SeqIdxQ)
		{
		const byte *MuSeqQ = MuSeqs[SeqIdxQ];
		const char *LabelQ = Input.GetLabel(SeqIdxQ).c_str();
		uint LQ = Input.GetSeqLength(SeqIdxQ);
		DH.SetQ(MuSeqQ, LQ);
		for (uint SeqIdxT = SeqIdxQ + 1; SeqIdxT < SeqCount; ++SeqIdxT)
			{
			const byte *MuSeqT = MuSeqs[SeqIdxT];
			const char *LabelT = Input.GetLabel(SeqIdxT).c_str();
			uint LT = Input.GetSeqLength(SeqIdxT);
			DH.SetT(MuSeqT, LT);

			diag dg(LQ, LT);
			int mind = dg.getmind();
			int maxd = dg.getmaxd();
			int BestDiagScore = 0;
			int BestDiag = 0;
			int BestLo = 0;
			int BestLen = 0;
			for (int d = mind; d <= maxd; ++d)
				{
				int Lo, Len;
				int DiagScore = DH.Search(d, Lo, Len);
				int CheckScore = DH.GetHSPScore(d, Lo, Len);
				assert(CheckScore == DiagScore);
				if (DiagScore > BestDiagScore)
					{
					BestDiagScore = DiagScore;
					BestDiag = d;
					BestLo = Lo;
					BestLen = Len;
					}
				}
			DoDiag(MuSeqQ, LQ, LabelQ,
				   MuSeqT, LT, LabelT,
				   BestDiag, BestLo, BestLen);
			}
		}
	}
