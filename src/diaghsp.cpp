#include "myutils.h"
#include "diaghsp.h"
#include "diag.h"
#include "alpha.h"

/***
F = highest-scoring suffix x_s .. x_k-1 for some s has score f
B = Best subsequence of x_0 .. x_k-1 has score b

If x_k-1 is included in B:
	if x_k > 0, add x_k to B
	if x_k <= 0 leave B unchanged
else: # X_k-1 is not included in B:
	if f + x_k > 0:
		Add x_k to F
	else:
		F := empty
		f = 0

M = f + x_k
if M > b:
	# add xk to F and replace B by F
	f = M
	b = M
else:
	if M > 0
		f = M
	else:
		# reset F to be empty.
		f = 0

For multiple subsequences:
W. L. Ruzzo and M. Tompa. A linear time algorithm for finding all maximal scoring subsequences.
In Proceedings of the Seventh International Conference on Intelligent Systems for Molecular Biology,
pages 234-241, Heidelberg, Germany, Aug. 1999. AAAI Press.
***/

int DiagHSP::SearchBrute(int d, int &Lo, int &Len) const
	{
	Lo = 0;
	Len = 0;
	diag dg(m_LQ, m_LT);
	int i = dg.getmini(d);
	int j = dg.getminj(d);
	int len = dg.getlen(d);

	int MaxSubseqScore = 0;
	for (int lo = 0; lo < len; ++lo)
		{
		for (int hi = lo; hi < len; ++hi)
			{
			int Score = GetHSPScore(d, lo, hi-lo+1);
			if (Score > MaxSubseqScore)
				{
				MaxSubseqScore = Score;
				Lo = lo;
				Len = hi - lo + 1;
				}
			}
		}
	return MaxSubseqScore;
	}

int DiagHSP::Search(int d, int &Lo, int &Len) const
	{
	diag dg(m_LQ, m_LT);
	int i = dg.getmini(d);
	int j = dg.getminj(d);
	int n = dg.getlen(d);

	int B = 0;
	int F = 0;
	int CurrLen = 0;
	Lo = 0;
	Len = 0;
	int SuffixLo = 0;
	for (int k = 0; k < n; ++k)
		{
		assert(i < m_LQ);
		assert(j < m_LT);
		byte q = m_Q[i++];
		byte t = m_T[j++];
		assert(q < m_AS);
		assert(t < m_AS);
		short Score = m_ScoreMx[q][t];
		F += Score;
		if (F > B)
			{
			B = F;
			Lo = SuffixLo;
			Len = ++CurrLen;
			}
		else if (F > 0)
			++CurrLen;
		else
			{
			F = 0;
			SuffixLo = k+1;
			CurrLen = 0;
			}
		}
	return B;
	}

int DiagHSP::Search_Trace(int d, int &Lo, int &Len) const
	{
	diag dg(m_LQ, m_LT);
	int i = dg.getmini(d);
	int j = dg.getminj(d);
	int n = dg.getlen(d);
	Log("Search_Trace(d=%d) i=%d, j=%d, n=%d\n", d, i, j, n);

	int B = 0;
	int F = 0;
	int CurrLen = 0;
	//int SuffixLo = 0;
	Lo = 0;
	Len = 0;
	int SuffixLo = 0;
	for (int k = 0; k < n; ++k)
		{
		assert(i < m_LQ);
		assert(j < m_LT);
		byte q = m_Q[i++];
		byte t = m_T[j++];
		assert(q < m_AS);
		assert(t < m_AS);
		short Score = m_ScoreMx[q][t];
		F += Score;
		if (F > B)
			{
			B = F;
			Lo = SuffixLo;
			Len = ++CurrLen;
			}
		else if (F > 0)
			++CurrLen;
		else
			{
			F = 0;
			SuffixLo = k+1;
			CurrLen = 0;
			}
		Log(" i=%d j=%d q=%c t=%c Score=%d F=%d B=%d Len=%d",
			i, j, g_LetterToCharMu[q], g_LetterToCharMu[t], Score, F, B, Len);
		Log(" => F=%d B=%d\n", F, B);
		}
	Log(" ====> return B=%d Lo=%d Len=%d\n", B, Lo, Len);
	return B;
	}

int DiagHSP::GetHSPScore_Trace(int d, int lo, int hi) const
	{
	Log("GetHSPScore_Trace(d=%d, lo=%d, hi=%d)\n", d, lo, hi);
	diag dg(m_LQ, m_LT);
	int i = dg.getmini(d) + lo;
	int j = dg.getminj(d) + lo;

	int Total = 0;
	int MaxSuffixScore = 0;
	for (int k = lo; k <= hi; ++k)
		{
		assert(i < m_LQ);
		assert(j < m_LT);
		byte q = m_Q[i];
		byte t = m_T[j];
		assert(q < m_AS);
		assert(t < m_AS);
		Total += m_ScoreMx[q][t];
		for (int k = lo; k <= hi; ++k)
			Log(" [%d]%c [%d]%c += %d\n",
				i, g_LetterToCharMu[q], j, g_LetterToCharMu[t], m_ScoreMx[q][t]);
		++i;
		++j;
		}
	Log("total=%+d\n", Total);
	return Total;
	}

int DiagHSP::GetHSPScore(int d, int lo, int len) const
	{
	diag dg(m_LQ, m_LT);
	int i = dg.getmini(d) + lo;
	int j = dg.getminj(d) + lo;

	int Total = 0;
	for (int k = 0; k < len; ++k)
		{
		assert(i < m_LQ);
		assert(j < m_LT);
		byte q = m_Q[i++];
		byte t = m_T[j++];
		assert(q < m_AS);
		assert(t < m_AS);
		Total += m_ScoreMx[q][t];
		}
	return Total;
	}

static void Test(DiagHSP &DH, const string sQ, const string sT)
	{
	uint LQ = SIZE(sQ);
	uint LT = SIZE(sT);

	byte *Q = myalloc(byte, LQ);
	byte *T = myalloc(byte, LT);

	for (uint i = 0; i < LQ; ++i)
		Q[i] = g_CharToLetterMu[sQ[i]];
	for (uint i = 0; i < LT; ++i)
		T[i] = g_CharToLetterMu[sT[i]];

	diag dg(LQ, LT);
	int mind = dg.getmind();
	int maxd = dg.getmaxd();

	DH.SetQ(Q, LQ);
	DH.SetT(T, LT);

	for (int d = mind; d <= maxd; ++d)
		{
		int BruteLo, BruteLen;
		int BruteScore = DH.SearchBrute(d, BruteLo, BruteLen);
		int CheckScoreBrute = DH.GetHSPScore(d, BruteLo, BruteLen);
		asserta(BruteScore == CheckScoreBrute);

#if 0
		int LoTrace, LenTrace;
		int ScoreTrace = DH.Search_Trace(d, LoTrace, LenTrace);
		int CheckScoreTrace = DH.GetHSPScore(d, LoTrace, LenTrace);
		asserta(ScoreTrace == CheckScoreTrace);
#endif
		int Lo, Len;
		int Score = DH.Search_Trace(d, Lo, Len);
		int CheckScore = DH.GetHSPScore(d, Lo, Len);

#if 0
		ProgressLog("d=%d  score=%d,%d,%d lo=%d,%d  len=%d,%d\n",
					d, BruteScore, Score, CheckScore, BruteLo, Lo, BruteLen, Len);
#endif
		if (Score != CheckScore)
			{
			DH.GetHSPScore_Trace(d, BruteLo, BruteLen);
			int lo_notused, len_notused;
			DH.Search_Trace(d, lo_notused, len_notused);
			Die("Score %d != CheckScore %d", Score, CheckScore);
			}
		}
	}

static void GetRandMuSeq(string &s, uint MinL, uint MaxL)
	{
	s.clear();
	uint L = MinL + randu32()%(MaxL - MinL);
	for (uint i = 0; i < L; ++i)
		{
		uint Letter = randu32()%36;
		s += g_LetterToCharMu[Letter];
		}
	}

void cmd_hsptest()
	{
	string sQ = "SOMEOTHERSTRING";
	string sT = "SEQVENCE";

	DiagHSP DH;

	const uint Iters = 100;
	for (uint Try = 0; Try < Iters; ++Try)
		{
		ProgressStep(Try, Iters, "Trying");
		GetRandMuSeq(sQ, 10, 20);
		GetRandMuSeq(sT, 10, 20);
		Test(DH, sQ, sT);
		}
	Progress("Test passed ok\n");
	}
