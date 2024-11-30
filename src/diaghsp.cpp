#include "myutils.h"
#include "diaghsp.h"
#include "diag.h"
#include "alpha.h"
#include "getticks.h"

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

void DiagHSP::FreeQueryProfile()
	{
	if (m_QueryProfile == 0)
		return;
	myfree(m_QueryProfile);
	m_QueryProfile = 0;
	}

void DiagHSP::SetQueryProfile()
	{
	assert(m_QueryProfile == 0);
	assert(m_LQ != 0);
	uint LP = m_LQ*m_AS;
	m_QueryProfile = myalloc(short, LP);
	uint k = 0;
	for (int PosQ = 0; PosQ < m_LQ; ++PosQ)
		{
		byte LetterQ = m_Q[PosQ];
		for (uint LetterT = 0; LetterT < m_AS; ++LetterT)
			m_QueryProfile[k++] = m_ScoreMx[LetterQ][LetterT];
		}
	assert(k == LP);
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

int DiagHSP::Search_Profile(int d, int &Lo, int &Len) const
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
	uint ProfOffset = i*m_AS;
	for (int k = 0; k < n; ++k)
		{
		assert(j < m_LT);
		byte t = m_T[j++];
		assert(t < m_AS);
		short Score = m_QueryProfile[ProfOffset + t];
		F += Score;
		ProfOffset += m_AS;
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
	DH.SetQueryProfile();
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
		int ProfileScore = DH.Search_Profile(d, Lo, Len);
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
		if (ProfileScore != CheckScore)
			Die("ProfileScore %d != CheckScore %d", Score, CheckScore);
		}
	DH.FreeQueryProfile();
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

static byte *GetRandMuByteSeq(uint MinL, uint MaxL, uint &L)
	{
	L = MinL + randu32()%(MaxL - MinL);
	byte *s = myalloc(byte, L);
	for (uint i = 0; i < L; ++i)
		s[i] = randu32()%36;
	return s;
	}

void cmd_hsptest()
	{
	string sQ = "SOMEOTHERSTRING";
	string sT = "SEQVENCE";

	DiagHSP DH;
	Test(DH, sQ, sT);

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

void cmd_hspspeedtest()
	{
	uint TargetCount = 10000;
	vector<const byte *> Ts;
	vector<uint> LTs;

	uint LQ, LT;
	const byte *Q = GetRandMuByteSeq(100, 200, LQ);
	for (uint i = 0; i < TargetCount; ++i)
		{
		const byte *T = GetRandMuByteSeq(100, 200, LT);
		Ts.push_back(T);
		LTs.push_back(LT);
		}

	DiagHSP DH;
	DH.SetQ(Q, LQ);
	TICKS t1 = GetClockTicks();
	for (uint i = 0; i < TargetCount; ++i)
		{
		const byte *T = Ts[i];
		LT = LTs[i];
		DH.SetT(T, LT);

		int iLQ = int(LQ);
		int iLT = int(LT);
		diag dg(iLQ, iLT);

		int dmin = dg.getmind();
		int dmax = dg.getmaxd();
		int d = dmin + int(randu32()%uint(dmax - dmin + 1));
		asserta(d >= dmin && d <= dmax);

		int Lo, Len;
		DH.Search(d, Lo, Len);
		}
	TICKS t2 = GetClockTicks();
	ProgressLog("%.0f ticks\n", double(t2 - t1));
	}
