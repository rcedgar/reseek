#include "kappa_hsp.h"
#include "flat_params.h"

void kappa_get_hsp_limits(int LQ, int LT, int Diag,
	int &mini, int &minj, int &n)
	{
	asserta(Diag >= 0);
	const int d = Diag;
	mini = LQ - d - 1;
	if (mini < 0)
		mini = 0;
	minj = d + 1 - LQ;
	if (minj < 0)
		minj = 0;
	int maxi = LQ + LT - d - 2;
	if (maxi >= LQ)
		maxi = LQ - 1;
	n = maxi - mini + 1;
	asserta(n > 0);
	}

int kappa_find_hsp(const byte *QSeq, const byte *TSeq,
	int LQ, int LT, int Diag)
	{
	int mini, minj, n;
	kappa_get_hsp_limits(LQ, LT, Diag, mini, minj, n);

	const byte *q = QSeq + mini;
	const byte *t = TSeq + minj;
	const int16_t *mx = kappa32_flat_logodds;
	int B = 0;
	int F = 0;

	int k = 0;
	for (; k + 4 <= n; k += 4)
		{
		for (int u = 0; u < 4; ++u)
			{
			const unsigned bq = q[u];
			const unsigned bt = t[u];
#if !defined(NDEBUG)
			assert(bq < KAPPA_AS);
			assert(bt < KAPPA_AS);
#endif
			const int Score = int(mx[bq * 32u + bt]);
			F += Score;
			if (F > B)
				B = F;
			else if (F < 0)
				F = 0;
			}
		q += 4;
		t += 4;
		}
	for (; k < n; ++k)
		{
		const unsigned bq = *q++;
		const unsigned bt = *t++;
#if !defined(NDEBUG)
		assert(bq < KAPPA_AS);
		assert(bt < KAPPA_AS);
#endif
		const int Score = int(mx[bq * 32u + bt]);
		F += Score;
		if (F > B)
			B = F;
		else if (F < 0)
			F = 0;
		}
	return B;
	}

int kappa_find_hsp2(const byte *QSeq, const byte *TSeq,
	int LQ, int LT, int Diag, int &Lo, int &Len)
	{
	int mini, minj, n;
	kappa_get_hsp_limits(LQ, LT, Diag, mini, minj, n);

	const byte *q = QSeq + mini;
	const byte *t = TSeq + minj;
	int B = 0;
	int F = 0;
	int CurrLen = 0;
	Lo = 0;
	Len = 0;
	int SuffixLo = 0;
	for (int k = 0; k < n; ++k)
		{
		const unsigned bq = *q++;
		const unsigned bt = *t++;
#if !defined(NDEBUG)
		assert(bq < KAPPA_AS);
		assert(bt < KAPPA_AS);
#endif
		const int Score = int(kappa32_flat_logodds[bq * 32u + bt]);
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
			SuffixLo = k + 1;
			CurrLen = 0;
			}
		}
	return B;
	}

int kappa_max_pos_logodds()
	{
	int mx = 0;
	for (uint i = 0; i < 32*32; ++i)
		{
		const int s = int(kappa32_flat_logodds[i]);
		if (s > mx)
			mx = s;
		}
	return mx;
	}
