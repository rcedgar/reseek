#include "myutils.h"
#include "diaghsp.h"
#include "idiag.h"

static void DiagTest(int LQ, int LT)
	{
	idiag D(LQ, LT);

	int mind = INT_MAX;
	int maxd = INT_MAX;
	for (int i = 0; i < LQ; ++i)
		{
		for (int j = 0; j < LT; ++j)
			{
			int d = D.getd(i, j);
			if (mind == INT_MAX || d < mind) mind = d;
			if (maxd == INT_MAX || d > maxd) maxd = d;
			int i2 = D.geti(d, j);
			int j2 = D.getj(d, i);
			asserta(j2 == j);
			asserta(i2 == i);
			}
		}

	int mind2 = D.getmind();
	int maxd2 = D.getmaxd();
	asserta(mind2 == mind);
	asserta(maxd2 == maxd);

	for (int d = mind; d <= maxd; ++d)
		{
		{
		int mini = D.getmini(d);
		int maxi = D.getmaxi(d);
		int minj = INT_MAX;
		int maxj = INT_MAX;
		for (int i = mini; i <= maxi; ++i)
			{
			int j = D.getj(d, i);
			if (minj == INT_MAX || j < minj) minj = j;
			if (maxj == INT_MAX || j > maxj) maxj = j;
			}
		int minj2 = D.getminj(d);
		int maxj2 = D.getmaxj(d);
		asserta(minj2 == minj);
		asserta(maxj2 == maxj);
		}
		{
		int minj = D.getminj(d);
		int maxj = D.getmaxj(d);
		int mini = INT_MAX;
		int maxi = INT_MAX;
		for (int j = minj; j <= maxj; ++j)
			{
			int i = D.geti(d, j);
			if (mini == INT_MAX || i < mini) mini = i;
			if (maxi == INT_MAX || i > maxi) maxi = i;
			}
		int minj2 = D.getminj(d);
		int maxj2 = D.getmaxj(d);
		asserta(minj2 == minj);
		asserta(maxj2 == maxj);
		}
		}
	ProgressLog("%d, %d ok\n", LQ, LT);
	}

void cmd_idiagtest()
	{
	DiagTest(6, 3);
	DiagTest(123, 456);
	DiagTest(1, 1);

	for (uint Try = 0; Try < 100; ++Try)
		{
		int LQ = 1 + randu32()%25;
		int LT = 1 + randu32()%25;
		DiagTest(LT, LQ);
		}

	for (int i = -10; i < 10; ++i)
		{
		uint16_t u = uint16_t(i);
		int i2 = int16_t(u);
		ProgressLog("i=%3d  u=%04x i2=%3d\n", i, u, i2);
		brk(i2 != i);
		}
	}
