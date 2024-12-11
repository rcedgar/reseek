#include "myutils.h"
#include "diag.h"

static void DiagTest(int LQ, int LT)
	{
	diag D(LQ, LT);

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

	asserta(D.getmind() == mind);
	asserta(D.getmaxd() == maxd);

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
		asserta(D.getminj(d) == minj);
		asserta(D.getmaxj(d) == maxj);
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
		asserta(D.getminj(d) == minj);
		asserta(D.getmaxj(d) == maxj);
		}
		}
	ProgressLog("%d, %d ok\n", LQ, LT);
	}

void cmd_diagtest()
	{
	DiagTest(6, 3);
	DiagTest(123, 456);
	DiagTest(1, 1);
	}
