#include "myutils.h"
#include "tracebit.h"
#include "xdpmem.h"
#include "pathinfo.h"

#define TRACE	0

void TraceBackBitMem(XDPMem &Mem, uint LA, uint LB, char State, string &Path)
	{
	Path.clear();
	byte **TB = Mem.GetTBBit();

#if	TRACE
	Log("\n");
	Log("TraceBackBit\n");
#endif

	size_t i = LA;
	size_t j = LB;
	for (;;)
		{
#if	TRACE
		Log("i=%3d  j=%3d  state=%c\n", (int) i, (int) j, State);
#endif
		if (i == 0 && j == 0)
			break;

		Path += State;

		byte t;
		switch (State)
			{
		case 'M':
			asserta(i > 0 && j > 0);
			t = TB[i-1][j-1];
			if (t & TRACEBITS_DM)
				State = 'D';
			else if (t & TRACEBITS_IM)
				State = 'I';
			else
				State = 'M';
			--i;
			--j;
			break;
		case 'D':
			asserta(i > 0);
			t = TB[i-1][j];
			if (t & TRACEBITS_MD)
				State = 'M';
			else
				State = 'D';
			--i;
			break;

		case 'I':
			asserta(j > 0);
			t = TB[i][j-1];
			if (t & TRACEBITS_MI)
				State = 'M';
			else
				State = 'I';
			--j;
			break;

		default:
			Die("TraceBackBit, invalid state %c", State);
			}
		}
	std::reverse(Path.begin(), Path.end());
	}
