#include "myutils.h"
#include "scop40bench.h"
#include "pdbchain.h"
#include "sort.h"

void cmd_scop40bit()
	{
	const string &FN = g_Arg1;
	asserta(optset_lookup);
	asserta(optset_output);

	SCOP40Bench SB;
	SB.ReadLookup(opt_lookup);
	SB.ReadHits(FN);
	SB.WriteBit(opt_output);
	}
