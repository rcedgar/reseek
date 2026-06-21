#if 0
#include "myutils.h"
#include "features.h"
#include "dssparams.h"

void cmd_dss_logodds()
	{
	asserta(optset_output);
	const string &feature_name = g_Arg1;
	FEATURE F = StrToFeature(g_Arg1.c_str());
	const float * const *scoremx = DSSParams::GetScoreMx(F);
	uint AS = DSSParams::GetAlphaSize(F);
	FILE *f = CreateStdioFile(opt(output));
	string s;
	GetCmdLine(s);
	fprintf(f, "# %s\n", s.c_str());
	fprintf(f, "# [%s]\n", GIT_HASH);
	fprintf(f, "logodds\t%u\n", AS);
	for (uint i = 0; i < AS; ++i)
		{
		fprintf(f, "%u", i);
		for (uint j = 0; j < AS; ++j)
			fprintf(f, "\t%.4g", scoremx[i][j]);
		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	}
#endif
