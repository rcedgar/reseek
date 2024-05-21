#include "myutils.h"
// #include "outputfiles.h"
#include "dss.h"
#include "timing.h"

int g_Frame = 0;
string g_Arg1;

void cmd_test() {}

int main(int argc, char **argv)
	{
	MyCmdLine(argc, argv);
	LogProgramInfoAndCmdLine();
	if (!opt_quiet)
		{
		PrintProgramInfo(stdout);
		PrintCopyright(stdout);
		}
	InitTiming();

	uint n = SIZE(g_Argv);
	asserta(n > 0);
	string ShortCmdLine;
	if (n > 1)
		ShortCmdLine = g_Argv[1];
	if (n > 2)
		{
		g_Arg1 = g_Argv[2];
		ShortCmdLine += " " + g_Argv[2];
		}
	if (n > 1)
		{
		ProgressPrefix(false);
		Progress("[%s]\n", ShortCmdLine.c_str() + 1);
		ProgressPrefix(true);
		}

	uint CmdCount = 0;
#define C(x)	if (optset_##x) ++CmdCount;
#include "cmds.h"
	if (CmdCount == 0)
		Die("No command specified");
	if (CmdCount > 1)
		Die("Two commands specified");

#define C(x)	if (optset_##x) { void cmd_##x(); cmd_##x(); }
#include "cmds.h"

	LogTiming();
	LogElapsedTimeAndRAM();

	return 0;
	}
