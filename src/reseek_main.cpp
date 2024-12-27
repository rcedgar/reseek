#include "myutils.h"
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
	if (optset_myalloc_summary_secs)
		mymalloc_set_summary_secs(opt_myalloc_summary_secs);
	if (optset_myalloc_dump_secs)
		mymalloc_set_dump_secs(opt_myalloc_dump_secs);
	if (optset_myalloc_trace)
		mymalloc_trace(true);

	uint CmdCount = 0;
#define C(x)	if (optset_##x) ++CmdCount;
#include "cmds.h"
	if (CmdCount == 0)
		Die("No command specified");
	if (CmdCount > 1)
		Die("Two commands specified");

#define C(x)	if (optset_##x) { void cmd_##x(); cmd_##x(); }
#include "cmds.h"

	if (optset_myalloc_exit_state)
		{
		mymalloc_write_state("myalloc_exit.state");
		mymalloc_write_map("myalloc_exit.map");
		mymalloc_print_summary("exit");
		}

	LogTiming();
	LogElapsedTimeAndRAM();
	MyutilsExit();
	return 0;
	}
