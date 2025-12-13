#include "myutils.h"
#include "ssslib.h"


void cmd_cluster_sss()
	{
	const string &ChainsFN = g_Arg1;
	asserta(optset_alpha_size);
	asserta(optset_fragl);

	uint AlphaSize = opt(alpha_size);
	uint K = AlphaSize - 1;
	uint FragL = opt(fragl);
	uint BandWidth = 2;
	uint DistN = 20;
	uint FragStep = 8;
	uint MinSize = 100;

	if (optset_bandwidth)
		BandWidth = opt(bandwidth);
	if (optset_fragstep)
		FragStep = opt(fragstep);
	if (optset_distn)
		DistN = opt(distn);
	if (optset_minsize)
		DistN = opt(minsize);

	SSSLib Lib;
	Lib.SetParams(K, FragL, BandWidth, DistN, FragStep);
	Lib.LoadChains(ChainsFN);
	Lib.SetTrainingFrags();
	Lib.AssertSames();
	Lib.Train();
	uint MinClusterSize = Lib.GetMinClusterSize();
	if (MinClusterSize < MinSize)
		Die("MinSize %u", MinClusterSize);
	Lib.AssignReps();
	Lib.LogClusters();
	ProgressLog("MinSize %u (%.2f%%)\n",
		MinClusterSize, GetPct(MinClusterSize, Lib.GetFragCount()));
	Lib.ToSpec(opt(output));
	Lib.AssignLetters();
	Lib.ToFasta(opt(fasta));
	}
