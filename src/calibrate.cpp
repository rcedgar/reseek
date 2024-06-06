#include "myutils.h"
#include "calibratesearcher.h"
#include "mx.h"
#include "pdbchain.h"
#include "xdpmem.h"
#include "cigar.h"
#include "timing.h"
#include <thread>

void CalibrateSearcher::OnSetup()
	{
	m_TestStatsVec.clear();
	m_TestStatsVec.resize(m_ChainCount);
	}

void CalibrateSearcher::OnAln(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA)
	{
	asserta(ChainIndex1 < SIZE(m_TestStatsVec));
	vector<float> &v = m_TestStatsVec[ChainIndex1];
	v.push_back(DA.m_TestStatisticAB);
	}

void cmd_calibrate()
	{
	const string &QCalFN = g_Arg1;
	const string &DBFN = opt_db;
	CalibrateSearcher DBS;
	DSSParams Params;
	Params.SetFromCmdLine();
	DBS.ReadChains(QCalFN, DBFN);
	Params.m_DBSize = (float) DBS.GetDBSize();
	if (optset_dbsize)
		Params.m_DBSize = (float) opt_dbsize;

	DBS.Setup(Params);
	DBS.Run();
	DSSAligner::Stats();
	FILE *f = CreateStdioFile(opt_output);
	const uint ChainCount = SIZE(DBS.m_QueryChains);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		vector<float> &v = DBS.m_TestStatsVec[ChainIndex];
		std::sort(v.begin(), v.end());
		const PDBChain &Chain = *DBS.m_Chains[ChainIndex];
		const char *Label = Chain.m_Label.c_str();
		const uint n = SIZE(v);
		fprintf(f, "%s\t%u", Label, n);
		for (uint i = 0; i < n; ++i)
			fprintf(f, "\t%.3g", v[i]);
		fprintf(f, "\n");
		}
	}
