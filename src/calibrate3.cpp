#include "myutils.h"
#include "dbsearcher.h"
#include "binner.h"
#include "timing.h"
#include <set>

static const uint SAMPLE = 32;
static const uint NBINS = 32;

//	Write small sample of per-query distributions for plotting
void DBSearcher::WriteCalibSample(FILE *f) const
	{
	if (f == 0)
		return;
	asserta(m_CollectTestStats);

// Mid-point TS values for bins
	vector<float> Mids;
	const vector<float> &TSs_0 = m_TestStatsVec[0];
	Binner<float> B(TSs_0, NBINS, -0.1f, 0.3f);
	vector<uint> Bins;
	B.GetAccumBinsReverse(Bins);
	for (uint Bin = 0; Bin < NBINS; ++Bin)
		Mids.push_back(B.GetBinMid(Bin));
 
	vector<vector<uint> > AccumBinsVec(SAMPLE);
	for (uint i = 0; i < SAMPLE; ++i)
		{
		const vector<float> &TSs_i = m_TestStatsVec[i];
		Binner<float> B(TSs_i, NBINS, -0.1f, 0.3f);
		B.GetAccumBinsReverse(AccumBinsVec[i]);
		}

	for (uint Bin = 0; Bin < NBINS; ++Bin)
		{
		fprintf(f, "%.4g", Mids[Bin]);
		for (uint i = 0; i < SAMPLE; ++i)
			{
			uint n = AccumBinsVec[i][Bin];
			if (n == 0)
				fprintf(f, "\t");
			else
				fprintf(f, "\t%u", n);
			}
		fprintf(f, "\n");
		}
	}

uint DBSearcher::GetQueryChainIdx(uint Idx) const
	{
	asserta(Idx < m_QueryChainCount);
	return Idx;
	}

uint DBSearcher::GetDBChainIdx(uint Idx) const
	{
	asserta(Idx < m_DBChainCount);
	if (m_QuerySelf)
		{
		asserta(Idx < m_QueryChainCount);
		return Idx;
		}
	else
		{
		asserta(m_QueryChainCount + Idx < SIZE(m_Chains));
		return m_QueryChainCount + Idx;
		}
	}

const PDBChain &DBSearcher::GetDBChain(uint Idx) const
	{
	uint ChainIdx = GetDBChainIdx(Idx);
	return GetChain(ChainIdx);
	}

const PDBChain &DBSearcher::GetQueryChain(uint Idx) const
	{
	uint ChainIdx = GetQueryChainIdx(Idx);
	return GetChain(ChainIdx);
	}

const char *DBSearcher::GetQueryLabel(uint Idx) const
	{
	const PDBChain &Chain = GetQueryChain(Idx);
	return Chain.m_Label.c_str();
	}

const char *DBSearcher::GetDBLabel(uint Idx) const
	{
	const PDBChain &Chain = GetDBChain(Idx);
	return Chain.m_Label.c_str();
	}

void DBSearcher::WriteCalibOutput(FILE *f) const
	{
	if (f == 0)
		return;
	asserta(m_CollectTestStats);

// Mid-point TS values for bins
	vector<float> Mids;
	const vector<float> &TSs_0 = m_TestStatsVec[0];
	Binner<float> B(TSs_0, NBINS, -0.1f, 0.3f);
	vector<uint> Bins;
	B.GetAccumBinsReverse(Bins);
	fprintf(f, "TS");
	for (uint Bin = 0; Bin < NBINS; ++Bin)
		fprintf(f, "\t%.3g", B.GetBinMid(Bin));
	fprintf(f, "\n");

	const uint NQ = GetQueryCount();
	for (uint i = 0; i < NQ; ++i)
		{
		const char *Label = GetQueryLabel(i);
		const vector<float> &TSs_i = m_TestStatsVec[i];
		Binner<float> B(TSs_i, NBINS, -0.1f, 0.3f);
		B.GetAccumBinsReverse(Bins);

		fprintf(f, "%s", Label);
		for (uint Bin = 0; Bin < NBINS; ++Bin)
			{
			uint n = Bins[Bin];
			fprintf(f, "\t%u", n);
			}
		fprintf(f, "\n");
		}
	}

/***
* Calibrate E-value by measuring score 
* distribution (unsupervised).
***/
void cmd_calibrate3()
	{
	optset_sensitive = true;
	opt_sensitive = true;
	asserta(!optset_fast);

	if (optset_output)
		Die("Use -calib_output");

	string &QCalFN = g_Arg1;
	if (g_Arg1 == ".")
#ifdef _MSC_VER
		QCalFN = "c:/src/reseek_scop40/reseek_db/scop40_family.cal";
#else
		QCalFN = "/c/src/reseek_scop40/reseek_db/scop40_family.cal";
#endif
	else
		QCalFN = g_Arg1;

	const string &DBFN = opt_db;
	DBSearcher DBS;
	DBS.ReadChains(QCalFN, DBFN);
	const uint NQ = DBS.GetQueryCount();
	ProgressLog("NQ=%u\n", NQ);
	set<string> Labels;
	for (uint i = 0; i < NQ; ++i)
		{
		const PDBChain &Chain = *DBS.m_Chains[i];
		const string &Label = Chain.m_Label;
		if (Labels.find(Label) != Labels.end())
			Die("Dupe label >%s", Label.c_str());
		Labels.insert(Label);
		}

	DSSParams Params;
	Params.SetFromCmdLine(DBS.GetDBSize());

	DBS.m_CollectTestStats = true;
	DBS.Setup(Params);
	DBS.Run();
	DSSAligner::Stats();
	DBS.WriteCalibSample(g_fLog);

	if (optset_calib_output)
		{
		FILE *fOut = CreateStdioFile(opt_calib_output);
		DBS.WriteCalibOutput(fOut);
		CloseStdioFile(fOut);
		}
	}
