#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"
#include "alpha.h"

char GetSSConsSymbol(double SSCons);

static const char *GetColor(uint Bin)
	{
	switch (Bin)
		{
	//case 0:	return "splitpea";
	//case 1: return "tv_green";
	//case 2: return "tv_yellow";
	//case 3: return "tv_red";
	//case 4: return "magenta";
	case 0:	return "tv_blue";
	case 1: return "lightblue";
	case 2: return "lightorange";
	case 3: return "salmon";
	case 4: return "tv_red";
		}
	return "white";
	}

static uint GetBin(double SSCons, const vector<double> &Thresholds)
	{
	const uint N = SIZE(Thresholds);
	for (uint i = 0; i < N; ++i)
		{
		if (SSCons <= Thresholds[i])
			return i;
		}
	return N;
	}

double DALIScorer::GetSSConsCol(uint Col, uint w) const
	{
	const uint SeqCount = GetSeqCount();
	vector<vector<vector<double> > > MxVec(SeqCount);
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		uint Pos = m_ColToPosVec[SeqIdx][Col];
		if (Pos == UINT_MAX)
			continue;
		GetDistMxWindow(SeqIdx, Pos, w, MxVec[SeqIdx]);
		}

	uint nr_seqs = GetSeqCount();
	uint nr_cols = GetColCount();
	double sum = 0;
	uint n = 0;
	for (uint seq_idxi = 0; seq_idxi < nr_seqs; ++seq_idxi)
		{
		const vector<vector<double> > &Mxi = MxVec[seq_idxi];
		if (Mxi.empty())
			continue;
		for (uint seq_idxj = seq_idxi + 1; seq_idxj < nr_seqs; ++seq_idxj)
			{
			const vector<vector<double> > &Mxj = MxVec[seq_idxj];
			if (Mxj.empty())
				continue;
			double score = GetLDDTScoreWindow(Mxi, Mxj, w);
			if (score == DBL_MAX)
				continue;
			sum += score;
			++n;
			}
		}
	uint PairCount = (SeqCount*(SeqCount-1))/2;
	double avg = sum/PairCount;
	return avg;
	}

// -msta_sscons MSA -input STRUCTS
void cmd_msta_sscons()
	{
	asserta(optset_input);

	uint w = 2;
	if (optset_window)
		w = opt_window;

	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);
	MSA.SetLabelToIndex();
	const uint ColCount = MSA.GetColCount();

	string Name;
	GetStemName(g_Arg1, Name);

	const bool DoCore = false;
	const bool MissingTestSeqOk = true;

	DALIScorer DS;
	DS.LoadChains(opt_input);
	bool Ok = DS.SetMSA(Name, MSA, DoCore, MissingTestSeqOk);
	if (!Ok)
		Die("SetMSA failed");

	vector<double> SSConsVec;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		double SSCons = DS.GetSSConsCol(Col, w);
		Log("%u", Col);
		if (SSCons == DBL_MAX)
			Log("  .");
		else
			Log("  %.4f", SSCons);
		Log("\n");
		SSConsVec.push_back(SSCons);
		}

	Log("SSCons");
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		double SSCons = SSConsVec[Col];
		Log(", %.3f", SSCons);
		}
	Log("\n");

	vector<double> Thresholds;
	Thresholds.push_back(0.2);
	Thresholds.push_back(0.4);
	Thresholds.push_back(0.6);
	Thresholds.push_back(0.8);

	vector<uint> ColBins;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		double SSCons = SSConsVec[Col];
		uint Bin = GetBin(SSCons, Thresholds);
		ColBins.push_back(Bin);
		}

	if (optset_label)
		{
		const string &Label = string(opt_label);
		const uint QuerySeqIdx = MSA.GetSeqIndex(Label);
		const string &QueryRow = MSA.GetSeq(QuerySeqIdx);
		asserta(SIZE(QueryRow) == ColCount);
		string QuerySeq;
		vector<uint> Bins;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char c = QueryRow[Col];
			if (!isgap(c))
				Bins.push_back(ColBins[Col]);
			}

		const uint LQ = SIZE(Bins);
		uint Start = 0;
		uint CurrentBin = Bins[0];
		Log("select tmp, all\n");
		Log("color %s, tmp\n", GetColor(0));
		for (uint Pos = 1; Pos < LQ; ++Pos)
			{
			uint Bin = Bins[Pos];
			if (Bin != CurrentBin)
				{
				Log("select tmp, resi %u-%u\n", Start+1, Pos);
				Log("color %s, tmp\n", GetColor(CurrentBin));
				Start = Pos;
				CurrentBin = Bin;
				}
			}
		Log("select tmp, resi %u-%u\n", Start+1, LQ);
		Log("color %s, tmp\n", GetColor(CurrentBin));
		}
	}
