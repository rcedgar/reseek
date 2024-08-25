#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"
#include "alpha.h"

char GetSSConsSymbol(double SSCons);

static const char *GetColor(uint Bin)
	{
	asserta(Bin >= 0 && Bin <= 10);
	static const char *Colors[11] =
		{ "br0","br1","br2","br3","br4","br5","br6","br7","br8","br9","br10" };
	return Colors[Bin];
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

static void GetSSVec(const vector<PDBChain *> &Chains,
  vector<string> &SSVec)
	{
	const uint N = SIZE(Chains);
	for (uint i = 0 ; i < N; ++i)
		{
		string SS;
		PDBChain &Chain = *Chains[i];
		Chain.GetSS(SS);
		SSVec.push_back(SS);
		}
	}

static void GetSSMSA(const SeqDB &MSA, const vector<uint> &SeqIdxToChainIdx, 
  const vector<string> &SSVec, vector<string> &SSMSA)
	{
	const uint SeqCount = MSA.GetSeqCount();
	const uint ColCount = MSA.GetColCount();
	const uint ChainCount = SIZE(SSVec);
	asserta(SIZE(SeqIdxToChainIdx) == SeqCount);
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		const string &Row = MSA.GetSeq(SeqIdx);
		asserta(SIZE(Row) == ColCount);
		uint ChainIdx = SeqIdxToChainIdx[SeqIdx];
		if (ChainIdx == UINT_MAX)
			continue;
		asserta(ChainIdx < ChainCount);
		const string &SS = SSVec[ChainIdx];
		const uint L = SIZE(SS);
		string SSRow;
		uint Pos = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char c = Row[Col];
			if (isgap(c))
				SSRow += c;
			else
				{
				asserta(Pos < L);
				SSRow += SS[Pos++];
				}
			}
		asserta(Pos == L);
		SSMSA.push_back(SSRow);
		}
	}

static char GetConsChar3(const vector<string> &SSMSA, uint Col)
	{
	const uint SeqCount = SIZE(SSMSA);
	vector<uint> Counts(4);
	for (uint SeqIdx = 0; SeqIdx  < SeqCount; ++SeqIdx)
		{
		char c = SSMSA[SeqIdx][Col];
		switch (c)
			{
		case 'h': ++Counts[0]; break;
		case 's': ++Counts[1]; break;
		case 't': ++Counts[2]; break;
		case '~': ++Counts[3]; break;
		case '-': break;
		default: asserta(false);
			}
		}
	uint MaxCount = 0;
	char MaxChar = '-';
	for (uint i = 0; i < 4; ++i)
		{
		uint n = Counts[i];
		if (n > MaxCount)
			{
			MaxCount = n;
			MaxChar = "hst~"[i];
			}
		}
	return MaxChar;
	}

static const char *GetCons3Color(char c3)
	{
	switch (c3)
		{
	case 'h': return "0,150,20";
	case 's': return "150,0,50";
	case 't': return "250,150,0";
	case '~': return "150,150,150";
	case '-': return "255,255,255";
	default: asserta(false);
		}
	return "0,0,0";
	}

static void SmoothS3(string &S3)
	{
	string s;
	const uint ColCount = SIZE(S3);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c3 = S3[Col];
		char Prev = 0;
		char Next = 0;
		for (int i = (int) Col - 1; i >= 0; --i)
			{
			if (!isgap(S3[i]))
				{
				Prev = S3[i];
				break;
				}
			}
		for (int i = (int) Col + 1; i < int(ColCount); ++i)
			{
			if (!isgap(S3[i]))
				{
				Next = S3[i];
				break;
				}
			}
		if (c3 != 's' && c3 != 'h')
			continue;
		if (Prev == 0 || Next == 0)
			continue;
		if (Prev != c3 && Next != c3)
			{
			if (Prev == Next)
				S3[Col] = Prev;
			else
				S3[Col] = '~';
			}
		}
	}

// -msta_sscons MSA -input STRUCTS
void cmd_msta_sscons()
	{
	asserta(optset_input);
	if (optset_sscons_pymol && !optset_label)
		Die("-sscons_pymol requires -label");

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
		SSConsVec.push_back(SSCons);
		}

	if (optset_sscons_jalview)
		{
		vector<string> SSVec;
		GetSSVec(DS.m_Chains, SSVec);
		vector<string> SSMSA;
		GetSSMSA(MSA, DS.m_SeqIdxToChainIdx, SSVec, SSMSA);
		FILE *f = CreateStdioFile(opt_sscons_jalview);
		fprintf(f, "JALVIEW_ANNOTATION\n");
		fprintf(f, "BAR_GRAPH\tSSCons\t");
		string S3;
		for (uint Col = 0; Col < ColCount; ++Col)
			S3 += GetConsChar3(SSMSA, Col);
		SmoothS3(S3);
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char c3 = S3[Col];
			const char *Color = GetCons3Color(c3);
			double SSCons = DS.GetSSConsCol(Col, w);
			if (Col > 0)
				fprintf(f, "|");
			fprintf(f, "%.3f[%s]", SSCons, Color);
			}
		fprintf(f, "\n");
		CloseStdioFile(f);
		}

	vector<double> Thresholds;
	Thresholds.push_back(0.1);
	Thresholds.push_back(0.2);
	Thresholds.push_back(0.3);
	Thresholds.push_back(0.4);
	Thresholds.push_back(0.5);
	Thresholds.push_back(0.6);
	Thresholds.push_back(0.7);
	Thresholds.push_back(0.8);
	Thresholds.push_back(0.9);

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
		if (optset_sscons_pymol)
			{
			FILE *f = CreateStdioFile(opt_sscons_pymol);
			uint Start = 0;
			uint CurrentBin = Bins[0];
			fprintf(f, "select tmp, all\n");
			fprintf(f, "color %s, tmp\n", GetColor(0));
			for (uint Pos = 1; Pos < LQ; ++Pos)
				{
				uint Bin = Bins[Pos];
				if (Bin != CurrentBin)
					{
					fprintf(f, "select tmp, resi %u-%u\n", Start+1, Pos);
					fprintf(f, "color %s, tmp\n", GetColor(CurrentBin));
					Start = Pos;
					CurrentBin = Bin;
					}
				}
			fprintf(f, "select tmp, resi %u-%u\n", Start+1, LQ);
			fprintf(f, "color %s, tmp\n", GetColor(CurrentBin));
			fprintf(f, "select none\n");
			}
		}
	}
