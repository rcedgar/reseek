#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "daliscorer.h"
#include "alpha.h"

char GetSSConsSymbol(double SSCons)
	{
	if (SSCons == 0)
		return ' ';
	if (SSCons < 0.2)
		return '.';
	if (SSCons < 0.5)
		return ':';
	if (SSCons < 0.75)
		return '|';
	return '@';
	}

double DALIScorer::GetLDDTScoreWindow(const vector<vector<double> > &DistMx1,
  const vector<vector<double> > &DistMx2, uint w) const
	{
	const uint nr_thresholds = SIZE(m_LDDT_thresholds);
	const uint n = 2*w + 1;
	asserta(SIZE(DistMx1) == n);
	asserta(SIZE(DistMx2) == n);

	double total = 0;
	uint nr_cols_considered = 0;
	int iw = int(w);
	for (int ii = -iw; ii <= iw; ++ii)
		{
		uint i = uint(ii+iw);
		const vector<double> &DistMxRow1i = DistMx1[i];
		const vector<double> &DistMxRow2i = DistMx2[i];
		for (int jj = i + 2; jj <= iw; ++jj)
			{
			uint j = uint(jj+iw);

			++nr_cols_considered;
			uint nr_considered = 0;
			uint nr_preserved = 0;
			double d1 = DistMxRow1i[j];
			double d2 = DistMxRow2i[j];
			if (d1 == DBL_MAX || d2 == DBL_MAX)
				continue;

			if (d1 > m_LDDT_R0)
				continue;

			for (uint k = 0; k < nr_thresholds; ++k)
				{
				double t = m_LDDT_thresholds[k];
				nr_considered += 1;
				double diff = abs(d1 - d2);
				if (diff <= t)
					nr_preserved += 1;
				}
			double score = 0;
			if (nr_considered > 0)
				score = double(nr_preserved)/nr_considered;
			total += score;
			}
		}

	if (nr_cols_considered == 0)
		return 0;
	double avg = total/nr_cols_considered;
	return avg;
	}

void DALIScorer::GetDistMxWindow(uint SeqIdx, uint Pos, uint w,
  vector<vector<double> > &Mx) const
	{
	asserta(SeqIdx < SIZE(m_SeqIdxToChainIdx));
	uint ChainId = m_SeqIdxToChainIdx[SeqIdx];
	const vector<vector<double> > &DistMx = m_DistMxVec[ChainId];
	const uint L = SIZE(DistMx);

	uint n = 2*w + 1;
	Mx.resize(n);
	for (uint i = 0; i < n; ++i)
		Mx[i].resize(n, DBL_MAX);

	for (uint i = 0; i < n; ++i)
		{
		Mx[i][i] = 0;
		int Posi = int(Pos) - int(w) + i;
		if (Posi < 0 || Posi >= int(L))
			continue;
		asserta(Posi >= int(Pos) - int(w) &&  Posi <= int(Pos) + int(w));
		for (uint j = i+1; j < n; ++j)
			{
			int Posj = int(Pos) - int(w) + j;
			if (Posj < 0 || Posj >= int(L))
				continue;
			asserta(Posj >= int(Pos) - int(w) &&  Posi <= int(Pos) + int(w));
			double d = DistMx[Posi][Posj];
			Mx[i][j] = d;
			Mx[j][i] = d;
			}
		}
	}

double DALIScorer::GetSSCons1(uint QuerySeqIdx, uint Col, uint w) const
	{
	asserta(QuerySeqIdx < SIZE(m_ColToPosVec));
	asserta(Col < SIZE(m_ColToPosVec[QuerySeqIdx]));
	uint QueryPos = m_ColToPosVec[QuerySeqIdx][Col];
	vector<vector<double> > QueryMx;
	GetDistMxWindow(QuerySeqIdx, QueryPos, w, QueryMx);
	asserta(QuerySeqIdx < SIZE(m_SeqIdxToChainIdx));
	uint QueryChainIdx = m_SeqIdxToChainIdx[QuerySeqIdx];

	uint nr_seqs = GetSeqCount();
	uint nr_cols = GetColCount();
	double sum = 0;
	uint n = 0;
	for (uint seq_idxi = 0; seq_idxi < nr_seqs; ++seq_idxi)
		{
		if (seq_idxi == QuerySeqIdx)
			continue;
		uint posi = m_ColToPosVec[seq_idxi][Col];
		if (posi == UINT_MAX)
			continue;

		vector<vector<double> > Mx;
		GetDistMxWindow(seq_idxi, posi, w, Mx);
		double score = GetLDDTScoreWindow(QueryMx, Mx, w);
		if (score == DBL_MAX)
			continue;
		sum += score;
		++n;
		}
	if (n == 0)
		return 0;
	double avg = sum/n;
	return avg;
	}

// -msta_sscons MSA -input STRUCTS -label QUERY
void cmd_msta_sscons1()
	{
	asserta(optset_input);
	asserta(optset_label);

	uint w = 2;
	if (optset_window)
		w = opt_window;

	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);
	MSA.SetLabelToIndex();
	const uint ColCount = MSA.GetColCount();
	const string &Label = string(opt_label);
	const uint QuerySeqIdx = MSA.GetSeqIndex(Label);
	const string &QueryRow = MSA.GetSeq(QuerySeqIdx);
	string QuerySeq;
	for (uint i = 0; i < SIZE(QueryRow); ++i)
		{
		char c = QueryRow[i];
		if (!isgap(c))
			QuerySeq += toupper(c);
		}

	//FILE *fOut = CreateStdioFile(opt_output);

	string Name;
	GetStemName(g_Arg1, Name);

	const bool DoCore = false;
	const bool MissingTestSeqOk = true;

	DALIScorer DS;
	DS.LoadChains(opt_input);
	bool Ok = DS.SetMSA(Name, MSA, DoCore, MissingTestSeqOk);
	if (!Ok)
		Die("SetMSA failed");

	const uint QueryChainIdx = DS.GetChainIdx(QuerySeq);

	const uint LQ = SIZE(QuerySeq);
	uint PosQ = 0;
	vector<double> SSConsVec;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char q = QueryRow[Col];
		if (isgap(q))
			continue;
		double SSCons = DS.GetSSCons1(QuerySeqIdx, PosQ, w);
		Log("%u", PosQ);
		Log("  %c", q);
		if (SSCons == DBL_MAX)
			Log("  .");
		else
			Log("  %.4f", SSCons);
		Log("\n");
		SSConsVec.push_back(SSCons);
		++PosQ;
		}

	Log("%s\n", QuerySeq.c_str());
	for (uint PosQ = 0; PosQ < LQ; ++PosQ)
		{
		double SSCons = SSConsVec[PosQ];
		char c = GetSSConsSymbol(SSCons);
		Log("%c", c);
		}
	Log("\n");
	}
