#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "quarts.h"

void GetColPosVecs(const string &Row,
  vector<uint> &ColToPos, vector<uint> &PosToCol);

static char GetMeanSymbol(double d)
	{
	if (d < 4)
		return '@';
	if (d < 8)
		return '*';
	if (d < 16)
		return ':';
	return '.';
	}

static char GetStdDevSymbol(double d)
	{
	if (d < 2)
		return ' ';
	if (d < 4)
		return '.';
	if (d < 8)
		return ':';
	return '*';
	}

// Generate contact map profile from MSA
// Input labels and lengths must match to PDB chains, alphabet of MSA is ignored.
void cmd_msa2cmp()
	{
	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);
	const uint SeqCount = MSA.GetSeqCount();

	SeqDB Input;
	Input.FromFasta(g_Arg1, false);

	vector<PDBChain *> Chains;
	ReadChains(opt_input, Chains);

	FILE *fTsv = CreateStdioFile(opt_output);

	map<string, uint> LabelToChainIdx;
	const uint ChainCount = SIZE(Chains);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const string &Label = Chains[ChainIdx]->m_Label;
		LabelToChainIdx[Label] = ChainIdx;
		}

	vector<vector<uint> > ColToPosVec(SeqCount);
	vector<vector<uint> > PosToColVec(SeqCount);
	vector<uint> ChainIdxs;
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		const string &Label = Input.GetLabel(SeqIdx);
		map<string, uint>::const_iterator iter = LabelToChainIdx.find(Label);
		if (iter == LabelToChainIdx.end())
			Die("Label not found in chains >%s", Label.c_str());
		uint ChainIdx = iter->second;
		ChainIdxs.push_back(ChainIdx);
		const PDBChain &Chain = *Chains[ChainIdx];
		uint L = Input.GetSeqLength(SeqIdx);
		uint L2 = Chain.GetSeqLength();
		if (L != L2)
			Die("Lengths disagree %u, %u > %s", L, L2, Label.c_str());

		const string &Row = MSA.GetSeq(SeqIdx);
		GetColPosVecs(Row, ColToPosVec[SeqIdx], PosToColVec[SeqIdx]);
		}
	double MaxGapFract = 0.2;
	if (optset_maxgappct)
		MaxGapFract = opt_maxgappct/100.0;
	const uint MSAColCount = MSA.GetColCount();
	vector<uint> ProfColToMSACol;
	uint GappyCount = 0;
	vector<uint> MSAColToProfCol;
	for (uint MSACol = 0; MSACol < MSAColCount; ++MSACol)
		{
		uint GapCount = MSA.GetGapCount(MSACol);
		double GapFract = double(GapCount)/double(SeqCount);
		if (GapFract <= MaxGapFract)
			{
			uint ProfCol = SIZE(ProfColToMSACol);
			MSAColToProfCol.push_back(ProfCol);
			ProfColToMSACol.push_back(MSACol);
			}
		else
			{
			MSAColToProfCol.push_back(UINT_MAX);
			++GappyCount;
			}
		}

	const uint ProfColCount = SIZE(ProfColToMSACol);
	ProgressLog("%u chains, %u / %u prof cols (%.1f%%)\n",
	  SeqCount, ProfColCount, MSAColCount, GetPct(ProfColCount, MSAColCount));

	vector<vector<double> > MeanMx(ProfColCount);
	vector<vector<double> > StdDevMx(ProfColCount);
	for (uint i1 = 0; i1 < ProfColCount; ++i1)
		{
		MeanMx[i1].resize(ProfColCount, DBL_MAX);
		StdDevMx[i1].resize(ProfColCount, DBL_MAX);
		}

	for (uint i1 = 0; i1 < ProfColCount; ++i1)
		{
		uint Col1 = ProfColToMSACol[i1];
		for (uint i2 = i1+1; i2 < ProfColCount; ++i2)
			{
			uint Col2 = ProfColToMSACol[i2];

			vector<float> Dists;
			for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
				{
				uint Pos1 = ColToPosVec[SeqIdx][Col1];
				uint Pos2 = ColToPosVec[SeqIdx][Col2];
				if (Pos1 != UINT_MAX && Pos2 != UINT_MAX)
					{
					const PDBChain &Chain = *Chains[ChainIdxs[SeqIdx]];
					float d = (float) Chain.GetDist(Pos1, Pos2);
					Dists.push_back(d);
					}
				}

			QuartsFloat QF;
			GetQuartsFloat(Dists, QF);

			MeanMx[i1][i2] = QF.Avg;
			MeanMx[i2][i1] = QF.Avg;
			StdDevMx[i1][i2] = QF.StdDev;
			StdDevMx[i2][i1] = QF.StdDev;
			}
		}

	if (fTsv)
		{
		FILE *f = fTsv;
		fprintf(f, "%u\t%u\t%u\n", SeqCount, MSAColCount, ProfColCount);

	// Full MSA
		for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
			{
			uint ChainIdx = ChainIdxs[SeqIdx];
			const vector<uint> &ColToPos = ColToPosVec[SeqIdx];
			asserta(SIZE(ColToPos) == MSAColCount);
			const string &Label = Chains[ChainIdx]->m_Label;
			const string &Seq = Chains[ChainIdx]->m_Seq;
			fprintf(f, "%u\t%s\t", SeqIdx, Label.c_str());
			for (uint MSACol = 0; MSACol < MSAColCount; ++MSACol)
				{
				uint Pos = ColToPos[MSACol];
				if (Pos == UINT_MAX)
					fprintf(f, "-");
				else
					{
					asserta(Pos < SIZE(Seq));
					fprintf(f, "%c", Seq[Pos]);
					}
				}
			fprintf(f, "\n");
			}

	// Profile MSA
		for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
			{
			uint ChainIdx = ChainIdxs[SeqIdx];
			const vector<uint> &ColToPos = ColToPosVec[SeqIdx];
			asserta(SIZE(ColToPos) == MSAColCount);
			const string &Label = Chains[ChainIdx]->m_Label;
			const string &Seq = Chains[ChainIdx]->m_Seq;
			fprintf(f, "%u\t%s\t", SeqIdx, Label.c_str());
			for (uint ProfCol = 0; ProfCol < ProfColCount; ++ProfCol)
				{
				uint MSACol = ProfColToMSACol[ProfCol];
				uint Pos = ColToPos[MSACol];
				if (Pos == UINT_MAX)
					fprintf(f, "-");
				else
					{
					asserta(Pos < SIZE(Seq));
					fprintf(f, "%c", Seq[Pos]);
					}
				}
			fprintf(f, "\n");
			}

		for (uint Col1 = 0; Col1 < ProfColCount; ++Col1)
			{
			fprintf(f, "%u", Col1);
			string Row;
			for (uint Col2 = 0; Col2 < ProfColCount; ++Col2)
				{
				if (Col2 == Col1)
					fprintf(f, "\t*");
				else if (Col1 > Col2)
					fprintf(f, "\t%.3g", MeanMx[Col1][Col2]);
				else if (Col1 < Col2)
					fprintf(f, "\t%.3g", StdDevMx[Col1][Col2]);
				else
					asserta(false);
				}
			fprintf(f, "\n");
			}
		}
	CloseStdioFile(fTsv);
	}
