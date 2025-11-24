#include "myutils.h"
#include "daliscorer.h"
#include "featuretrainer2.h"
#include "sfasta.h"

static char GetLDDTChar(double LDDT)
	{
	asserta(LDDT >= 0 && LDDT <= 1);
	if (LDDT >= 0.9) return 'A';
	if (LDDT >= 0.8) return 'C';
	if (LDDT >= 0.7) return 'D';
	if (LDDT >= 0.6) return 'E';
	if (LDDT >= 0.5) return 'F';
	if (LDDT >= 0.4) return 'g';
	if (LDDT >= 0.3) return 'h';
	if (LDDT >= 0.2) return 'i';
	if (LDDT >= 0.1) return 'k';
	return 'l';
	}

void cmd_fa2_lddt()
	{
	const uint LDDT_W = 8;
	const double MIN_LDDT = 0.1;

	asserta(optset_db);
	const string ChainFN = opt(db);
	const string &Fa2FN = g_Arg1; // "src/2025-10_reseek_tune/big_fa2/tp.mints05.maxts25.fa2";

	FILE *fOut = CreateStdioFile(opt(output));
	FILE *fOut2 = CreateStdioFile(opt(output2));

	vector<PDBChain *> Chains;
	DALIScorer DS;
	map<string, uint> LabelToChainIdx;
	DS.LoadChains(ChainFN);
	const uint ChainCount = SIZE(DS.m_Chains);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		string Label = DS.m_Chains[ChainIdx]->m_Label;
		FeatureTrainer2::TruncLabel(Label);
		LabelToChainIdx[Label] = ChainIdx;
		}

	vector<string> Rows;
	vector<string> Labels;
	vector<string> FullLabels;
	vector<uint> ChainIdxs;

	SFasta SF;
	SF.Open(Fa2FN);
	SF.m_AllowGaps = true;

	bool Row1 = true;
	string FirstLabel = "";
	string FirstRow;
	for (;;)
		{
		const char* Seq = SF.GetNextSeq();
		if (Seq == 0)
			break;
		string Label = SF.GetLabel();
		FullLabels.push_back(Label);
		FeatureTrainer2::TruncLabel(Label);
		map<string, uint>::const_iterator iter = LabelToChainIdx.find(Label);
		if (iter == LabelToChainIdx.end())
			Die("Label not found >%s", Label.c_str());
		uint ChainIdx = iter->second;
		ChainIdxs.push_back(ChainIdx);
		const unsigned L = SF.GetSeqLength();
		asserta(L != 0);
		string Row;
		Row.reserve(L);
		for (uint i = 0; i < L; ++i)
			Row.push_back(Seq[i]);
		if (!Row1)
			asserta(SIZE(Row) == SIZE(Rows.back()));
		FeatureTrainer2::TruncLabel(Label, Label);
		Labels.push_back(Label);
		Rows.push_back(Row);
		Row1 = !Row1;
		}

	uint N = SIZE(Rows);
	asserta(N%2 == 0);
	asserta(SIZE(Labels) == N);
	asserta(SIZE(FullLabels) == N);
	asserta(SIZE(ChainIdxs) == N);

	SeqDB MSA;
	const uint PairCount = N/2;
	uint TotalColCount = 0;
	uint KeptColCount = 0;
	for (uint PairIdx = 0; PairIdx < PairCount; ++PairIdx)
		{
		ProgressStep(PairIdx, PairCount, "Processing");
		uint RowIdx1 = PairIdx*2;
		uint RowIdx2 = RowIdx1 + 1;
		const string &FullLabel1 = FullLabels[RowIdx1];
		const string &FullLabel2 = FullLabels[RowIdx2];
		const string &Label1 = Labels[RowIdx1];
		const string &Label2 = Labels[RowIdx2];
		string Row1 = Rows[RowIdx1];
		string Row2 = Rows[RowIdx2];
		const uint ColCount = SIZE(Row1);
		asserta(SIZE(Row2) == ColCount);

		MSA.m_Labels.clear();
		MSA.m_Seqs.clear();
		MSA.m_LabelToIndex.clear();

		MSA.m_IsAligned = true;
		MSA.m_IsNucleo = false;
		MSA.m_IsNucleoSet = true;
		MSA.m_ColCount = ColCount;

		MSA.m_Labels.push_back(Label1);
		MSA.m_Labels.push_back(Label2);

		MSA.m_Seqs.push_back(Row1);
		MSA.m_Seqs.push_back(Row2);

		DS.SetMSA("MSA", MSA, false, false);

		vector<double> LDDTs;
		string LDDTRow;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			double LDDT = DS.GetLDDTMuWCol(Col, LDDT_W);
			LDDTs.push_back(LDDT);
			char c = GetLDDTChar(LDDT);
			char c1 = Row1[Col];
			char c2 = Row2[Col];
			if (c1 == '-' || c2 == '-')
				c = '-';
			else if (islower(c1) || islower(c2) || c1 == '.' || c2 == '.')
				c = '.';
			LDDTRow += c;

			if (c1 == '-')
				Row2[Col] = tolower(c2);
			if (c2 == '-')
				Row1[Col] = tolower(c1);

			if (isupper(c1) && isupper(c2))
				{
				++TotalColCount;
				if (LDDT < MIN_LDDT)
					{
					Row1[Col] = tolower(c1);
					Row2[Col] = tolower(c2);
					}
				else
					++KeptColCount;
				}
			}

		if (fOut2 != 0)
			{
			fprintf(fOut2, ">%s\n", FullLabel1.c_str());
			fprintf(fOut2, "%s\n", Row1.c_str());
			fprintf(fOut2, ">%s\n", FullLabel2.c_str());
			fprintf(fOut2, "%s\n", Row2.c_str());
			fprintf(fOut2, "\n");
			}

		if (fOut != 0)
			{
			//fprintf(fOut, "\n");
			fprintf(fOut, ">%s\n", FullLabel1.c_str());
			fprintf(fOut, "%s\n", Row1.c_str());
			fprintf(fOut, ">%s\n", FullLabel2.c_str());
			fprintf(fOut, "%s\n", Row2.c_str());
			fprintf(fOut, ">%s\n", "LDDT");
			fprintf(fOut, "%s\n", LDDTRow.c_str());
			}
		}

	CloseStdioFile(fOut);
	CloseStdioFile(fOut2);

	ProgressLog("%u / %u columns > MIN_LDDT (%.1f%%)\n",
		KeptColCount,
		TotalColCount,
		GetPct(KeptColCount, TotalColCount));
	}
