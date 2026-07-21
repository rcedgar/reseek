#include "myutils.h"
#include "flat_chain.h"
#include "flat_distmx.h"
#include "flat_helpers.h"
#include "chaq.h"
#include "cigar.h"

/***
Input tsv:
query+target+qlo+qhi+tlo+thi+cigar+dali
    0      1   2   3   4   5     6    7

Output tsv:
query+target+dali2+dali_greedy
***/

static uint GetChainIdx(
	const unordered_map<string, uint> &LabelToIdx,
	const string &Label,
	uint LineNr)
	{
	auto p = LabelToIdx.find(Label);
	if (p == LabelToIdx.end())
		Die("Line %u, chain label not found: '%s'", LineNr, Label.c_str());
	return p->second;
	}

static void CigarToMatchedPositions(
	const string &Cigar,
	uint QLo, uint QHi, uint LQ,
	uint TLo, uint THi, uint LT,
	uint LineNr,
	vector<uint> &PosQs,
	vector<uint> &PosTs)
	{
	if (QLo == 0 || QHi == 0 || TLo == 0 || THi == 0)
		Die("Line %u, coordinates must be one-based", LineNr);
	if (QLo > QHi || QHi > LQ)
		Die("Line %u, invalid query coordinates %u-%u (length %u)",
			LineNr, QLo, QHi, LQ);
	if (TLo > THi || THi > LT)
		Die("Line %u, invalid target coordinates %u-%u (length %u)",
			LineNr, TLo, THi, LT);

	string Ops;
	vector<uint> Lengths;
	CIGARGetOps(Cigar, Ops, Lengths);
	if (Ops.empty())
		Die("Line %u, empty CIGAR", LineNr);

	PosQs.clear();
	PosTs.clear();
	uint QPos = QLo - 1;
	uint TPos = TLo - 1;
	bool Started = false;
	bool GotS = false;
	bool GotT = false;
	const uint OpCount = SIZE(Ops);
	for (uint OpIdx = 0; OpIdx < OpCount; ++OpIdx)
		{
		const char Op = Ops[OpIdx];
		const uint Length = Lengths[OpIdx];
		if (Op == 'S')
			{
			if (Started || GotS || Length != QLo - 1)
				Die("Line %u, invalid query soft clip in CIGAR '%s'",
					LineNr, Cigar.c_str());
			GotS = true;
			continue;
			}
		if (Op == 'T')
			{
			if (Started || GotT || Length != TLo - 1)
				Die("Line %u, invalid target soft clip in CIGAR '%s'",
					LineNr, Cigar.c_str());
			GotT = true;
			continue;
			}

		Started = true;
		switch (Op)
			{
		case 'M':
		case '=':
		case 'X':
			for (uint i = 0; i < Length; ++i)
				{
				if (QPos >= LQ || TPos >= LT)
					Die("Line %u, CIGAR exceeds chain length", LineNr);
				PosQs.push_back(QPos++);
				PosTs.push_back(TPos++);
				}
			break;

		case 'I':
			if (Length > LQ - QPos)
				Die("Line %u, CIGAR exceeds query length", LineNr);
			QPos += Length;
			break;

		case 'D':
			if (Length > LT - TPos)
				Die("Line %u, CIGAR exceeds target length", LineNr);
			TPos += Length;
			break;

		default:
			Die("Line %u, unsupported operation '%c' in CIGAR '%s'",
				LineNr, Op, Cigar.c_str());
			}
		}

	if (QPos != QHi || TPos != THi)
		Die("Line %u, CIGAR endpoints are qhi=%u, thi=%u; expected %u, %u",
			LineNr, QPos, TPos, QHi, THi);
	}

void cmd_dali_roundtrip()
	{
	asserta(optset_input);

	vector<flat_chain_t *> Chains;
	read_flat_chains(g_Arg1, Chains);
	const uint ChainCount = SIZE(Chains);
	unordered_map<string, uint> LabelToIdx;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const string &Label = Chains[ChainIdx]->m_label;
		if (LabelToIdx.find(Label) != LabelToIdx.end())
			Die("Duplicate chain label '%s'", Label.c_str());
		LabelToIdx[Label] = ChainIdx;
		}

	vector<sid_t *> DistMxs(ChainCount, 0);
	FILE *fIn = OpenStdioFile(opt(input));
	FILE *fOut = CreateStdioFile(opt(output));
	string Line;
	vector<string> Fields;
	vector<uint> PosQs;
	vector<uint> PosTs;
	vector<uint8_t> DaliTermsScratch;
	vector<uint8_t> DaliWorkScratch;
	vector<uint> RetainedCols;
	uint LineNr = 0;
	while (ReadLineStdioFile(fIn, Line))
		{
		++LineNr;
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 8)
			Die("Line %u, expected 8 tab-separated fields, got %u",
				LineNr, SIZE(Fields));

		const uint QIdx = GetChainIdx(LabelToIdx, Fields[0], LineNr);
		const uint TIdx = GetChainIdx(LabelToIdx, Fields[1], LineNr);
		const flat_chain_t *Q = Chains[QIdx];
		const flat_chain_t *T = Chains[TIdx];
		const uint QLo = StrToUint(Fields[2]);
		const uint QHi = StrToUint(Fields[3]);
		const uint TLo = StrToUint(Fields[4]);
		const uint THi = StrToUint(Fields[5]);
		CigarToMatchedPositions(
			Fields[6],
			QLo, QHi, Q->m_L,
			TLo, THi, T->m_L,
			LineNr, PosQs, PosTs);

		if (DistMxs[QIdx] == 0)
			{
			DistMxs[QIdx] = myalloc(
				sid_t, Q->m_L*flat_params::m_distmx_bandwidth);
			chaq::fill_distmx(Q, DistMxs[QIdx]);
			}
		if (DistMxs[TIdx] == 0)
			{
			DistMxs[TIdx] = myalloc(
				sid_t, T->m_L*flat_params::m_distmx_bandwidth);
			chaq::fill_distmx(T, DistMxs[TIdx]);
			}

		const float Dali2 = flat_get_dali4(
			PosQs.data(), Q->m_L,
			PosTs.data(), T->m_L, SIZE(PosQs),
			DistMxs[QIdx], DistMxs[TIdx]);
		const uint NMatch = SIZE(PosQs);
		const size_t TermsBytes = dali_greedy_terms_bytes(NMatch);
		const size_t WorkBytes = dali_greedy_work_bytes(NMatch);
		if (DaliTermsScratch.size() < TermsBytes)
			DaliTermsScratch.resize(TermsBytes);
		if (DaliWorkScratch.size() < WorkBytes)
			DaliWorkScratch.resize(WorkBytes);
		if (RetainedCols.size() < NMatch)
			RetainedCols.resize(NMatch);
		uint NRetained = 0;
		const float DaliGreedy = dali_greedy(
			PosQs.data(), Q->m_L,
			PosTs.data(), T->m_L, NMatch,
			DistMxs[QIdx], DistMxs[TIdx],
			DaliTermsScratch.data(), DaliTermsScratch.size(),
			DaliWorkScratch.data(), DaliWorkScratch.size(),
			RetainedCols.data(), SIZE(RetainedCols), NRetained);
		// fprintf(fOut, "%s\t%.3g\n", Line.c_str(), Dali2);
		if (fOut != 0)
			{
			fprintf(fOut, "%s", Fields[0].c_str());
			fprintf(fOut, "\t%s", Fields[1].c_str());
			fprintf(fOut, "\t%.3g", DaliGreedy);
			fprintf(fOut, "\t%.3g", Dali2);
			fprintf(fOut, "\n");
			}
		}

	CloseStdioFile(fIn);
	CloseStdioFile(fOut);
	//for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
	//	{
	//	myfree(DistMxs[ChainIdx]);
	//	delete Chains[ChainIdx];
	//	}
	}
