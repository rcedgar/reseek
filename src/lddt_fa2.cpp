#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "flat_chain.h"
#include "chaq.h"

void trunc_label(string &Label);

static const uint M = 64;

double GetLDDT_muscle_fast(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs);

float flat_getlddt_muscle_some_floats2(
	const sid_t *distmxQ,
	const sid_t *distmxT,
	uint LQ, uint LT,
	const vector<uint32_t> &posQs,
	const vector<uint32_t> &posTs,
	const uint M);

static void rows2posvecs(
	const string &rowQ,
	const string &rowT, 
	uint LQ,
	uint LT,
	vector<uint> &posQs,
	vector<uint> &posTs)
	{
	posQs.clear();
	posTs.clear();

	const uint ncol = uint(rowQ.size());

	posQs.reserve(ncol);
	posTs.reserve(ncol);

	asserta(rowT.size() == ncol);
	uint posQ = 0;
	uint posT = 0;
	for (uint col = 0; col < ncol; ++col)
		{
		char q = rowQ[col];
		char t = rowT[col];

		if (isupper(q) && isupper(t))
			{
			asserta(posQ < LQ);
			asserta(posT < LT);
			posQs.push_back(posQ);
			posTs.push_back(posT);
			}

		if (isalpha(q)) ++posQ;
		if (isalpha(t)) ++posT;
		}
	asserta(posQ == LQ);
	asserta(posT == LT);
	}

void cmd_lddt_fa2()
	{
	asserta(optset_output);

	SeqDB fa2;
	fa2.FromFasta(g_Arg1, true);

	FILE *fout = CreateStdioFile(opt(output));

	vector<PDBChain *> OldChains;
	ReadChains(opt(input), OldChains);
	const uint ChainCount = SIZE(OldChains);
	unordered_map<string, uint> label2oldchainidx;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		string Label = OldChains[ChainIdx]->m_Label;
		trunc_label(Label);
		label2oldchainidx[Label] = ChainIdx;
		}

	vector<flat_chain_t *> flat_chains;
	unordered_map<string, uint> label2flatchainidx;
	read_flat_chains_idx_trunclabel(
		opt(input), flat_chains, label2flatchainidx);

	vector<uint16_t *> distmxs(ChainCount);

	for (unordered_map<string, uint>::const_iterator iter =
		label2flatchainidx.begin();
		iter != label2flatchainidx.end();
		++iter)
		{
		const string &label = iter->first;
		uint flatidx = iter->second;
		const flat_chain_t *flat_chain = flat_chains[flatidx];
		const uint L = flat_chain->get_length();
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(flat_chain->m_xyz->m_data, L, M, distmx);
		distmxs[flatidx] = distmx;
		}

	const uint nseq = fa2.GetSeqCount();
	asserta(nseq%2 == 0);
	const uint npair = nseq/2;
	for (uint pairidx = 0; pairidx < npair; ++pairidx)
		{
		ProgressStep(pairidx, npair, "Working");

		string labelQ = fa2.GetLabel(2*pairidx);
		string labelT = fa2.GetLabel(2*pairidx+1);
		trunc_label(labelQ);
		trunc_label(labelT);

		unordered_map<string, uint>::const_iterator flat_iterQ =
			label2flatchainidx.find(labelQ);
		unordered_map<string, uint>::const_iterator flat_iterT =
			label2flatchainidx.find(labelT);

		asserta(flat_iterQ != label2flatchainidx.end());
		asserta(flat_iterT != label2flatchainidx.end());

		uint flatidxQ = flat_iterQ->second;
		uint flatidxT = flat_iterT->second;

		const flat_chain_t *chainQ = flat_chains[flatidxQ];
		const flat_chain_t *chainT = flat_chains[flatidxT];
		const uint LQ = chainQ->get_length();
		const uint LT = chainT->get_length();

		const sid_t *distmxQ = distmxs[flatidxQ];
		const sid_t *distmxT = distmxs[flatidxT];

		unordered_map<string, uint>::const_iterator old_iterQ =
			label2oldchainidx.find(labelQ);
		unordered_map<string, uint>::const_iterator old_iterT =
			label2oldchainidx.find(labelT);

		asserta(old_iterQ != label2oldchainidx.end());
		asserta(old_iterT != label2oldchainidx.end());

		uint oldidxQ = old_iterQ->second;
		uint oldidxT = old_iterT->second;

		const PDBChain &OldChainQ = *OldChains[oldidxQ];
		const PDBChain &OldChainT = *OldChains[oldidxT];
		asserta(OldChainQ.GetSeqLength() == LQ);
		asserta(OldChainT.GetSeqLength() == LT);

		const string &rowQ = fa2.GetSeq(2*pairidx);
		const string &rowT = fa2.GetSeq(2*pairidx+1);

		const uint ncol = uint(rowQ.size());
		asserta(rowT.size() == ncol);
		vector<uint32_t> posQs;
		vector<uint32_t> posTs;
		rows2posvecs(rowQ, rowT, LQ, LT, posQs, posTs);

		float flat_lddt = flat_getlddt_muscle_some_floats2(
			distmxQ, distmxT, LQ, LT, posQs, posTs, M);

		double old_lddt = GetLDDT_muscle_fast(
			OldChainQ, OldChainT, posQs, posTs);

		fprintf(fout, "%s", labelQ.c_str());
		fprintf(fout, "\t%s", labelT.c_str());
		fprintf(fout, "\t%.4g", flat_lddt);
		fprintf(fout, "\t%.4g", old_lddt);
		fprintf(fout, "\n");
		}
	CloseStdioFile(fout);
	}