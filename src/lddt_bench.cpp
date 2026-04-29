#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"
#include "flat_chain.h"
#include "chaq.h"

static const uint M = 64;

#define COMPARE	0
#define	FAST	1

double GetLDDT_muscle(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  bool DaliScorerCompatible);
double GetLDDT_muscle_fast(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs);

float flat_getlddt_muscle_some_floats(
	const uint32_t *posQs,
	const uint32_t LQ,
	const uint32_t *posTs,
	const uint32_t LT,
	const uint ncol,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	uint32_t *nr_considered_vec,
	uint32_t *nr_preserved_vec);

void cmd_lddt_bench()
	{
	const bool MissingTestSeqOk = opt(missingtestseqok);

	string Name;
	GetStemName(g_Arg1, Name);

	FILE *fout = CreateStdioFile(opt(output));

	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);
	const uint SeqCount = MSA.GetSeqCount();

	vector<PDBChain *> Chains;
	ReadChains(opt(input), Chains);

	vector<flat_chain_t *> flat_chains;
	unordered_map<string, uint> label2flatchainidx;
	read_flat_chains_idx(opt(input), flat_chains, label2flatchainidx);

	const uint ChainCount = SIZE(Chains);
	map<string, uint> LabelToChainIdx;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const string &Label = Chains[ChainIdx]->m_Label;
		LabelToChainIdx[Label] = ChainIdx;
		}

	vector<uint> ChainIdxs;
	vector<uint> flat_chainidxs;
	vector<uint16_t *> distmxs;
	const uint ColCount = MSA.GetColCount();
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		const string &Label = MSA.GetLabel(SeqIdx);
		map<string, uint>::const_iterator iter = LabelToChainIdx.find(Label);
		asserta(iter != LabelToChainIdx.end());
		uint ChainIdx = iter->second;
		ChainIdxs.push_back(ChainIdx);

		unordered_map<string, uint>::const_iterator flat_iter =
			label2flatchainidx.find(Label);
		asserta(flat_iter != label2flatchainidx.end());
		uint flat_chainidx = flat_iter->second;
		flat_chainidxs.push_back(flat_chainidx);

		const flat_chain_t *flat_chain = flat_chains[flat_chainidx];
		const uint L = flat_chain->get_length();
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(flat_chain->m_xyz->m_data, L, distmx);
		distmxs.push_back(distmx);
		}

	double SumLDDT = 0;
	uint PairCount = 0;
	for (uint SeqIdx1 = 0; SeqIdx1 < SeqCount; ++SeqIdx1)
		{
		const string &Row1 = MSA.GetSeq(SeqIdx1);

		const uint ChainIdx1 = ChainIdxs[SeqIdx1];
		const PDBChain &Chain1 = *Chains[ChainIdx1];

		const uint flat_chainidx1 = flat_chainidxs[SeqIdx1];
		const flat_chain_t &flat_chain1 = *flat_chains[flat_chainidx1];
		asserta(flat_chain1.m_label == Chain1.m_Label);
		const sid_t *distmx1 = distmxs[SeqIdx1];

		const char *Label1 = MSA.GetLabel(SeqIdx1).c_str();
		const uint L1 = flat_chain1.get_length();

		for (uint SeqIdx2 = SeqIdx1 + 1; SeqIdx2 < SeqCount; ++SeqIdx2)
			{
			const string &Row2 = MSA.GetSeq(SeqIdx2);
			const uint ChainIdx2 = ChainIdxs[SeqIdx2];
			const PDBChain &Chain2 = *Chains[ChainIdx2];

			const uint flat_chainidx2 = flat_chainidxs[SeqIdx2];
			const flat_chain_t &flat_chain2 = *flat_chains[ChainIdx2];
			asserta(flat_chain2.m_label == Chain2.m_Label);

			const sid_t *distmx2 = distmxs[SeqIdx2];

			const char *Label2 = MSA.GetLabel(SeqIdx2).c_str();
			const uint L2 = flat_chain2.get_length();

			vector<uint> Pos1s;
			vector<uint> Pos2s;
			uint Pos1 = 0;
			uint Pos2 = 0;
			for (uint Col = 0; Col < ColCount; ++Col)
				{
				char c1 = Row1[Col];
				char c2 = Row2[Col];
				if (c1 != '-' && c2 != '-')
					{
					Pos1s.push_back(Pos1);
					Pos2s.push_back(Pos2);
					}
				if (c1 != '-')
					Pos1++;
				if (c2 != '-')
					Pos2++;
				}

			const uint ncol = uint(Pos1s.size());
			asserta(Pos2s.size() == ncol);
			uint32_t *nr_considered_vec = myalloc(uint32_t, ColCount);
			uint32_t *nr_preserved_vec = myalloc(uint32_t, ColCount);
			double PairLDDT_flat = flat_getlddt_muscle_some_floats(
				Pos1s.data(), L1,
				Pos2s.data(), L2,
				ncol,
				distmx1,
				distmx2,
				nr_considered_vec,
				nr_preserved_vec);
			myfree(nr_considered_vec);
			myfree(nr_preserved_vec);

#if FAST || COMPARE
			double PairLDDT_fast = GetLDDT_muscle_fast(Chain1, Chain2, Pos1s, Pos2s);
#endif
#if !FAST || COMPARE
			double PairLDDT = GetLDDT_muscle(Chain1, Chain2, Pos1s, Pos2s, false);
#endif
#if COMPARE
			//ProgressLog("%.4f %.4f\n", PairLDDT, PairLDDT_fast);
			asserta(feq(PairLDDT, PairLDDT_fast));
#endif
			if (fout != 0)
				{
				fprintf(fout, "%s", Label1);
				fprintf(fout, "\t%s", Label2);
				fprintf(fout, "\t%.4f", PairLDDT_fast);
				fprintf(fout, "\t%.4f", PairLDDT_flat);
				fprintf(fout, "\n");
				}
			++PairCount;
			if (optset_iters && PairCount >= opt(iters))
				break;
			}
		}
	}
