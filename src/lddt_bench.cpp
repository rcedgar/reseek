#include "myutils.h"
#include "seqdb.h"
#include "pdbchain.h"

#define COMPARE	0
#define	FAST	1

double GetLDDT_mu(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  bool DaliScorerCompatible);
double GetLDDT_mu_fast(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs);

void cmd_lddt_bench()
	{
	asserta(optset_input);
	const bool MissingTestSeqOk = opt(missingtestseqok);

	string Name;
	GetStemName(g_Arg1, Name);

	SeqDB MSA;
	MSA.FromFasta(g_Arg1, true);
	const uint SeqCount = MSA.GetSeqCount();

	vector<PDBChain *> Chains;
	ReadChains(opt(input), Chains);

	const uint ChainCount = SIZE(Chains);
	map<string, uint> LabelToChainIdx;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const string &Label = Chains[ChainIdx]->m_Label;
		LabelToChainIdx[Label] = ChainIdx;
		}

	vector<uint> ChainIdxs;
	const uint ColCount = MSA.GetColCount();
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		const string &Label = MSA.GetLabel(SeqIdx);
		map<string, uint>::const_iterator iter = LabelToChainIdx.find(Label);
		asserta(iter != LabelToChainIdx.end());
		uint ChainIdx = iter->second;
		ChainIdxs.push_back(ChainIdx);
		}

#if COMPARE
	const uint NITERS = 1;
#else
	const uint NITERS = 20;
#endif
	double SumLDDT = 0;
	uint PairCount = 0;
	for (uint Iter = 0; Iter < NITERS; ++Iter)
		{
		ProgressStep(Iter, NITERS, "Working");
		for (uint SeqIdx1 = 0; SeqIdx1 < SeqCount; ++SeqIdx1)
			{
			const string &Row1 = MSA.GetSeq(SeqIdx1);
			const uint ChainIdx1 = ChainIdxs[SeqIdx1];
			const PDBChain &Chain1 = *Chains[ChainIdx1];
			const char *Label1 = MSA.GetLabel(SeqIdx1).c_str();
			for (uint SeqIdx2 = SeqIdx1 + 1; SeqIdx2 < SeqCount; ++SeqIdx2)
				{
				const string &Row2 = MSA.GetSeq(SeqIdx2);
				const uint ChainIdx2 = ChainIdxs[SeqIdx2];
				const PDBChain &Chain2 = *Chains[ChainIdx2];
				const char *Label2 = MSA.GetLabel(SeqIdx2).c_str();
				if (ChainIdx1 == UINT_MAX || ChainIdx2 == UINT_MAX)
					Die("structure_not_found %s %s\n", Label1, Label2);
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

#if FAST || COMPARE
				double PairLDDT_fast = GetLDDT_mu_fast(Chain1, Chain2, Pos1s, Pos2s);
#endif
#if !FAST || COMPARE
				double PairLDDT = GetLDDT_mu(Chain1, Chain2, Pos1s, Pos2s, false);
#endif
#if COMPARE
				//ProgressLog("%.4f %.4f\n", PairLDDT, PairLDDT_fast);
				asserta(feq(PairLDDT, PairLDDT_fast));
#endif
				if (Iter == 0)
					{
					++PairCount;
#if FAST
					SumLDDT += PairLDDT_fast;
#else
					SumLDDT += PairLDDT;
#endif
					}
				}
			}
		}

	double LDDT = 0;
	if (PairCount > 0)
		LDDT = SumLDDT/PairCount;
	ProgressLog("LDDT=%.4f MSA=%s\n", LDDT, Name.c_str());
	}
