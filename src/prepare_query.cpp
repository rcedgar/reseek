#include "myutils.h"
#include "dssparams.h"
#include "chainreader2.h"
#include "xdpmem.h"
#include <set>

float ViterbiFastMem(XDPMem &Mem, const char *A, uint LA,
					 const char *B, uint LB, string &Path);

static double GetPctId(const string &Seq_i, const string &Seq_j)
	{
	const uint LA = SIZE(Seq_i);
	const uint LB = SIZE(Seq_j);
	if (LA == LB && Seq_i == Seq_j)
		return 100;
	const char *A = Seq_i.c_str();
	const char *B = Seq_j.c_str();
	string Path;
	XDPMem Mem;
	ViterbiFastMem(Mem, A, LA, B, LB, Path);
	const uint ColCount = SIZE(Path);
	uint PosA = 0;
	uint PosB = 0;
	uint Ids = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			if (Seq_i[PosA] == Seq_j[PosB])
				++Ids;
			++PosA;
			++PosB;
			break;
		case 'D':
			++PosA;
			break;
		case 'I':
			++PosB;
			break;
			}
		}
	double PctId = (100.0*Ids)/ColCount;
	return PctId;
	}

void cmd_prepare_query()
	{
	DSSParams Params;
	Params.SetDSSParams(DM_AlwaysFast);
	asserta(optset_bca);
	BCAData BCA;
	if (optset_bca)
		BCA.Create(opt(bca));


	FILE *fOut = CreateStdioFile(opt(output));
	vector<PDBChain *> InputChains;
	ReadChains(g_Arg1, InputChains);
	const uint InputChainCount = SIZE(InputChains);
	const double MinPctId = 90;
	const uint MinLen = (optset_minchainlength ? opt(minchainlength) : 30);
	const uint MaxChains = (optset_minchainlength ? opt(n) : 4);
	vector<PDBChain *> OutputChains;

	set<uint> DeletedChainIdxs;
	uint TooShort = 0;
	uint NrQueries = 0;
	for (uint i = 0; i < InputChainCount; ++i)
		{
		if (DeletedChainIdxs.find(i) != DeletedChainIdxs.end())
			continue;
		const PDBChain &Chain_i = *InputChains[i];
		const string &Label_i = Chain_i.m_Label;
		const string &Seq_i = Chain_i.m_Seq;
		uint L = SIZE(Seq_i);
		Pf(fOut, "%u\t%s\t%u", i, Label_i.c_str(), L);
		if (L < MinLen)
			{
			Pf(fOut, "\tshort\n", L);
			++TooShort;
			continue;
			}

		bool del = false;
		if (NrQueries >= MaxChains)
			{
			del = true;
			Pf(fOut, "\ttoomany\n");
			continue;
			}
		for (uint j = 0; j < i; ++j)
			{
			if (DeletedChainIdxs.find(j) != DeletedChainIdxs.end())
				continue;
			const PDBChain &Chain_j = *InputChains[j];
			const string &Label_j = Chain_j.m_Label;
			const string &Seq_j = Chain_j.m_Seq;
			if (SIZE(Seq_j) < MinLen)
				continue;
			double PctId = GetPctId(Seq_i, Seq_j);
			if (PctId >= MinPctId)
				{
				Pf(fOut, "\t%.1f%%%u\n", PctId, j);
				DeletedChainIdxs.insert(i);
				del = true;
				break;
				}
			}
		if (!del)
			{
			OutputChains.push_back(InputChains[i]);
			++NrQueries;
			Pf(fOut, "\tquery\n");
			}
		}
	uint OutputChainCount = SIZE(OutputChains);
	for (uint i = 0; i < OutputChainCount; ++i)
		BCA.WriteChain(*OutputChains[i]);
	BCA.Close();
	CloseStdioFile(fOut);
	}
