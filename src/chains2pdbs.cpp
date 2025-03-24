#include "myutils.h"
#include "pdbchain.h"

void GetThreeFromOne(char aa, string &AAA);

void PDBChain::ToPDB(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	ToPDB(f);
	CloseStdioFile(f);
	}

void PDBChain::ToPDB(FILE *f, bool TruncateAtZ) const
	{
	string ChainStr = "A";
	if (m_HasChainStr)
		ChainStr = m_ChainStr;
	const uint L = GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		char aa = m_Seq[i];
		string sAAA;
		GetThreeFromOne(aa, sAAA);
		const char *AAA = sAAA.c_str();

		fprintf(f, "ATOM  ");			//  1 -  6        Record name   "ATOM  "
		fprintf(f, "%5u", i+1);			//  7 - 11        Integer       serial       Atom serial number.
		fprintf(f, " ");				// 12
		fprintf(f, " CA ");				// 13 - 16        Atom          name         Atom name.
		fprintf(f, " ");				// 17             Character     altLoc       Alternate location indicator.
		fprintf(f, "%3.3s", AAA);		// 18 - 20        Residue name  resName      Residue name.
		fprintf(f, " ");				// 21
		fprintf(f, "%c", ChainStr[0]);	// 22             Character     chainID      Chain identifier.
		fprintf(f, "%4u", i+1);			// 23 - 26        Integer       resSeq       Residue sequence number.
		fprintf(f, " ");				// 27             AChar         iCode        Code for insertion of residues.
		fprintf(f, "   ");				// 28 - 30
		fprintf(f, "%8.3f", m_Xs[i]);	// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		fprintf(f, "%8.3f", m_Ys[i]);	// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		fprintf(f, "%8.3f", m_Zs[i]);	// 47 - 57        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		if (!TruncateAtZ)
			{
			fprintf(f, "%6.2f", 1.0);		// 55 - 60        Real(6.2)     occupancy    Occupancy.
			fprintf(f, "%6.2f", 0.0);		// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
			fprintf(f, "          ");		// 67 - 76
			fprintf(f, " C");				// 77 - 78        LString(2)    element      Element symbol, right-justified.
			fprintf(f, "  ");				// 79 - 80        LString(2)    charge       Charge on the atom.
			}

		fprintf(f, "\n");
		}
	}

void cmd_chains2pdbs()
	{
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	const uint N = SIZE(Chains);
	string fn;
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Working %s", fn.c_str());
		PDBChain &Chain = *Chains[i];
		Ps(fn, "chain%u.pdb", i+1);
		Chain.ToPDB(fn);
		}
	}
