#include "myutils.h"
#include "pdbchain.h"

// Based on sec_str() in TMalign.cpp
//  Zhang & Skolnick 2005
static char GetSSChar(
  double dis13, double dis14, double dis15,
  double dis24, double dis25, double dis35)
	{
// Helix
	const double DH = 2.1;
	if (fabs(dis15 - 6.37) < DH && fabs(dis14 - 5.18) < DH &&
		fabs(dis25 - 5.18) < DH && fabs(dis13 - 5.45) < DH &&
		fabs(dis24 - 5.45) < DH && fabs(dis35 - 5.45) < DH)
		return 'h';

// Strand
	const double DS = 1.42;
	if (fabs(dis15 - 13) < DS && fabs(dis14 - 10.4) < DS &&
		fabs(dis25 - 10.4) < DS && fabs(dis13 - 6.1) < DS &&
		fabs(dis24 - 6.1) < DS && fabs(dis35 - 6.1) < DS)
		return 's';

// Turn
//	if (dis15 < 8.0)
	if (dis15 < 8.2)
		return 't';

// Default to loop, no well-defined ss
	return '~';
	}

void PDBChain::GetSS(string &SS) const
	{
	SS.clear();
	uint L = int(SIZE(m_Seq));
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		if (Pos < 2 || Pos + 2 >= L)
			{
			SS += '~';
			continue;
			}

		uint Pos_minus_2 = Pos - 2;
		uint Pos_minus_1 = Pos - 1;
		uint Pos_plus_1 = Pos + 1;
		uint Pos_plus_2 = Pos + 2;

		double d13 = GetDist(Pos_minus_2, Pos);
		double d14 = GetDist(Pos_minus_2, Pos_plus_1);
		double d15 = GetDist(Pos_minus_2, Pos_plus_2);
		double d24 = GetDist(Pos_minus_1, Pos_plus_1);
		double d25 = GetDist(Pos_minus_1, Pos_plus_2);
		double d35 = GetDist(Pos, Pos_plus_2);

		char c = GetSSChar(d13, d14, d15, d24, d25, d35);
		SS += c;
		}
	}

void cmd_pdb2ss()
	{
	const string &QueryFN = opt(pdb2ss);

	vector<PDBChain *> Qs;
	ReadChains(QueryFN, Qs);

	const uint N = SIZE(Qs);
	for (uint i = 0; i < N; ++i)
		{
		PDBChain &Q = *Qs[i];
		uint QL = Q.GetSeqLength();
		string Sketch;
		string SS;
		Q.GetSS(SS);
		Log("%s   SecStr  %s\n", Q.m_Label.c_str(), SS.c_str());
		}
	}
