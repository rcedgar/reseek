#include "myutils.h"
#include "cigar.h"

const int SIZE_32 = 32;

/***
https://samtools.github.io/hts-specs/SAMv1.pdf
 
6. CIGAR: CIGAR string. The CIGAR operations are given in the following table (set '*' if unavailable):

Op BAM	Description Consumes query Consumes reference
M	0	alignment match (can be a sequence match or mismatch) yes yes
I	1	insertion to the reference yes no
D	2	deletion from the reference no yes
N	3	skipped region from the reference no yes
S	4	soft clipping (clipped sequences present in SEQ) yes no
H	5	hard clipping (clipped sequences NOT present in SEQ) no no
P	6	padding (silent deletion from padded reference) no no
=	7	sequence match yes yes
X	8	sequence mismatch yes yes

(*) "Consumes query" and "consumes reference" indicate whether the CIGAR operation causes the
	  alignment to step along the query sequence and the reference sequence respectively.
(*) H can only be present as the first and/or last operation.
(*) S may only have H operations between them and the ends of the CIGAR string.
(*) For mRNA-to-genome alignment, an N operation represents an intron. For other types of alignments,
	  the interpretation of N is not defined.
(*) Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
***/

// Query ("read") gets Snnn for soft-clip at start
// In SAM, reference position is tabbed field 4, here
//   T is used for number of clipped letters at start
//   in the reference.
void LocalPathToCIGAR(const char *Path, uint LoQ, uint LoR,
  string &CIGAR, bool FlipDI)
	{
	CIGAR.clear();
	char LastC = *Path;
	uint n = 1;
	char Tmp[SIZE_32];

	if (LoQ > 0)
		{
		snprintf(Tmp, SIZE_32, "%uS", LoQ);
		CIGAR += string(Tmp);
		}
	if (LoR > 0)
		{
		snprintf(Tmp, SIZE_32, "%uT", LoR);
		CIGAR += string(Tmp);
		}

	for (uint i = 1; ; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (FlipDI)
			{
			if (c == 'D')
				c = 'I';
			else if (c == 'I')
				c = 'D';
			}
		if (c == LastC)
			{
			++n;
			continue;
			}
		else
			{
			assert(n > 0);
			if (LastC == 'D')
				LastC = 'I';
			else if (LastC == 'I')
				LastC = 'D';
			snprintf(Tmp, SIZE_32, "%u%c", n, LastC);
			CIGAR += string(Tmp);
			LastC = c;
			n = 1;
			}
		}
	if (n > 0)
		{
		if (LastC == 'D')
			LastC = 'I';
		else if (LastC == 'I')
			LastC = 'D';
		snprintf(Tmp,  SIZE_32, "%u%c", n, LastC);
		CIGAR += string(Tmp);
		}
	}

void PathToCIGAR(const char *Path, string &CIGAR, bool FlipDI)
	{
	CIGAR.clear();
	char LastC = *Path;
	uint n = 1;
	char Tmp[SIZE_32];
	for (uint i = 1; ; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == LastC)
			{
			++n;
			continue;
			}
		else
			{
			assert(n > 0);
			if (!FlipDI)
				{
				if (LastC == 'D')
					LastC = 'I';
				else if (LastC == 'I')
					LastC = 'D';
				}
			snprintf(Tmp, SIZE_32, "%u%c", n, LastC);
			CIGAR += string(Tmp);
			LastC = c;
			n = 1;
			}
		}
	if (n > 0)
		{
		if (!FlipDI)
			{
			if (LastC == 'D')
				LastC = 'I';
			else if (LastC == 'I')
				LastC = 'D';
			}
		snprintf(Tmp, SIZE_32, "%u%c", n, LastC);
		CIGAR += string(Tmp);
		}
	}

void CIGARGetOps(const string &CIGAR, string &Ops, vector<uint> &Lengths)
	{
	Ops.clear();
	Lengths.clear();
	if (CIGAR.empty())
		return;

	uint L = SIZE(CIGAR);
	uint n = 0;
	for (uint i = 0; i < L; ++i)
		{
		char c = CIGAR[i];
		if (isdigit(c))
			n = n*10 + (c - '0');
		else if (isupper(c) || c == '=')
			{
			if (n == 0)
				Die("Operation '%c' has zero length in CIGAR '%s'", c, CIGAR.c_str());
			Ops.push_back(c);
			Lengths.push_back(n);
			n = 0;
			}
		else
			Die("Invalid char '%c' in CIGAR '%s'", c, CIGAR.c_str());
		}
	if (n > 0)
		Die("Missing operation at end of CIGAR '%s'", CIGAR.c_str());
	}

const char *LocalCIGARToPath(const string &CIGAR, string &Path,
  uint &LoQ, uint &LoR, bool FlipDI)
	{
	Path.clear();
	LoQ = 0;
	LoR = 0;

	string Ops;
	vector<uint> OpLengths;
	CIGARGetOps(CIGAR, Ops, OpLengths);

	const uint n = SIZE(Ops);
	asserta(SIZE(OpLengths) == n);
	for (uint i = 0; i < n; ++i)
		{
		char Op = Ops[i];
		if (FlipDI)
			{
			if (Op == 'D')
				Op = 'I';
			else if (Op == 'I')
				Op = 'D';
			}
		uint OpLength = OpLengths[i];
		if (Op == 'S')
			LoQ = OpLength;
		else if (Op == 'T')
			LoR = OpLength;
		else
			{
			for (uint j = 0; j < OpLength; ++j)
				Path += Op;
			}
		}
	return Path.c_str();
	}

const char *CIGARToPath(const string &CIGAR, string &Path)
	{
	Path.clear();

	string Ops;
	vector<uint> OpLengths;
	CIGARGetOps(CIGAR, Ops, OpLengths);

	const uint n = SIZE(Ops);
	asserta(SIZE(OpLengths) == n);
	for (uint i = 0; i < n; ++i)
		{
		char Op = Ops[i];
		uint OpLength = OpLengths[i];
		if (Op == 'M' || Op == 'D' || Op == 'I')
			{
			for (uint j = 0; j < OpLength; ++j)
				Path += Op;
			}
		else if (Op == 'S')
			{
			for (uint j = 0; j < OpLength; ++j)
				Path += 'I';
			}
		else if (Op == 'T')
			{
			for (uint j = 0; j < OpLength; ++j)
				Path += 'D';
			}
		else
			Die("Bad op in CIGAR '%c'", Op);
		}
	return Path.c_str();
	}

uint CIGARToQL(const string &CIGAR)
	{
	string Ops;
	vector<uint> Lengths;
	CIGARGetOps(CIGAR, Ops, Lengths);
	const uint N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	uint QL = 0;
	for (uint i = 0; i < N; ++i)
		{
		char Op = Ops[i];
		switch (Op)
			{
		case 'M':
		case 'I':
			QL += Lengths[i];
			break;

		case 'D':
			break;

		default:
			Die("Unsupported op '%c' in CIGAR '%s'", Op, CIGAR.c_str());
			}
		}
	return QL;
	}

void CIGAROpsToLs(const string &Ops, const vector<uint> &Lengths,
  uint &QL, uint &TL)
	{
	QL = 0;
	TL = 0;
	const uint N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	for (uint i = 0; i < N; ++i)
		{
		char Op = Ops[i];
		switch (Op)
			{
		case 'M':
			QL += Lengths[i];
			TL += Lengths[i];
			break;

	// CIGAR D&I reverse of my usual convention
		case 'I':
			QL += Lengths[i];
			break;

		case 'D':
			TL += Lengths[i];
			break;

		default:
			Die("Unsupported op '%c' in CIGAR", Op);
			}
		}
	}

void PathToLs(const string &Path, uint &QL, uint &TL)
	{
	QL = 0;
	TL = 0;
	for (uint i = 0; i < SIZE(Path); ++i)
		{
		switch (Path[i])
			{
		case 'M': ++QL; ++TL; break;
		case 'D': ++QL; break;
		case 'I': ++TL; break;
		default: asserta(false);
			}
		}
	}

void CIGARToLs(const string &CIGAR, uint &QL, uint &TL)
	{
	string Ops;
	vector<uint> Lengths;
	CIGARGetOps(CIGAR, Ops, Lengths);
	CIGAROpsToLs(Ops, Lengths, QL, TL);
	}

void CIGAROpsFixDanglingMs(string &Ops, vector<uint> &Lengths)
	{
	const uint N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	if (N < 3)
		return;

// 1M 6I 100M
	if (Ops[0] == 'M' && Lengths[0] <= 2 && Lengths[1] > 4 && Ops[2] == 'M')
		{
		uint OldQL;
		uint OldTL;
		CIGAROpsToLs(Ops, Lengths, OldQL, OldTL);

		string NewOps;
		vector<uint> NewLengths;
		for (uint i = 1; i < N; ++i)
			{
			NewOps.push_back(Ops[i]);
			NewLengths.push_back(Lengths[i]);
			}
		NewLengths[1] += Lengths[0];

		uint NewQL;
		uint NewTL;
		CIGAROpsToLs(NewOps, NewLengths, NewQL, NewTL);
		asserta(NewQL == OldQL);
		asserta(NewTL == OldTL);

		Ops = NewOps;
		Lengths = NewLengths;
		}

// 100M 6D M1
	if (Ops[N-1] == 'M' && Lengths[N-1] <= 2 && Lengths[N-2] > 4 && Ops[N-3] == 'M')
		{
		uint OldQL;
		uint OldTL;
		CIGAROpsToLs(Ops, Lengths, OldQL, OldTL);

		string NewOps;
		vector<uint> NewLengths;
		for (uint i = 0; i < N-1; ++i)
			{
			NewOps.push_back(Ops[i]);
			NewLengths.push_back(Lengths[i]);
			}
		NewLengths[N-3] += Lengths[N-1];

		uint NewQL;
		uint NewTL;
		CIGAROpsToLs(NewOps, NewLengths, NewQL, NewTL);
		asserta(NewQL == OldQL);
		asserta(NewTL == OldTL);

		Ops = NewOps;
		Lengths = NewLengths;
		}
	}

void OpsToCIGAR(const string &Ops, const vector<uint> &Lengths,
  string &CIGAR)
	{
	CIGAR.clear();
	const uint N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	for (uint i = 0; i < N; ++i)
		Psa(CIGAR, "%u%c", Lengths[i], Ops[i]);
	}

void ExpandParaCigar(const string &s, string &Path)
	{
	string Ops;
	vector<uint> ns;
	CIGARGetOps(s, Ops, ns);
	const uint N = SIZE(Ops);
	asserta(SIZE(ns) == N);
	for (uint i = 0; i < N; ++i)
		{
		char Op = Ops[i];
		if (Op == 'X' || Op == '=')
			Op = 'M';
		asserta(Op == 'M' || Op == 'D' || Op == 'I');
		for (uint j = 0; j < ns[i]; ++j)
			Path += Op;
		}
	}

void ExpandParaCigar_reverseDI(const string &s, string &Path)
	{
	string Ops;
	vector<uint> ns;
	CIGARGetOps(s, Ops, ns);
	const uint N = SIZE(Ops);
	asserta(SIZE(ns) == N);
	for (uint i = 0; i < N; ++i)
		{
		char Op = Ops[i];
		if (Op == 'X' || Op == '=')
			Op = 'M';
		if (Op == 'I')
			Op = 'D';
		else if (Op == 'D')
			Op = 'I';
		asserta(Op == 'M' || Op == 'D' || Op == 'I');
		for (uint j = 0; j < ns[i]; ++j)
			Path += Op;
		}
	}
