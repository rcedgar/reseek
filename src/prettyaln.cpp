#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"

double Kabsch(const PDBChain &ChainA, const PDBChain &ChainB,
  uint LoA, uint LoB, const string &Path,
  vector<double> &t, vector<vector<double> > &R);

static char GetAnnotChar(char a, char b)
	{
	float GetBlosum62Score(char a, char b);
	if (a == b)
		return '|';
	float Score = GetBlosum62Score(a, b);
	if (Score >= 2.0f)
		return ':';
	if (Score > 0)
		return '.';
	return ' ';
	}

void PrettyAln(FILE *f, const PDBChain &A, const PDBChain &B,
  uint LoA, uint LoB, const string &Path, float Evalue)
	{
	if (f == 0)
		return;
	const string &SeqA = A.m_Seq;
	const string &SeqB = B.m_Seq;
	const uint LA = SIZE(SeqA);
	const uint LB = SIZE(SeqB);

	string RowA;
	string RowB;
	string AnnotRow;
	const uint ColCount = SIZE(Path);
	uint PosA = LoA;
	uint PosB = LoB;
	uint Ids = 0;
	uint Gaps = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			{
			asserta(PosA < LA);
			asserta(PosB < LB);
			char a = SeqA[PosA];
			char b = SeqB[PosB];
			RowA += a;
			RowB += b;
			AnnotRow += GetAnnotChar(a, b);
			++PosA;
			++PosB;
			if (a == b) ++Ids;
			break;
			}

		case 'D':
			asserta(PosA < LA);
			RowA += SeqA[PosA];
			RowB += '-';
			AnnotRow += ' ';
			++PosA;
			++Gaps;
			break;

		case 'I':
			asserta(PosB < LB);
			RowA += '-';
			RowB += SeqB[PosB];
			AnnotRow += ' ';
			++PosB;
			++Gaps;
			break;

		default:
			asserta(false);
			}
		}
	const uint ROWLEN = 100;
	uint ColLo = 0;
	fprintf(f, "\n");
	fprintf(f, "_____________________________________________________________________________________________________________\n");
	for (;;)
		{
		uint ColHi = ColLo + ROWLEN;
		if (ColHi >= ColCount)
			ColHi = ColCount;
		for (uint Col = ColLo; Col < ColHi; ++Col)
			fprintf(f, "%c", RowA[Col]);
		fprintf(f, "  %s\n", A.m_Label.c_str());
		for (uint Col = ColLo; Col < ColHi; ++Col)
			fprintf(f, "%c", AnnotRow[Col]);
		fprintf(f, "\n");
		for (uint Col = ColLo; Col < ColHi; ++Col)
			fprintf(f, "%c", RowB[Col]);
		fprintf(f, "  %s\n\n", B.m_Label.c_str());
		if (ColHi == ColCount)
			break;
		asserta(ColHi < ColCount);
		ColLo = ColHi;
		}
	double PctId = GetPct(Ids, ColCount);
	double PctGaps = GetPct(Gaps, ColCount);

	fprintf(f, "%s %u-%u length %u\n",
	  A.m_Label.c_str(), LoA + 1, PosA, LA);

	fprintf(f, "%s %u-%u length %u\n",
	  B.m_Label.c_str(), LoB + 1, PosB, LB);

	fprintf(f, "Cols %u, gaps %u (%.1f%%), ids %u (%.1f%%), E-value %.3g\n",
	  ColCount, Gaps, PctGaps, Ids, PctId, Evalue);
	}
