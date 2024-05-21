#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"

double Kabsch(const PDBChain &ChainA, const PDBChain &ChainB,
  uint LoA, uint LoB, const string &Path,
  vector<double> &t, vector<vector<double> > &R);

static char GetDistSymbol(double d, bool Id)
	{
	if (d <= 5 && Id)
		return '|';
	if (d <= 1)
		{
		return '*';
		}
	if (d <= 2.5)
		return '+';
	if (d <= 5)
		return ':';
	if (d <= 10)
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

	vector<double> t;
	vector<vector<double> > R;
	double RMS = Kabsch(A, B, LoA, LoB, Path, t, R);
	PDBChain RB;
	if (RMS > 0)
		{
		t[0] = -t[0];
		t[1] = -t[1];
		t[2] = -t[2];
		B.GetXFormChain_tR(t, R, RB);
		}

	string RowA;
	string RowB;
	string RowDist;
	const uint ColCount = SIZE(Path);
	uint PosA = LoA;
	uint PosB = LoB;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			{
			double d = 0;
			if (RMS > 0)
				{
				vector<double> PtA;
				vector<double> PtRB;
				A.GetPt(PosA, PtA);
				RB.GetPt(PosB, PtRB);
				d = GetDist(PtA, PtRB);
				}
			asserta(PosA < LA);
			asserta(PosB < LB);
			char a = SeqA[PosA];
			char b = SeqB[PosB];
			RowA += a;
			RowB += b;
			RowDist += GetDistSymbol(d, a==b);
			++PosA;
			++PosB;
			break;
			}

		case 'D':
			asserta(PosA < LA);
			RowA += SeqA[PosA];
			RowB += '-';
			RowDist += ' ';
			++PosA;
			break;

		case 'I':
			asserta(PosB < LB);
			RowA += '-';
			RowB += SeqB[PosB];
			RowDist += ' ';
			++PosB;
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
		if (RMS > 0)
			{
			for (uint Col = ColLo; Col < ColHi; ++Col)
				fprintf(f, "%c", RowDist[Col]);
			fprintf(f, "\n");
			}
		for (uint Col = ColLo; Col < ColHi; ++Col)
			fprintf(f, "%c", RowB[Col]);
		fprintf(f, "  %s\n\n", B.m_Label.c_str());
		if (ColHi == ColCount)
			break;
		asserta(ColHi < ColCount);
		ColLo = ColHi + 1;
		}
	if (Evalue >= 0)
		fprintf(f, "Evalue %.3g\n", Evalue);
	}
