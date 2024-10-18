#include "myutils.h"

void WriteAnnotRow(FILE *f, const byte *A, const byte *B, const char *Path,
  unsigned i, unsigned j, unsigned ColLo, unsigned ColHi)
	{
	fprintf(f, "%5.5s ", "");
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M')
			{
			byte a = A[i++];
			byte b = B[j++];
			if (toupper(a) == toupper(b))
				fprintf(f, "|");
			else
				fprintf(f, " ");
			}
		else
			{
			if (c == 'D')
				++i;
			else if (c == 'I')
				++j;
			else
				asserta(false);
			fprintf(f, " ");
			}
		}
	fprintf(f, "\n");
	}

void WriteBRow(FILE *f, const byte *B, const char *Path,
  unsigned &j, unsigned ColLo, unsigned ColHi, const string &LabelB)
	{
	fprintf(f, "%5u ", j+1);
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M' || c == 'I')
			fprintf(f, "%c", B[j++]);
		else
			fprintf(f, "-");
		}
	fprintf(f, " %u  %s\n", j, LabelB.c_str());
	}

void WriteARow(FILE *f, const byte *A, const char *Path,
  unsigned &i, unsigned ColLo, unsigned ColHi, const string &LabelA)
	{
	fprintf(f, "%5u ", i+1);
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M' || c == 'D')
			fprintf(f, "%c", A[i++]);
		else
			fprintf(f, "-");
		}
	fprintf(f, " %u  %s\n", i, LabelA.c_str());
	}

void WriteLocalAln(FILE *f, const string &LabelA, const byte *A,
  const string &LabelB, const byte *B,
  uint Loi, uint Loj, const char *Path)
	{
	unsigned BLOCK_SIZE = 80;
	if (optset_rowlen)
		BLOCK_SIZE = opt_rowlen;
	uint ColLo = 0;
	uint ColHi = (unsigned) strlen(Path) - 1;

	asserta(ColHi >= ColLo);

	unsigned PosA = Loi;
	unsigned PosB = Loj;
	unsigned ColFrom = ColLo;
	for (;;)
		{
		if (ColFrom > ColHi)
			break;
		unsigned ColTo = ColFrom + BLOCK_SIZE - 1;
		if (ColTo > ColHi)
			ColTo = ColHi;

		unsigned i0 = PosA;
		unsigned j0 = PosB;
		WriteARow(f, A, Path, PosA, ColFrom, ColTo, LabelA);
		WriteAnnotRow(f, A, B, Path, i0, j0, ColFrom, ColTo);
		WriteBRow(f, B, Path, PosB, ColFrom, ColTo, LabelB);
		fprintf(f, "\n");

		ColFrom += BLOCK_SIZE;
		}
	}
