#include "myutils.h"

void WriteAnnotRow(
	FILE *f,
	const char *A,
	const char *B,
	const char *Path,
	unsigned i, unsigned j,
	unsigned ColLo, unsigned ColHi)
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

void WriteBRow(
	FILE *f,
	const char *B,
	const char *Path,
	unsigned &j,
	unsigned ColLo,
	unsigned ColHi,
	const char *LabelB)
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
	fprintf(f, " %u  %s\n", j, LabelB);
	}

void WriteARow(
	FILE *f,
	const char *A,
	const char *Path,
	unsigned &i,
	unsigned ColLo,
	unsigned ColHi,
	const char *LabelA)
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
	fprintf(f, " %u  %s\n", i, LabelA);
	}

void WriteLocalAln(
	FILE *f,
	const char *LabelA, const char *A, uint LA,
	const char *LabelB, const char *B, uint LB,
	uint loA, uint loB,
	const char *Path, uint ncol)
	{
	if (f == 0) return;
	if (ncol == 0) return;
	const unsigned BLOCK_SIZE = 80;
	uint ColLo = 0;
	uint ColHi = ncol - 1;

	asserta(ColHi >= ColLo);

	unsigned PosA = loA;
	unsigned PosB = loB;
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

// human-readable blast-like alignment
// with aa sequence rows
void human_aln(FILE *f,
	const char *labelA, const char *seqA, uint LA,
	const char *labelB, const char *seqB, uint LB,
	uint LoA, uint LoB, const char *path, uint ncol, 
	float pvalue)
	{
	if (f == 0)
		return;

	uint PosA = LoA;
	uint PosB = LoB;
	uint Ids = 0;
	uint Gaps = 0;
	for (uint Col = 0; Col < ncol; ++Col)
		{
		char c = path[Col];
		switch (c)
			{
		case 'M':
			{
			asserta(PosA < LA);
			asserta(PosB < LB);
			char a = seqA[PosA];
			char b = seqB[PosB];
			++PosA;
			++PosB;
			if (a == b) ++Ids;
			break;
			}

		case 'D':
			asserta(PosA < LA);
			++PosA;
			++Gaps;
			break;

		case 'I':
			asserta(PosB < LB);
			++PosB;
			++Gaps;
			break;

		default:
			asserta(false);
			}
		}
	double PctId = GetPct(Ids, ncol);
	double PctGaps = GetPct(Gaps, ncol);

	uint seglenA = PosA - LoA;
	uint seglenB = PosB - LoB;
	double pctcovA = GetPct(seglenA, LA);
	double pctcovB = GetPct(seglenB, LB);

	static mutex s_lock;
	s_lock.lock();

	fprintf(f, "\n");
	fprintf(f, "_____________________________________________________________________________________________________________\n");

	WriteLocalAln(f,
		labelA, seqA, LA,
		labelB, seqB, LB,
		LoA, LoB, path, ncol);

	fprintf(f, "Qry %u-%u/%u (%.1f%%) >%s\n",
		LoA + 1, PosA, LA, pctcovA, labelA);
	fprintf(f, " DB %u-%u/%u (%.1f%%) >%s\n",
		LoB + 1, PosB, LB, pctcovB, labelB);

	if (pvalue != FLT_MAX)
		fprintf(f, "P-value %.3g, ", pvalue);
	fprintf(f, "cols %u, gaps %u (%.1f%%), ids %u (%.1f%%)\n",
	  ncol, Gaps, PctGaps, Ids, PctId);
	fprintf(f, "\n");
	s_lock.unlock();
	}
