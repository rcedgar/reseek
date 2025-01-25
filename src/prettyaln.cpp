#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "dssaligner.h"

double Kabsch(const PDBChain &ChainA, const PDBChain &ChainB,
  uint LoA, uint LoB, const string &Path,
  double t[3], double u[3][3]);
void WriteLocalAln(FILE *f, const string &LabelA, const byte *A,
  const string &LabelB, const byte *B,
  uint Loi, uint Loj, const char *Path);

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

void DSSAligner::PrettyAln(FILE *f,
  const PDBChain &A, const PDBChain &B,
  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB,
  uint LoA, uint LoB, const string &Path, float AQ, float Evalue) const
	{
	if (f == 0)
		return;
	const string &LabelA = A.m_Label;
	const string &LabelB = B.m_Label;
	const string &SeqA = A.m_Seq;
	const string &SeqB = B.m_Seq;
	const uint LA = SIZE(SeqA);
	const uint LB = SIZE(SeqB);

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
	uint ColLo = 0;
	fprintf(f, "\n");
	fprintf(f, "_____________________________________________________________________________________________________________\n");
	double PctId = GetPct(Ids, ColCount);
	double PctGaps = GetPct(Gaps, ColCount);

	const byte *ByteSeqA = (const byte *) A.m_Seq.c_str();
	const byte *ByteSeqB = (const byte *) B.m_Seq.c_str();
	WriteLocalAln(f, LabelA, ByteSeqA, LabelB, ByteSeqB, LoA, LoB, Path.c_str());

	fprintf(f, "%s %u-%u length %u\n",
	  A.m_Label.c_str(), LoA + 1, PosA, LA);

	fprintf(f, "%s %u-%u length %u\n",
	  B.m_Label.c_str(), LoB + 1, PosB, LB);

	fprintf(f, "AQ %.4f, cols %u, gaps %u (%.1f%%), ids %u (%.1f%%)",
	  AQ, ColCount, Gaps, PctGaps, Ids, PctId);
	if (Evalue != FLT_MAX)
		fprintf(f, ", E-value %.3g", Evalue);
	fprintf(f, "\n");
	}
