#include "myutils.h"
#include "dss.h"
#include "mumx.h"
#include "mermx.h"
#include "featuretrainer.h"

static int8_t int8round(float f)
	{
	int i = int(round(f));
	asserta(i >= INT8_MIN && i <= INT8_MAX);
	return int8_t(i);
	}

static MerMx *s_ptrMuMerMx = 0;

const MerMx &GetMuMerMx(uint k)
	{
	if (s_ptrMuMerMx != 0)
		return *s_ptrMuMerMx;
	s_ptrMuMerMx = new MerMx;

	short **MxPtrs = myalloc(short *, 36);
	for (uint i = 0; i < 36; ++i)
		{
		short *Row = myalloc(short, 36);
		for (uint j = 0; j < 36; ++j)
			Row[j] = Mu_S_ij_i8[i][j];
		MxPtrs[i] = Row;
		}
	(*s_ptrMuMerMx).Init(MxPtrs, k, 36, 2);
	return *s_ptrMuMerMx;
	}

void cmd_musubstmx()
	{
	FeatureTrainer FT;
	FT.FromTsv(g_Arg1);
	const uint AS = FT.m_AlphaSize;

	vector<vector<float> > MuMx;
	FT.GetLogOddsMx(MuMx);

	FILE *f = CreateStdioFile(opt(output));

	fprintf(f, "\n");
	fprintf(f, "float ScoreMx_%s[%u][%u] = {\n", "Mu", AS, AS);
	for (uint i = 0; i < AS; ++i)
		{
		fprintf(f, "  {");
		for (uint j = 0; j < AS; ++j)
			fprintf(f, " %5.2ff,", MuMx[i][j]);
		fprintf(f, "  }, // %u\n", i);
		}
	fprintf(f, "};\n");

	fprintf(f, "\n");
	fprintf(f, "\n");
	fprintf(f, "int8_t IntScoreMx_%s[%u][%u] = {\n", "Mu", AS, AS);
	for (uint i = 0; i < AS; ++i)
		{
		fprintf(f, "  {");
		for (uint j = 0; j < AS; ++j)
			fprintf(f, "%3d,", int8round(MuMx[i][j]));
		fprintf(f, "  }, // %u\n", i);
		}
	fprintf(f, "};\n");

	fprintf(f, "\n");
	fprintf(f, "\n");
	fprintf(f, "int8_t int8_t Mu_S_ij_i8[%u][%u] = {\n", AS, AS);
	for (uint i = 0; i < AS; ++i)
		{
		fprintf(f, "  {");
		for (uint j = 0; j < AS; ++j)
			fprintf(f, "%3d,", int8round(3*MuMx[i][j]));
		fprintf(f, "  }, // %u\n", i);
		}
	fprintf(f, "};\n");

	const char *MuAlphaStr = "ABCDEFGHIJKLMNOPQRSTUVWZYZabcdefghij";
	size_t n = strlen(MuAlphaStr);
	asserta(n == 36);

	fprintf(f, "\n");
	fprintf(f, "\n");
	fprintf(f, "static const int parasail_mu_[%u*%u] = {\n", AS, AS);
	int MinScore = INT_MAX;
	int MaxScore = INT_MIN;
	for (uint i = 0; i < AS; ++i)
		{
		for (uint j = 0; j < AS; ++j)
			{
			int Score = int8round(MuMx[i][j]);
			MinScore = min(MinScore, Score);
			MaxScore = max(MaxScore, Score);
			fprintf(f, "%2d,", Score);
			}
		fprintf(f, "  // %u\n", i);
		}
	fprintf(f, "};\n");

	int imap[256];
	for (uint i = 0; i < 256; ++i)
		imap[i] = 0;
	for (uint i = 0; i < 36; ++i)
		{
		char c = MuAlphaStr[i];
		imap[c] = i;
		}

	fprintf(f, "\n");
	fprintf(f, "static const int parasail_mu_map[256] = {");
	for (uint i = 0; i < 256; ++i)
		{
		fprintf(f, " %2d,", imap[i]);
		if (isalnum(i))
			fprintf(f, " // %c\n", byte(i));
		else
			fprintf(f, " // %02X\n", byte(i));
		}
	fprintf(f, "};\n");

	fprintf(f, "\n");
	fprintf(f, "static const parasail_matrix_t parasail_mu = {\n");
	fprintf(f, "	\"mu\",\n");
	fprintf(f, "	parasail_mu_,\n");
	fprintf(f, "	parasail_mu_map,\n");
	fprintf(f, "	%u,\n", AS);
	fprintf(f, "	%d,\n", MaxScore);
	fprintf(f, "	%d,\n", MinScore);
	fprintf(f, "	NULL,\n");
	fprintf(f, "	PARASAIL_MATRIX_TYPE_SQUARE,\n");
	fprintf(f, "	%u,\n", AS);
	fprintf(f, "	\"%s\",\n", MuAlphaStr);
	fprintf(f, "	NULL\n");
	fprintf(f, "};\n");

	CloseStdioFile(f);
	}
