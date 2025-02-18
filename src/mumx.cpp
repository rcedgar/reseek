#include "myutils.h"
#include "dss.h"
#include "mumx.h"
#include "mermx.h"

static int8_t intround(float f)
	{
	int i = int(round(f));
	asserta(i >= INT8_MIN && i <= INT8_MAX);
	return int8_t(i);
	}

static MerMx *s_ptrMuMerMx = 0;

const MerMx &GetMuMerMx()
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
	(*s_ptrMuMerMx).Init(MxPtrs, 5, 36, 2);
	return *s_ptrMuMerMx;
	}

void cmd_combomx()
	{
	FILE *f = CreateStdioFile(g_Arg1);

	vector<FEATURE> Fs;
	Fs.push_back(FEATURE_SS3);
	Fs.push_back(FEATURE_NENSS3);
	Fs.push_back(FEATURE_RENDist4);
	const uint NF = SIZE(Fs);

	vector<float **> ScoreMxs;
	ScoreMxs.push_back(g_ScoreMxs2[FEATURE_SS3]);
	ScoreMxs.push_back(g_ScoreMxs2[FEATURE_NENSS3]);
	ScoreMxs.push_back(g_ScoreMxs2[FEATURE_RENDist4]);

	DSS D;
	DSSParams Params;
	Params.SetMuFeatures(Fs);
	D.SetParams(Params);
	uint AS = D.GetAlphaSize(FEATURE_Mu);
	vector<vector<float> > MuMx(AS);
	for (uint i = 0; i < AS; ++i)
		{
		MuMx[i].resize(AS);

		vector<uint> Lettersi;
		D.GetMuLetters(i, Lettersi);
		asserta(SIZE(Lettersi) == NF);

		for (uint j = 0; j < AS; ++j)
			{
			vector<uint> Lettersj;
			D.GetMuLetters(j, Lettersj);
			asserta(SIZE(Lettersj) == NF);

			float Score = 0;
			for (uint k = 0; k < NF; ++k)
				{
				uint Letteri = Lettersi[k];
				uint Letterj = Lettersj[k];
				Score += ScoreMxs[k][Letteri][Letterj];
				}
			MuMx[i][j] = Score;
			}
		}

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
	fprintf(f, "int IntScoreMx_%s[%u][%u] = {\n", "Mu", AS, AS);
	for (uint i = 0; i < AS; ++i)
		{
		fprintf(f, "  {");
		for (uint j = 0; j < AS; ++j)
			fprintf(f, " %3d,", intround(MuMx[i][j]));
		fprintf(f, "  }, // %u\n", i);
		}
	fprintf(f, "};\n");

	fprintf(f, "\n");
	fprintf(f, "\n");
	fprintf(f, "int IntScoreMx_%s[%u][%u] = {\n", "Mu_x2", AS, AS);
	for (uint i = 0; i < AS; ++i)
		{
		fprintf(f, "  {");
		for (uint j = 0; j < AS; ++j)
			fprintf(f, " %3d,", intround(2*MuMx[i][j]));
		fprintf(f, "  }, // %u\n", i);
		}
	fprintf(f, "};\n");

	const char *MuAlphaStr = "ABCDEFGHIJKLMNOPQRSTUVWZYZabcdefghij";
	size_t n = strlen(MuAlphaStr);
	asserta(n == 36);

	fprintf(f, "\n");
	fprintf(f, "\n");
	fprintf(f, "static const int parasail_combo_[%u*%u] = {\n", AS, AS);
	int MinScore = INT_MAX;
	int MaxScore = INT_MIN;
	for (uint i = 0; i < AS; ++i)
		{
		for (uint j = 0; j < AS; ++j)
			{
			int Score = intround(MuMx[i][j]);
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
	fprintf(f, "static const int parasail_combo_map[256] = {");
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
	fprintf(f, "static const parasail_matrix_t parasail_combo = {\n");
	fprintf(f, "	\"combo\",\n");
	fprintf(f, "	parasail_combo_,\n");
	fprintf(f, "	parasail_combo_map,\n");
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
