#include "myutils.h"
#include "dss.h"
#include "mumx.h"
#include "mermx.h"

#pragma warning(disable:4305) // double -> float

/********************************
FEATURE_SS3
FEATURE_NENSS3
FEATURE_RENDist4

git checkout 8e89307
scoremxs2.cpp

Expected scores
0.3531   4  RevNbrDist4
0.4655   3  SS3
0.2556   3  NbrSS3

Freqs SS3/3
   [ 0]  0.338479
   [ 1]  0.209130
   [ 2]  0.452391

Freqs NbrSS3/3
   [ 0]  0.281821
   [ 1]  0.283820
   [ 2]  0.434359

Freqs RevNbrDist4/4
   [ 0]  0.191469
   [ 1]  0.452078
   [ 2]  0.151362
   [ 3]  0.205091
********************************/

float ScoreMx_SS3[3][3] = {
  {  0.8913, -2.8877, -1.1241,  }, // 0
  { -2.8877,  1.4429, -0.5733,  }, // 1
  { -1.1241, -0.5733,  0.4515,  }, // 2
};

float ScoreMx_NbrSS3[3][3] = {
  {  0.7132, -1.4900, -0.6652,  }, // 0
  { -1.4900,  1.0304, -0.3502,  }, // 1
  { -0.6652, -0.3502,  0.3824,  }, // 2
};

float ScoreMx_RevNbrDist4[4][4] = {
  {  1.4279, -0.4181, -1.4284, -2.4223,  }, // 0
  { -0.4181,  0.4804, -0.4363, -0.9914,  }, // 1
  { -1.4284, -0.4363,  1.0964, -0.3005,  }, // 2
  { -2.4223, -0.9914, -0.3005,  0.9618,  }, // 3
};

static int8_t intround(float f)
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
	FILE *f = CreateStdioFile(g_Arg1);

	//vector<FEATURE> Fs;
	//Fs.push_back(FEATURE_SS3);
	//Fs.push_back(FEATURE_NENSS3);
	//Fs.push_back(FEATURE_RENDist4);
	const uint NF = 3;
	const uint AS = 36;

	//vector<float **> ScoreMxs;
	//ScoreMxs.push_back(g_ScoreMxs2[FEATURE_SS3]);
	//ScoreMxs.push_back(g_ScoreMxs2[FEATURE_NENSS3]);
	//ScoreMxs.push_back(g_ScoreMxs2[FEATURE_RENDist4]);

	vector<vector<vector<float> > > ScoreMxs(3);
	vector<vector<float> > &SS3 = ScoreMxs[0];
	vector<vector<float> > &NENSS3 = ScoreMxs[1];
	vector<vector<float> > &RENDist4 = ScoreMxs[2];

	SS3.resize(3);
	NENSS3.resize(3);
	RENDist4.resize(4);

	for (uint i = 0; i < 3; ++i)
		{
		SS3[i].resize(3);
		NENSS3[i].resize(3);
		for (uint j = 0; j < 3; ++j)
			{
			SS3[i][j] = (float) ScoreMx_SS3[i][j];
			NENSS3[i][j] = (float) ScoreMx_NbrSS3[i][j];
			}
		}

	for (uint i = 0; i < 4; ++i)
		{
		RENDist4[i].resize(4);
		for (uint j = 0; j < 3; ++j)
			RENDist4[i][j] = (float) ScoreMx_RevNbrDist4[i][j];
		}

	DSS D;
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
	fprintf(f, "static const int parasail_mu_[%u*%u] = {\n", AS, AS);
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
