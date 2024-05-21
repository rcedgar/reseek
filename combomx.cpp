#include "myutils.h"
#include "dss.h"
#include "combomx.h"

void cmd_combomx()
	{
	FILE *f = CreateStdioFile(g_Arg1);

	vector<FEATURE> Fs;
	Fs.push_back(FEATURE_SS3);
	Fs.push_back(FEATURE_NbrSS3);
	Fs.push_back(FEATURE_RevNbrDist4);
	const uint NF = SIZE(Fs);

	vector<float **> ScoreMxs;
	ScoreMxs.push_back(g_ScoreMxs2[FEATURE_SS3]);
	ScoreMxs.push_back(g_ScoreMxs2[FEATURE_NbrSS3]);
	ScoreMxs.push_back(g_ScoreMxs2[FEATURE_RevNbrDist4]);

	DSS D;
	D.SetComboFeatures(Fs);
	uint AS = D.GetAlphaSize(FEATURE_Combo);
	vector<vector<float> > ComboMx(AS);
	for (uint i = 0; i < AS; ++i)
		{
		ComboMx[i].resize(AS);

		vector<uint> Lettersi;
		D.GetComboLetters(i, Lettersi);
		asserta(SIZE(Lettersi) == NF);

		for (uint j = 0; j < AS; ++j)
			{
			vector<uint> Lettersj;
			D.GetComboLetters(j, Lettersj);
			asserta(SIZE(Lettersj) == NF);

			float Score = 0;
			for (uint k = 0; k < NF; ++k)
				{
				uint Letteri = Lettersi[k];
				uint Letterj = Lettersj[k];
				Score += ScoreMxs[k][Letteri][Letterj];
				}
			ComboMx[i][j] = Score;
			}
		}

	fprintf(f, "\n");
	fprintf(f, "float ScoreMx_%s[%u][%u] = {\n", "Combo", AS, AS);
	for (uint i = 0; i < AS; ++i)
		{
		fprintf(f, "  {");
		for (uint j = 0; j < AS; ++j)
			fprintf(f, " %5.2ff,", ComboMx[i][j]);
		fprintf(f, "  }, // %u\n", i);
		}
	fprintf(f, "};\n");

	fprintf(f, "\n");
	fprintf(f, "\n");
	fprintf(f, "int IntScoreMx_%s[%u][%u] = {\n", "Combo", AS, AS);
	for (uint i = 0; i < AS; ++i)
		{
		fprintf(f, "  {");
		for (uint j = 0; j < AS; ++j)
			fprintf(f, " %3d,", int(2*ComboMx[i][j]));
		fprintf(f, "  }, // %u\n", i);
		}
	fprintf(f, "};\n");

	CloseStdioFile(f);
	}
