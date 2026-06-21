#if 0
#include "myutils.h"
#include "dssaligner.h"
#include "chainreader2.h"

static void DumpInt(FILE *f, uint i)
	{
	if (i == UINT_MAX)
		fprintf(f, "\t.");
	else
		fprintf(f, "\t%u", i);
	}

static void Dump(FILE *f, DSS &D)
	{
	const uint L = D.m_Chain->GetSeqLength();

	fprintf(f, "pos");
	fprintf(f, "\taa");
	fprintf(f, "\tNEN");
	fprintf(f, "\tREN");
	fprintf(f, "\tPEN");
	fprintf(f, "\tMEN");
	for (uint fi = 0; fi < FEATURE_COUNT; ++fi)
		{
		if (fi == FEATURE_SSSA || fi == FEATURE_SSSB) continue;
		fprintf(f, "\t%s", FeatureToStr(fi));
		}
	fprintf(f, "\tf_NENDist");
	fprintf(f, "\tf_RENDist");
	fprintf(f, "\tf_PENDist");
	fprintf(f, "\tf_MENDist");
	fprintf(f, "\n");
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		fprintf(f, "%u", Pos);
		fprintf(f, "\t%c", D.m_Chain->m_Seq[Pos]);
		DumpInt(f, D.GetNEN(Pos));
		DumpInt(f, D.GetREN(Pos));
		DumpInt(f, D.GetPEN(Pos));
		DumpInt(f, D.GetMEN(Pos));
		for (uint fi = 0; fi < FEATURE_COUNT; ++fi)
			{
			if (fi == FEATURE_SSSA || fi == FEATURE_SSSB) continue;
			uint letter = D.GetFeature(fi, Pos);
			fprintf(f, "\t%u", letter);
			}
		fprintf(f, "\t%.3g", D.GetFloatFeature(FEATURE_NENDist, Pos));
		fprintf(f, "\t%.3g", D.GetFloatFeature(FEATURE_RENDist, Pos));
		fprintf(f, "\t%.3g", D.GetFloatFeature(FEATURE_PENDist, Pos));
		fprintf(f, "\t%.3g", D.GetFloatFeature(FEATURE_MENDist, Pos));
		fprintf(f, "\n");
		}
	}

static void LogNEN(DSS &D, uint Pos)
	{
	// int DSSParams::m_NEN_w = 12;
	const int w = 12;

	const uint L = D.GetSeqLength();
	uint NEN = D.GetNEN(Pos);
	uint NEN2 = UINT_MAX;
	float mind = FLT_MAX;
	for (uint i = 0; i < L; ++i)
		{
		if (abs(int(Pos) - int(i)) < w)
			{
			Log("[%3u]  (too close)\n", i);
			continue;
			}
		float d = D.m_Chain->GetDist(Pos, i);
		Log("[%3u]  %7.1f\n", i, d);
		if (d < mind)
			{
			NEN = i;
			mind = d;
			}
		}
	Log("\nNEN=%u, d=%.1f\n", NEN, mind);
	}

void cmd_dump_dss()
	{
	const string &QFN = g_Arg1;
	FILE *fOut = CreateStdioFile(opt(output));

	DSSParams::Init(DM_AlwaysSensitive);
	if (optset_varstr)
		DSSParams::SetParamsFromStr(opt(varstr));

	DSSParams::m_Omega8 = 0;
	DSSParams::m_Omega16 = 0;

	DSSAligner DA;
	DSS D;

	vector<vector<byte> > Profile;
	ChainReader2 CR;
	CR.Open(QFN);
	uint N = 0;
	for (;;)
		{
		PDBChain *Chain = CR.GetNext();
		if (Chain == 0)
			break;
		++N;
		if (N%100 == 0)
			Progress("%u\r", N);

		DSS D;
		D.Init(*Chain);
		LogNEN(D, 16);
		break;
		Dump(fOut, D);
		}
	CloseStdioFile(fOut);
	}

#endif
