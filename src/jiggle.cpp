#include "myutils.h"
#include "sweeper.h"
#include "sort.h"

static float randf(float lo, float hi)
	{
	const uint M = 1000003;
	uint r = randu32()%M;
	float f = float(r)/(M-1);
	float v = lo + f*(hi - lo);
	return v;
	}

static void Jiggle(float &v)
	{
	if (v == 0)
		{
		if (randu32()%16 == 0)
			{
			uint r = randu32()%100;
			v = randf(0.1f, 1.5f);
			}
		}
	else
		{
		if (randu32()%16 == 0)
			{
			v = 0;
			return;
			}
		float factor = randf(0.7f, 1.3f);
		float adder = randf(-0.03f, 0.03f);
		v = v*factor + adder;
		if (v < 0)
			v = -v;
		}
	}

void cmd_jiggle()
	{
	Die("TODO");

	//asserta(optset_tsv);
	//Sweeper::m_ftsv = CreateStdioFile(opt(tsv));
	//setbuf(Sweeper::m_ftsv, 0);

	//const string &CalFN = g_Arg1;

	//Sweeper S;

	//DSSParams Params;
	//Params.SetDefaults();
	//S.m_SB.Setup(Params, CalFN);

	//const uint N = Params.GetFeatureCount();
	//vector<float> Values = Params.m_Weights;

	//S.m_StepName = "Init";
	//uint BestScore = S.Run(Params, Values);
	//vector<float> BestValues = Values;
	//for (;;)
	//	{
	//	Values = BestValues;
	//	for (uint Idx = 0; Idx < N; ++Idx)
	//		Jiggle(Values[Idx]);
	//	uint Score = S.Run(Params, Values);
	//	if (Score > BestScore)
	//		{
	//		BestScore = Score;
	//		BestValues = Values;
	//		}
	//	}
	}
