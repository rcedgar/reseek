#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"
#include "scop40bench.h"

void TrainFeature(const SeqDB &Input, const vector<PDBChain *> &Chains,
				  FEATURE F, uint AlphaSize, 
				  vector<float> &BinTs, LogOdds &LO);
void WriteLO(FILE *f, const string &strName, const LogOdds &LO);
void WriteLOInt8(FILE *f, const string &strName, const LogOdds &LO, int8_t MaxAbsi8);

void cmd_train_features()
	{
	optset_fast = true;
	opt(fast) = true;
	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast, SCOP40_DBSIZE);
	int8_t MaxAbsi8 = 20;
	if (optset_maxi8)
		{
		uint Max = opt(maxi8);
		MaxAbsi8 = uint8_t(Max);
		asserta(uint(MaxAbsi8) == Max);
		}

	vector<string> FeatureNames;
	vector<FEATURE> Fs;
	if (!optset_features)
		{
		for (uint i = 0; i < Params.GetFeatureCount(); ++i)
			{
			FEATURE F = Params.m_Features[i];
			string Name = FeatureToStr(F);
			FeatureNames.push_back(Name);
			Fs.push_back(F);
			}
		}
	else
		Split(opt(features), FeatureNames, '_');

	vector<PDBChain *> Chains;
	ReadChains(opt(train_cal), Chains);

	SeqDB Input;
	Input.FromFasta(g_Arg1, true);

	const uint N = SIZE(FeatureNames);
	FILE *fOut = CreateStdioFile(opt(output));
	FILE *fOut2 = CreateStdioFile(opt(output2));
	LogOdds LO;
	for (uint i = 0; i < N; ++i)
		{
		FEATURE F = Fs[i];
		uint AlphaSize = GetAlphaSize(F);
		vector<float> BinTs;
		TrainFeature(Input, Chains, F, AlphaSize, BinTs, LO);
		WriteLO(fOut, FeatureNames[i], LO);
		WriteLOInt8(fOut2, FeatureNames[i], LO, MaxAbsi8);
		}
	CloseStdioFile(fOut);
	CloseStdioFile(fOut2);
	}
