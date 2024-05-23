#include "myutils.h"
#include "chainreader.h"
#include "dss.h"
#include "alpha.h"

extern float *g_FeatureFreqsVec[FEATURE_COUNT];

void cmd_pdb2mega()
	{
	if (!optset_output)
		Die("-output not set");
	optset_verysensitive = true;
	opt_verysensitive = true;

	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	const uint ChainCount = SIZE(Chains);

	DSSParams Params;
	Params.SetFromCmdLine();
	const uint FeatureCount = Params.GetFeatureCount();
	asserta(FeatureCount > 0);
	asserta(Params.m_Features[0] == FEATURE_AA);

	FILE *fOut = CreateStdioFile(opt_output);

	fprintf(fOut, "mega\t%u\t%u\n",
	  FeatureCount, ChainCount);

	for (uint i = 0; i < FeatureCount; ++i)
		{
		FEATURE F = Params.m_Features[i];
		uint AlphaSize = g_AlphaSizes2[F];
		fprintf(fOut, "%u\t%s\t%u\t%.6g\n",
		  i, FeatureToStr(F), AlphaSize, Params.m_Weights[i]);
		fprintf(fOut, "freqs");
		const float *Freqs = g_FeatureFreqsVec[F];
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			fprintf(fOut, "\t%.3f", Freqs[Letter]);
		fprintf(fOut, "\n");
		float **ScoreMx = Params.m_ScoreMxs[F];
		for (uint Letter1 = 0; Letter1 < AlphaSize; ++Letter1)
			{
			fprintf(fOut, "%u", Letter1);
			for (uint Letter2 = 0; Letter2 <= Letter1; ++Letter2)
				{
				float Score = ScoreMx[Letter1][Letter2];
				fprintf(fOut, "\t%.4f", Score);
				}
			fprintf(fOut, "\n");
			}
		}

	DSS D;
	D.m_Params = &Params;
	vector<vector<byte> > Profile;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		PDBChain &Chain = *Chains[ChainIndex];
		const uint L = Chain.GetSeqLength();
		const char *Label = Chain.m_Label.c_str();
		fprintf(fOut, "chain\t%u\t%s\t%u\n", ChainIndex, Label, L);

		D.Init(Chain);
		D.GetProfile(Profile);
		asserta(SIZE(Profile) == FeatureCount);
		for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
			asserta(SIZE(Profile[FeatureIdx]) == L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			string s;
			for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
				{
				byte Letter = Profile[FeatureIdx][Pos];
				if (FeatureIdx == 0)
					{
					byte c = g_LetterToCharAmino[Letter];
					if (c == INVALID_CHAR)
						c = 'X';
					s += c;
					}
				else
					s += 'A' + Letter;
				}
			fprintf(fOut, "%u\t%u\t%s\n", ChainIndex, Pos, s.c_str());
			}
		}

	CloseStdioFile(fOut);
	}
