#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"
#include "chainreader2.h"

static char GetFeatureChar(byte Letter, uint AlphaSize)
	{
	asserta(AlphaSize <= 36);
	if (Letter == UINT_MAX)
		return '*';
	if (Letter < 26)
		return 'A' + Letter;
	else if (Letter < 36)
		return 'a' + (Letter - 26);
	asserta(false);
	return '!';
	}

static FEATURE GetFeatureFromCmdLine()
	{
	string Alpha = "Mu";
	FEATURE Feat = FEATURE_Combo;
	if (optset_alpha)
		Alpha = opt_alpha;

#define c(x, y)	if (Alpha == #x) Alpha = #y;
	c(Conf3, SS3);
	c(Conf4, SS);
	c(Conf16, MySS);
	c(NENConf16, NbrMySS);
	c(RENConf16, RevNbrMySS);
	c(NENDist16, NbrDist);
	c(RENDist16, RevNbrDist);
#undef c

#define F(x) if (Alpha == #x) Feat = FEATURE_##x;
#include "intfeatures.h"
#undef F

	return Feat;
	}

void cmd_convert()
	{
	ChainReader2 CR;
	CR.Open(g_Arg1);

	FILE *fCal = CreateStdioFile(opt_cal);
	FILE *fFasta = CreateStdioFile(opt_fasta);
	FILE *fFeatureFasta = CreateStdioFile(opt_feature_fasta);

	DSSParams Params;
	DSS D;
	uint AlphaSize = 0;
	FEATURE Feat = FEATURE(0);
	if (fFeatureFasta != 0)
		{
		Params.SetFromCmdLine(10000);
		DSS::GetAlphaSize(Feat);
		Feat = GetFeatureFromCmdLine();
		}

	uint Count = 0;
	for (;;)
		{
		PDBChain *ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		PDBChain &Chain = *ptrChain;

		Chain.ToCal(fCal);
		Chain.ToFasta(fFasta);
		if (fFeatureFasta != 0)
			{
			D.Init(Chain);
			const uint L = Chain.GetSeqLength();
			string Seq;
			for (uint Pos = 0; Pos < L; ++Pos)
				{
				uint Letter = D.GetFeature(Feat, Pos);
				char c = GetFeatureChar(Letter, AlphaSize);
				Seq += c;
				}
			SeqToFasta(fFeatureFasta, Chain.m_Label, Seq);
			}
		delete ptrChain;
		}
	}
