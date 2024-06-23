#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"
#include "chainreader2.h"

static char GetChar(byte Letter, uint AlphaSize)
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

void cmd_pdb2alpha()
	{
	ChainReader2 CR;
	CR.Open(g_Arg1);

	DSSParams Params;
	Params.SetFromCmdLine(10000);

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

	uint AlphaSize = DSS::GetAlphaSize(Feat);

	FILE *fOut = CreateStdioFile(opt_output);

	uint Count = 0;
	DSS D;
	for (;;)
		{
		PDBChain *ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		PDBChain &Chain = *ptrChain;

		if (++Count%100 == 0)
			Progress("%u converted >%s\r",
			  Count, Chain.m_Label.c_str());
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		string Seq;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint Letter = D.GetFeature(Feat, Pos);
			char c = GetChar(Letter, AlphaSize);
			Seq += c;
			}
		SeqToFasta(fOut, Chain.m_Label, Seq);
		}
	Progress("%u chains converted\n", Count);
	CloseStdioFile(fOut);
	}
