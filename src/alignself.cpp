#include "myutils.h"
#include "dssaligner.h"
#include "chainreader2.h"

void cmd_alignself()
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
		D.GetProfile(Profile);
		DA.SetQuery(*Chain, &Profile, 0, 0, FLT_MAX);
		DA.SetTarget(*Chain, &Profile, 0, 0, FLT_MAX);
		DA.AlignQueryTarget();
		DA.WriteDPMx(fOut, true);
		}
	CloseStdioFile(fOut);
	}
