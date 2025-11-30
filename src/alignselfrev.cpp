#include "myutils.h"
#include "dssaligner.h"
#include "chainreader2.h"

void cmd_alignselfrev()
	{
	const string &QFN = g_Arg1;
	FILE *fOut = CreateStdioFile(opt(output));

	DSSParams::Init(DM_AlwaysSensitive);
	DSSParams::m_Omega8 = 0;
	DSSParams::m_Omega16 = 0;

	DSSAligner DA;
	DSS D;

	vector<vector<byte> > Profile;
	vector<vector<byte> > RevProfile;
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
		PDBChain Rev = *Chain;
		Rev.Reverse();

		D.Init(*Chain);
		D.GetProfile(Profile);

		D.Init(Rev);
		D.GetProfile(RevProfile);

		DA.SetQuery(*Chain, &Profile, 0, 0, FLT_MAX);
		DA.SetTarget(Rev, &RevProfile, 0, 0, FLT_MAX);
		DA.AlignQueryTarget();
		DA.ToTsv(fOut, true);
		}
	CloseStdioFile(fOut);
	}
