#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"
#include "chainreader2.h"
#include "bcadata.h"

static FILE *s_fCal;
static FILE *s_fFasta;
static FILE *s_fFeatureFasta;

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

	uint MinChainLength = 0;
	if (optset_minchainlength)
		MinChainLength = opt_minchainlength;

	DSSParams Params;
	DSS D;
	uint AlphaSize = 0;
	FEATURE Feat = FEATURE(0);
	if (s_fFeatureFasta != 0)
		{
		Params.SetFromCmdLine(10000);
		Feat = GetFeatureFromCmdLine();
		}

	BCAData BCA;
	if (optset_bca)
		BCA.Create(opt_bca);

	s_fCal = CreateStdioFile(opt_cal);
	s_fFasta = CreateStdioFile(opt_fasta);
	s_fFeatureFasta = CreateStdioFile(opt_feature_fasta);

	set<string> Labels;
	uint Count = 0;
	uint DupeLabelCount = 0;
	time_t LastTime = 0;

	for (;;)
		{
		PDBChain *ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		PDBChain &Chain = *ptrChain;
		if (Chain.GetSeqLength() == 0)
			{
			delete ptrChain;
			continue;
			}
		if (Labels.find(Chain.m_Label) != Labels.end())
			{
			Log("Dupe >%s\n", Chain.m_Label.c_str());
			++DupeLabelCount;
			delete ptrChain;
			continue;
			}
		time_t Now = time(0);
		++Count;
		if (Now - LastTime > 0)
			{
			Progress("%u chains converted, %u dupe labels\r",
			  Count, DupeLabelCount);
			LastTime = Now;
			}
		Labels.insert(Chain.m_Label);
		Chain.ToCal(s_fCal);
		Chain.ToFasta(s_fFasta);
		Chain.ToFeatureFasta(s_fFeatureFasta, D, Feat);
		if (optset_bca)
			BCA.WriteChain(Chain);
		delete ptrChain;
		}

	Log("%u dupe labels\n", DupeLabelCount);
	ProgressLog("\n%u chains converted\n", Count);

	CloseStdioFile(s_fCal);
	CloseStdioFile(s_fFasta);
	CloseStdioFile(s_fFeatureFasta);
	if (optset_bca)
		BCA.Close();
	}
