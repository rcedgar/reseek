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
	FEATURE Feat = FEATURE(FEATURE_COUNT);
	if (optset_alpha)
		Alpha = opt_alpha;
	if (Alpha == "NENConf3")
		Alpha = "NbrSS3";
	else if (Alpha == "RENDist4")
		Alpha = "RevNbrDist4";

	//ComboFeatures.push_back(FEATURE_SS3);
	//ComboFeatures.push_back(FEATURE_NbrSS3);
	//ComboFeatures.push_back(FEATURE_RevNbrDist4);
#define c(x, y)	if (stricmp(Alpha.c_str(), #x) == 0) Alpha = #y;
	c(Mu, COMBO);
	c(Conf3, SS3);
	c(NbrSS3, NbrSS3);
	c(RevNbrDist4, RevNbrDist4);
	c(Conf4, SS);
	c(Conf16, MySS);
	c(NENConf16, NbrMySS);
	c(RENConf16, RevNbrMySS);
	c(NENDist16, NbrDist);
	c(RENDist16, RevNbrDist);
	c(RENDist4, RevNbrDist4);
	c(NormDens16, NormDens);
	c(StrandDens16, StrandDens);
#undef c

#define F(x) if (stricmp(Alpha.c_str(), #x) == 0) Feat = FEATURE_##x;
#include "intfeatures.h"
#undef F

	if (Feat == FEATURE(FEATURE_COUNT))
		Die("Invalid -alpha %s", opt_alpha);

	return Feat;
	}

void cmd_convert()
	{
	ChainReader2 CR;
	CR.Open(g_Arg1);

	uint MinChainLength = 1;
	if (optset_minchainlength)
		MinChainLength = opt_minchainlength;

	DSSParams Params;
	DSS D;
	uint AlphaSize = 0;
	FEATURE Feat = FEATURE(FEATURE_COUNT);
	if (optset_feature_fasta)
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

	uint InputCount = 0;
	uint Converted = 0;
	uint TooShort = 0;
	uint Shortest = UINT_MAX;
	time_t LastTime = 0;
	for (;;)
		{
		PDBChain *ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		if (opt_reverse)
			ptrChain->Reverse();
		time_t Now = time(0);
		if (Now - LastTime > 0)
			{
			if (TooShort > 0)
				Progress("%s chains, %.1f%% too short (min %u, shortest %u)\r",
				  IntToStr(Converted), GetPct(TooShort, Converted), MinChainLength, Shortest);
			else
				Progress("%s chains\r", IntToStr(Converted));
			LastTime = Now;
			}

		const uint L = ptrChain->GetSeqLength();
		if (L == 0)
			{
			delete ptrChain;
			continue;
			}

		++InputCount;
		Shortest = min(L, Shortest);
		if (L < MinChainLength)
			{
			++TooShort;
			delete ptrChain;
			continue;
			}

		ptrChain->ToCal(s_fCal);
		ptrChain->ToFasta(s_fFasta);
		ptrChain->ToFeatureFasta(s_fFeatureFasta, D, Feat);
		if (optset_bca)
			BCA.WriteChain(*ptrChain);
		++Converted;
		delete ptrChain;
		}

	if (InputCount > 10000)
		{
		ProgressLog("\n");
		ProgressLog("%10u converted (%s)\n", InputCount, IntToStr(InputCount));
		if (TooShort > 0)
			ProgressLog("%10u too short (%s, %.1f%%) min length %u\n",
			  TooShort, IntToStr(TooShort), GetPct(TooShort, InputCount), MinChainLength);
		}
	else
		{
		ProgressLog("\n");
		ProgressLog("%10u converted\n", InputCount);
		if (TooShort > 0)
			ProgressLog("%10u too short (%.1f%%), min length %u, shortest %u\n",
			  TooShort, GetPct(TooShort, InputCount), MinChainLength, Shortest);
		}

	CloseStdioFile(s_fCal);
	CloseStdioFile(s_fFasta);
	CloseStdioFile(s_fFeatureFasta);
	if (optset_bca)
		BCA.Close();
	}
