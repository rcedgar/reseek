#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"
#include "chainreader2.h"
#include "bcadata.h"
#include <thread>

static FILE *s_fCal;
static FILE *s_fFasta;
static FILE *s_fFeatureFasta;

static PDBFileScanner *s_ptrFS = 0;
static time_t s_Now;
static time_t s_LastTime;
static uint s_TooShort;
static uint s_Converted;
static uint s_Shortest;
static uint s_InputCount;
static uint s_OutputCount;
static uint MinChainLength;
static BCAData *s_ptrBCA = 0;
static DSS *s_ptrD = 0;
static FEATURE s_Feat;
static mutex s_LockStats;
static mutex s_LockCal;
static mutex s_LockFasta;
static mutex s_LockBCA;
static set<string> *s_ptrLabelSet;
static uint s_LabelSetSize;

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

static void ThreadBody(uint ThreadIndex)
	{
	ChainReader2 CR;
	CR.Open(*s_ptrFS);
	if (optset_pdboutdir)
		CR.m_SaveLines = true;
	for (;;)
		{
		s_LockStats.lock();
		s_Now = time(0);
		if (s_Now - s_LastTime > 0)
			{
			if (s_LabelSetSize > 0)
				Progress("%u / %u chains found", s_OutputCount, s_LabelSetSize);
			else if (s_TooShort > 0)
				Progress("%s chains, %.1f%% too short (min %u, shortest %u)",
				  IntToStr(s_Converted), GetPct(s_TooShort, s_Converted),
				  MinChainLength, s_Shortest);
			else
				Progress("%s chains", IntToStr(s_Converted));
			uint ne = ChainReader2::m_CRGlobalFormatErrors;
			if (ne > 0)
				Progress(" %u format errors", ne);
			Progress("\r");
			s_LastTime = s_Now;
			}
		s_LockStats.unlock();

		PDBChain *ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		const uint L = ptrChain->GetSeqLength();
		asserta(L > 0);

		if (opt_reverse)
			ptrChain->Reverse();

		s_LockStats.lock();
		++s_InputCount;
		s_Shortest = min(L, s_Shortest);
		s_LockStats.unlock();

		if (s_ptrLabelSet != 0)
			{
			if (s_ptrLabelSet->find(ptrChain->m_Label) == s_ptrLabelSet->end())
				{
				delete ptrChain;
				continue;
				}
			}

		if (L < MinChainLength)
			{
			s_LockStats.lock();
			++s_TooShort;
			s_LockStats.unlock();
			delete ptrChain;
			continue;
			}

		s_LockStats.lock();
		++s_OutputCount;
		if (s_ptrLabelSet != 0)
			s_ptrLabelSet->erase(ptrChain->m_Label);
		s_LockStats.unlock();

		if (s_fCal != 0)
			{
			s_LockCal.lock();
			ptrChain->ToCal(s_fCal);
			s_LockCal.unlock();
			}

		if (s_fFasta != 0)
			{
			s_LockFasta.lock();
			ptrChain->ToFasta(s_fFasta);
			s_LockFasta.unlock();
			}

		if (s_fFeatureFasta != 0 && uint(s_Feat) < FEATURE_COUNT)
			{
			s_LockFasta.lock();
			ptrChain->ToFeatureFasta(s_fFeatureFasta, *s_ptrD, s_Feat);
			s_LockFasta.unlock();
			}

		if (s_ptrBCA != 0)
			{
			s_LockBCA.lock();
			s_ptrBCA->WriteChain(*ptrChain);
			s_LockBCA.unlock();
			}

		if (optset_pdboutdir)
			{
			const vector<string> &Lines = ptrChain->m_Lines;
			asserta(!Lines.empty());
			const string &FN = ptrChain->m_Label;
			string PathN = opt_pdboutdir;
			Dirize(PathN);
			PathN += FN;
			PathN += ".pdb";
			FILE *fOut = CreateStdioFile(PathN);
			for (uint i = 0; i < SIZE(Lines); ++i)
				{
				fputs(Lines[i].c_str(), fOut);
				fputc('\n', fOut);
				}
			CloseStdioFile(fOut);
			}

		s_LockStats.lock();
		++s_Converted;
		s_LockStats.unlock();

		delete ptrChain;
		}
	}

void cmd_convert()
	{
	PDBFileScanner FS;
	FS.Open(g_Arg1);

	uint MinChainLength = 1;
	if (optset_minchainlength)
		MinChainLength = opt_minchainlength;

	DSSParams Params;
	DSS D;
	uint AlphaSize = 0;
	s_Feat = FEATURE(FEATURE_COUNT);
	if (optset_feature_fasta)
		{
		Params.SetFromCmdLine(10000);
		s_Feat = GetFeatureFromCmdLine();
		}

	vector<string> Labels;
	if (optset_labels)
		{
		ReadLinesFromFile(opt_labels, Labels);
		if (Labels.empty())
			Die("No labels found in '%s'", opt_labels);
		}
	else if (optset_label)
		Labels.push_back(opt_label);

	set<string> LabelSet;
	for (uint i = 0; i < SIZE(Labels); ++i)
		LabelSet.insert(Labels[i]);
	if (!LabelSet.empty())
		{
		s_ptrLabelSet = &LabelSet;
		s_LabelSetSize = SIZE(LabelSet);
		}

	BCAData BCA;
	if (optset_bca)
		{
		BCA.Create(opt_bca);
		s_ptrBCA = &BCA;
		}

	s_ptrFS = &FS;
	s_ptrD = &D;

	s_fCal = CreateStdioFile(opt_cal);
	s_fFasta = CreateStdioFile(opt_fasta);
	s_fFeatureFasta = CreateStdioFile(opt_feature_fasta);

	s_InputCount = 0;
	s_Converted = 0;
	s_TooShort = 0;
	s_Shortest = UINT_MAX;
	s_LastTime = 0;
	vector<thread *> ts;

	uint ThreadCount = GetRequestedThreadCount();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	if (s_InputCount > 10000)
		{
		ProgressLog("\n");
		ProgressLogPrefix("%u / %u converted (%s, %.1f%%)\n",
		  s_OutputCount, s_InputCount, IntToStr(s_InputCount),
		  GetPct(s_OutputCount, s_InputCount));
		if (s_TooShort > 0)
			ProgressLogPrefix("%u too short (%s, %.1f%%) min length %u\n",
			  s_TooShort, IntToStr(s_TooShort), GetPct(s_TooShort, s_InputCount), MinChainLength);
		}
	else
		{
		ProgressLog("\n");
		ProgressLogPrefix("%u converted\n", s_InputCount);
		if (s_TooShort > 0)
			ProgressLogPrefix("%u too short (%.1f%%), min length %u, shortest %u\n",
			  s_TooShort, GetPct(s_TooShort, s_InputCount), MinChainLength, s_Shortest);
		}
	uint ne = ChainReader2::m_CRGlobalFormatErrors;
	if (ne > 0)
		ProgressLogPrefix("%u format errors\n", ne);

	if (optset_label || optset_labels)
		ProgressLogPrefix("Searched for %u labels, %u found (%.1f%%)\n", 
			s_LabelSetSize, s_OutputCount, GetPct(s_OutputCount, s_LabelSetSize));

	CloseStdioFile(s_fCal);
	CloseStdioFile(s_fFasta);
	CloseStdioFile(s_fFeatureFasta);
	if (optset_bca)
		BCA.Close();
	}
