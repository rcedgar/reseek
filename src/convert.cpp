#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"
#include "chainreader2.h"
#include "bcadata.h"
#include <thread>

static FILE *s_fCal;
static FILE *s_fFasta;
static FILE *s_fMultiPDB;
static FILE *s_fFeatureFasta;

static PDBFileScanner *s_ptrFS = 0;
static time_t s_Now;
static time_t s_LastTime;
static uint s_TooShort;
static uint s_Converted;
static uint s_Shortest;
static uint s_InputCount;
static uint s_OutputCount;
static uint s_MinChainLength;
static BCAData *s_ptrBCA = 0;
static FEATURE s_Feat;
static mutex s_LockStats;
static mutex s_LockCal;
static mutex s_LockMultiPDB;
static uint s_MultiPDBIdx;
static mutex s_LockFasta;
static mutex s_LockBCA;
static set<string> *s_ptrLabelSet;
static uint s_LabelSetSize;

static FEATURE GetFeatureFromCmdLine()
	{
	string Alpha = "Mu";
	FEATURE Feat = FEATURE(FEATURE_COUNT);
	if (optset_alpha)
		Alpha = opt(alpha);
	if (Alpha == "NENConf3")
		Alpha = "NENSS3";
	else if (Alpha == "RENDist4")
		Alpha = "RENDist4";

	//MuFeatures.push_back(FEATURE_SS3);
	//MuFeatures.push_back(FEATURE_NENSS3);
	//MuFeatures.push_back(FEATURE_RENDist4);
#define c(x, y)	if (stricmp(Alpha.c_str(), #x) == 0) Alpha = #y;
	c(Mu, Mu);
	c(Conf3, SS3);
	c(NENSS3, NENSS3);
	c(RENDist4, RENDist4);
	c(Conf4, SS);
	c(Conf16, Conf);
	c(NENConf16, NENConf);
	c(RENConf16, RENConf);
	c(NENDist16, NENDist);
	c(RENDist16, RENDist);
	c(RENDist4, RENDist4);
	c(NormDens16, NormDens);
	c(StrandDens16, StrandDens);
#undef c

#define F(x) if (stricmp(Alpha.c_str(), #x) == 0) Feat = FEATURE_##x;
#include "intfeatures.h"
#undef F

	if (Feat == FEATURE(FEATURE_COUNT))
		Die("Invalid -alpha %s", opt(alpha));

	return Feat;
	}

static void ThreadBody(uint ThreadIndex)
	{
	DSS *ptrD = 0;
	if (s_fFeatureFasta != 0 && uint(s_Feat) < FEATURE_COUNT)
		ptrD = new DSS;

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
				Progress("%u / %u chains found (%s searched)",
						 s_OutputCount, s_LabelSetSize, IntToStr(s_InputCount));
			else if (s_TooShort > 0)
				Progress("%s chains, %.1f%% too short (min %u, shortest %u)",
				  IntToStr(s_Converted), GetPct(s_TooShort, s_Converted),
				  s_MinChainLength, s_Shortest);
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

		if (opt(reverse))
			ptrChain->Reverse();
		if (opt(flip))
			ptrChain->Flip();

		s_LockStats.lock();
		++s_InputCount;
		s_Shortest = min(L, s_Shortest);
		s_LockStats.unlock();

		if (s_ptrLabelSet != 0)
			{
			string UpperLabel = ptrChain->m_Label;
			ToUpper(UpperLabel);
			set<string>::iterator iter = s_ptrLabelSet->find(UpperLabel);
			if (iter == s_ptrLabelSet->end())
				{
				delete ptrChain;
				continue;
				}
			}

		if (L < s_MinChainLength)
			{
			s_LockStats.lock();
			++s_TooShort;
			s_LockStats.unlock();
			delete ptrChain;
			continue;
			}

		if (optset_subsample)
			{
			s_LockStats.lock();
			uint n = s_InputCount;
			s_LockStats.unlock();
			if (n%opt(subsample) != 0)
				{
				delete ptrChain;
				continue;
				}
			}

		s_LockStats.lock();
		++s_OutputCount;
		if (s_ptrLabelSet != 0)
			{
			string UpperLabel = ptrChain->m_Label;
			ToUpper(UpperLabel);
			s_ptrLabelSet->erase(UpperLabel);
			}
		s_LockStats.unlock();

		if (s_fMultiPDB != 0)
			{
			s_LockMultiPDB.lock();
			fprintf(s_fMultiPDB, "MODEL%10u\n", s_MultiPDBIdx);
			if (ptrChain->m_Label != "")
				fprintf(s_fMultiPDB, "TITLE     %s\n", ptrChain->m_Label.c_str());
			else
				fprintf(s_fMultiPDB, "TITLE     _blank_%u\n", s_MultiPDBIdx);
			ptrChain->ToPDB(s_fMultiPDB, true);
			fprintf(s_fMultiPDB, "ENDMDL\n");
			++s_MultiPDBIdx;
			s_LockMultiPDB.unlock();
			}

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
			DSS &D = *ptrD;
			D.Init(*ptrChain);
			const uint L = ptrChain->GetSeqLength();
			const uint AlphaSize = DSSParams::GetAlphaSize(s_Feat);
			string Seq;
			for (uint Pos = 0; Pos < L; ++Pos)
				{
				uint Letter = D.GetFeature(s_Feat, Pos);
				char GetFeatureChar(byte Letter, uint AlphaSize);
				char c = GetFeatureChar(Letter, AlphaSize);
				Seq += c;
				}
			s_LockFasta.lock();
			//ptrChain->ToFeatureFasta(s_fFeatureFasta, *s_ptrD, s_Feat);
			SeqToFasta(s_fFeatureFasta, ptrChain->m_Label, Seq);
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
			string PathN = opt(pdboutdir);
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

		if (optset_pdbcaoutdir)
			{
			string &FN = ptrChain->m_Label;
			string PathN = opt(pdbcaoutdir);
			Dirize(PathN);
			PathN += FN;
			PathN += ".pdb";
			ptrChain->ToPDB(PathN);
			}

		s_LockStats.lock();
		++s_Converted;
		s_LockStats.unlock();

		delete ptrChain;
		}
	}
	
void cmd_subset()
	{
	Die("Use -convert -subsample");
	}

void cmd_convert()
	{
	PDBFileScanner FS;
	FS.Open(g_Arg1);

	optset_fast = true;
	opt(fast) = true;

	s_MinChainLength = 1;
	if (optset_minchainlength)
		s_MinChainLength = opt(minchainlength);

	uint AlphaSize = 0;
	s_Feat = FEATURE(FEATURE_COUNT);
	if (optset_feature_fasta)
		{
		optset_fast = true;
		opt(fast) = true;
		DSSParams::Init(DM_AlwaysFast);
		s_Feat = GetFeatureFromCmdLine();
		}

	vector<string> Labels;
	if (optset_labels)
		{
		ReadLinesFromFile(opt(labels), Labels);
		if (Labels.empty())
			Die("No labels found in '%s'", opt(labels));
		}
	else if (optset_label)
		Labels.push_back(opt(label));

	set<string> LabelSet;
	for (uint i = 0; i < SIZE(Labels); ++i)
		{
		ToUpper(Labels[i]);
		LabelSet.insert(Labels[i]);
		}
	if (!LabelSet.empty())
		{
		s_ptrLabelSet = &LabelSet;
		s_LabelSetSize = SIZE(LabelSet);
		}

	BCAData BCA;
	if (optset_bca)
		{
		BCA.Create(opt(bca));
		s_ptrBCA = &BCA;
		}

	s_ptrFS = &FS;

	s_fCal = CreateStdioFile(opt(cal));
	s_fFasta = CreateStdioFile(opt(fasta));
	s_fFeatureFasta = CreateStdioFile(opt(feature_fasta));
	s_fMultiPDB = CreateStdioFile(opt(multipdb));

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
			  s_TooShort, IntToStr(s_TooShort), GetPct(s_TooShort, s_InputCount), s_MinChainLength);
		}
	else
		{
		ProgressLog("\n");
		ProgressLogPrefix("%u converted\n", s_InputCount);
		if (s_TooShort > 0)
			ProgressLogPrefix("%u too short (%.1f%%), min length %u, shortest %u\n",
			  s_TooShort, GetPct(s_TooShort, s_InputCount), s_MinChainLength, s_Shortest);
		}
	uint ne = ChainReader2::m_CRGlobalFormatErrors;
	if (ne > 0)
		ProgressLogPrefix("%u format errors\n", ne);

	if (optset_label || optset_labels)
		{
		ProgressLogPrefix("Searched for %u labels, %u found (%.1f%%)\n", 
			s_LabelSetSize, s_OutputCount, GetPct(s_OutputCount, s_LabelSetSize));
		uint NotFound = SIZE(*s_ptrLabelSet);
		if (NotFound > 0)
			{
			Progress("%u not found", NotFound);
			Log("%u not found\n", NotFound);
			uint Counter = 0;
			for (set<string>::const_iterator iter = s_ptrLabelSet->begin();
				 iter != s_ptrLabelSet->end(); ++iter)
				{
				if (Counter++ < 3)
					Progress(" %s", iter->c_str());
				Log(">%s\n", iter->c_str());
				}
			Progress("\n");
			}
		}

	CloseStdioFile(s_fCal);
	CloseStdioFile(s_fFasta);
	CloseStdioFile(s_fFeatureFasta);
	CloseStdioFile(s_fMultiPDB);
	if (optset_bca)
		BCA.Close();
	}
