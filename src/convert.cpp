#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"
#include "chainreader2.h"

static FILE *s_fCal;
static FILE *s_fFasta;
static FILE *s_fFeatureFasta;

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

static void ThreadBody(
  uint ThreadIndex,
  ChainReader2 *ptrCR,
  set<string> *ptrLabels,
  time_t *ptrLastTime,
  uint *ptrCount,
  uint *ptrDupeLabelCount,
  mutex *ptrLock)
	{
// Using reference arguments to ThreadBody gives compile errors
	ChainReader2 &CR = *ptrCR;
	set<string> &Labels = *ptrLabels;
	time_t &LastTime = *ptrLastTime;
	uint &Count = *ptrCount;
	uint &DupeLabelCount = *ptrDupeLabelCount;
	mutex &Lock = *ptrLock;

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
		DSS::GetAlphaSize(Feat);
		Feat = GetFeatureFromCmdLine();
		}

	for (;;)
		{
		PDBChain *ptrChain = ptrCR->GetNext();
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
			Lock.lock();
			Log("Dupe >%s\n", Chain.m_Label.c_str());
			++DupeLabelCount;
			Lock.unlock();
			delete ptrChain;
			continue;
			}
		time_t Now = time(0);
		Lock.lock();
		++Count;
		if (Now - LastTime > 0)
			{
			Progress("%u chains converted, %u dupe labels, %u dupe seqs\r",
			  Count, DupeLabelCount, ChainReader2::m_DupeSeqCount);
			LastTime = Now;
			}
		Labels.insert(Chain.m_Label);
		Lock.unlock();

		Lock.lock();
		Chain.ToCal(s_fCal);
		Lock.unlock();

		Lock.lock();
		Chain.ToFasta(s_fFasta);
		Lock.unlock();

		if (s_fFeatureFasta != 0)
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
			Lock.lock();
			SeqToFasta(s_fFeatureFasta, Chain.m_Label, Seq);
			Lock.unlock();
			}
		delete ptrChain;
		}
	}

void cmd_convert()
	{
	PDBFileScanner FS;
	FS.Open(g_Arg1);

	ChainReader2 CR;
	CR.Open(FS);

	s_fCal = CreateStdioFile(opt_cal);
	s_fFasta = CreateStdioFile(opt_fasta);
	s_fFeatureFasta = CreateStdioFile(opt_feature_fasta);

	set<string> Labels;
	uint Count = 0;
	uint DupeLabelCount = 0;
	time_t LastTime = 0;
	mutex Lock;

	ThreadBody(
			0, 
			&CR,
			&Labels,
			&LastTime,
			&Count,
			&DupeLabelCount,
			&Lock);

	Log("%u dupe labels\n", DupeLabelCount);
	ProgressLog("\n%u chains converted\n", Count);
	CloseStdioFile(s_fCal);
	CloseStdioFile(s_fFasta);
	CloseStdioFile(s_fFeatureFasta);
	}
