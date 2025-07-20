#include "myutils.h"
#include "pdbchain.h"
#include "bcadata.h"
#include <thread>

static BCAData *s_ptrBCAOut = 0;
static time_t s_Now;
static time_t s_LastTime;
static uint s_TooShort;
static uint s_Converted;
static uint s_Shortest;
static uint s_BCASize;
static uint s_ChainsPerSplit;
static uint s_NextSplitIdx;
static uint s_NextBCAIdx;
static uint s_OutputCount;
static uint s_MinChainLength;
static BCAData *s_ptrBCAIn = 0;
static mutex s_LockStats;
static mutex s_LockBCAOut;
static string s_SplitFN;

static void CreateSplit()
	{
	if (s_NextSplitIdx > 0)
		{
		asserta(s_ptrBCAOut != 0);
		s_ptrBCAOut->Close();
		delete s_ptrBCAOut;
		}
	else
		asserta(s_ptrBCAOut == 0);
	s_ptrBCAOut = new BCAData;
	s_SplitFN = opt(output);
	string sn;
	Ps(sn, "%u", s_NextSplitIdx+1);
	++s_NextSplitIdx;
	uint n = Replace(s_SplitFN, "@", sn);
	if (n == 0)
		Die("Missing @ in -output");
	s_ptrBCAOut->Create(s_SplitFN);
	}

static void ThreadBody(uint ThreadIndex)
	{
	for (;;)
		{
		uint MyIdx = UINT_MAX;
		s_LockStats.lock();
		s_Now = time(0);
		if (s_Now - s_LastTime > 0)
			{
			if (s_TooShort > 0)
				Progress("%s chains, %.1f%% too short (min %u, shortest %u) split %u",
				  IntToStr(s_Converted), GetPct(s_TooShort, s_Converted),
				  s_MinChainLength, s_Shortest,
				  s_NextSplitIdx);
			else
				Progress("%s chains, split %u", IntToStr(s_Converted), s_NextSplitIdx);
			Progress("\r");
			s_LastTime = s_Now;
			}
		if (s_NextBCAIdx < s_BCASize)
			{
			MyIdx = s_NextBCAIdx++;
			if (MyIdx >= s_NextSplitIdx*s_ChainsPerSplit)
				CreateSplit();
			}
		s_LockStats.unlock();
		if (MyIdx == UINT_MAX)
			return;
		asserta(MyIdx < s_BCASize);

		PDBChain Chain;
		s_ptrBCAIn->ReadChain(MyIdx, Chain);

		const uint L = Chain.GetSeqLength();
		asserta(L > 0);

		s_LockStats.lock();
		s_Shortest = min(L, s_Shortest);
		s_LockStats.unlock();

		if (L < s_MinChainLength)
			{
			s_LockStats.lock();
			++s_TooShort;
			s_LockStats.unlock();
			continue;
			}

		s_LockStats.lock();
		++s_OutputCount;
		s_LockStats.unlock();

		asserta(s_ptrBCAOut != 0);
		s_LockBCAOut.lock();
		s_ptrBCAOut->WriteChain(Chain);
		s_LockBCAOut.unlock();

		s_LockStats.lock();
		++s_Converted;
		s_LockStats.unlock();
		}
	}

void cmd_split()
	{
	if (optset_threads)
		Die("-threads not supported");
	asserta(optset_n);
	BCAData bca;
	bca.Open(g_Arg1);
	s_ptrBCAIn = &bca;
	s_BCASize = bca.GetChainCount();
	const uint SplitCount = opt(n);
	s_ChainsPerSplit = (s_BCASize + SplitCount - 1)/SplitCount;
	ProgressLog("%s chains/split\n", IntToStr(s_ChainsPerSplit));
	asserta(s_ChainsPerSplit*SplitCount >= s_BCASize);

	s_MinChainLength = 1;
	if (optset_minchainlength)
		s_MinChainLength = opt(minchainlength);

	s_Converted = 0;
	s_TooShort = 0;
	s_Shortest = UINT_MAX;
	s_LastTime = 0;

	ThreadBody(0);

	if (s_NextBCAIdx > 10000)
		{
		ProgressLog("\n");
		ProgressLogPrefix("%u / %u converted (%s, %.1f%%)\n",
		  s_OutputCount, s_NextBCAIdx, IntToStr(s_NextBCAIdx),
		  GetPct(s_OutputCount, s_NextBCAIdx));
		if (s_TooShort > 0)
			ProgressLogPrefix("%u too short (%s, %.1f%%) min length %u\n",
			  s_TooShort, IntToStr(s_TooShort), GetPct(s_TooShort, s_BCASize), s_MinChainLength);
		}
	else
		{
		ProgressLog("\n");
		ProgressLogPrefix("%u converted\n", s_NextBCAIdx);
		if (s_TooShort > 0)
			ProgressLogPrefix("%u too short (%.1f%%), min length %u, shortest %u\n",
			  s_TooShort, GetPct(s_TooShort, s_NextBCAIdx), s_MinChainLength, s_Shortest);
		}

	asserta(s_ptrBCAOut != 0);
	s_ptrBCAOut->Close();
	delete s_ptrBCAOut;
	s_ptrBCAIn->Close();
	}
