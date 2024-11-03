#include "myutils.h"
#include "chainreader2.h"
#include "dss.h"
#include "dssaligner.h"
#include "profileloader.h"

float GetSelfRevScore(DSSAligner &DA, DSS &D, 
  const PDBChain &Chain, const vector<vector<byte> > &Profile);

static ChainReader2 CR;
static DSSParams Params;
static uint g_Count = 0;

static void Test2ThreadBody(uint ThreadIndex)
	{
	DSS D;
	D.SetParams(Params);

	DSSAligner DA;
	DSSParams DA_Params = Params;
	DA_Params.m_UsePara = false;
	DA_Params.m_Omega = 0;
	DA.m_Params = &DA_Params;

	time_t Last = time(0);
	for (;;)
		{
		++g_Count;
		PDBChain *Chain = CR.GetNext();
		if (Chain == 0)
			break;
		time_t Now = time(0);
		if (ThreadIndex == 0 && Now - Last > 0)
			{
			Last = Now;
			static mutex m;
			m.lock();
			Progress("Testing %s %s\r", IntToStr(g_Count),
					MemBytesToStr(GetMemUseBytes()));
			m.unlock();
			}
		vector<vector<byte> > Profile;
		D.Init(*Chain);
		D.GetProfile(Profile);
		float SelfRevScore = GetSelfRevScore(DA, D, *Chain, Profile);
		}
	}

static void Test2()
	{
	CR.Open(g_Arg1);
	Params.SetFromCmdLine(10000);
	const uint ThreadCount = GetRequestedThreadCount();

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(Test2ThreadBody, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();

	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	}

static void Test1()
	{
	ProfileLoader PL;
	ChainReader2 CR;
	CR.Open(g_Arg1);

	DSSParams Params;
	Params.SetFromCmdLine(10000);

	static vector<vector<byte> *> m_ComboLettersVec;
	static vector<vector<uint> *> m_KmerBitsVec;
	static vector<PDBChain *> m_DBChains;
	static vector<vector<vector<byte> > *> m_Profiles;
	static vector<float> m_DBSelfRevScores;

	uint ThreadCount = GetRequestedThreadCount();

	vector<vector<byte> *> *ptrComboLetters = &m_ComboLettersVec;
	vector<vector<uint> *> *ptrKmerBitsVec = &m_KmerBitsVec;
	vector<float> *ptrSelfRevScores = &m_DBSelfRevScores;
	//ptrSelfRevScores = 0;

	if (Params.m_MinU <= 0)
		ptrKmerBitsVec = 0;
	if (Params.m_Omega <= 0)
		ptrComboLetters = 0;

	ProgressLogPrefix("Before Load %.3g cap=%.3g\n",
	  GetMemUseBytes(), double(m_DBSelfRevScores.capacity()));
	PL.Load(Params, CR, &m_DBChains, &m_Profiles, ptrComboLetters,
	  ptrKmerBitsVec, ptrSelfRevScores, ThreadCount);
	ProgressLogPrefix("Before return %.3g cap=%.3g\n",
	  GetMemUseBytes(), double(m_DBSelfRevScores.capacity()));
	}

void cmd_test()
	{
	switch (opt_n)
		{
	case 1:	Test1(); break;
	case 2:	Test2(); break;
	default: asserta(false);
		}
	ProgressLogPrefix("After return %.3g\n", GetMemUseBytes());
	}
