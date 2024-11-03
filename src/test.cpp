#include "myutils.h"
#include "chainreader2.h"
#include "dss.h"
#include "dssaligner.h"

float GetSelfRevScore(DSSAligner &DA, DSS &D, 
  const PDBChain &Chain, const vector<vector<byte> > &Profile);

static ChainReader2 CR;
static DSSParams Params;
static uint g_Count = 0;

static void BodyVec()
	{
	DSSParams DA_Params = Params;
	DA_Params.m_UsePara = false;
	DA_Params.m_Omega = 0;

	vector<DSS *> Ds;
	vector<DSSAligner *> DAs;
	for (uint i = 0; i < 32; ++i)
		{
		DSS *ptrD = new DSS;
		DSSAligner *ptrDA = new DSSAligner;
		ptrD->SetParams(Params);
		ptrDA->m_Params = &DA_Params;
		Ds.push_back(ptrD);
		DAs.push_back(ptrDA);
		}

	time_t Last = time(0);
	uint Counter = 0;
	for (;;)
		{
		++g_Count;
		PDBChain *Chain = CR.GetNext();
		if (Chain == 0)
			break;
		time_t Now = time(0);
		if (Now - Last > 0)
			{
			Last = Now;
			static mutex m;
			m.lock();
			Progress("Testing %s %s\r", IntToStr(g_Count),
					MemBytesToStr(GetMemUseBytes()));
			m.unlock();
			}
		vector<vector<byte> > Profile;
		DSS &D = *Ds[++Counter%32];
		DSSAligner &DA = *DAs[++Counter%32];
		D.Init(*Chain);
		D.GetProfile(Profile);
		float SelfRevScore = GetSelfRevScore(DA, D, *Chain, Profile);
		}
	}

static void ThreadBody(uint ThreadIndex)
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

void cmd_test()
	{
	CR.Open(g_Arg1);
	Params.SetFromCmdLine(10000);
	BodyVec();
	return;
	const uint ThreadCount = GetRequestedThreadCount();

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();

	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	}
