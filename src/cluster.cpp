#include "myutils.h"
#include "dbsearcher.h"
#include "chainreader2.h"
#include "mx.h"
#include "pdbchain.h"
#include "xdpmem.h"
#include "cigar.h"
#include "timing.h"
#include "sort.h"
#include <thread>

static void USort1(DBSearcher &DBS,  uint MinU,
  const vector<uint> &QueryKmerBits, vector<uint> &DBChainIndexes,
  vector<uint> &Order)
	{
	uint DBSize = DBS.GetDBSize();
	if (MinU == 0)
		{
		Order.clear();
		DBChainIndexes.clear();
		for (uint Idx = 0; Idx < DBSize; ++Idx)
			{
			Order.push_back(Idx);
			DBChainIndexes.push_back(Idx);
			}
		return;
		}

	StartTimer(USort);
	DBChainIndexes.clear();

	vector<uint> Us;
	for (uint Idx = 0; Idx < DBSize; ++Idx)
		{
		const vector<uint> &DBKmerBits = *DBS.m_DBKmerBitsVec[Idx];
		uint U = GetUBits(QueryKmerBits, DBKmerBits);
		if (U < MinU)
			continue;
		Us.push_back(U);
		DBChainIndexes.push_back(Idx);
		}
	EndTimer(USort);
	const uint N = SIZE(Us);
	if (N > 0)
		{
		Order.resize(N);
		QuickSortOrderDesc(Us.data(), N, Order.data());
		}
	else
		Order.clear();
	}

static uint ClusterQuery(DBSearcher &DBS, DSSAligner &DA,
  const PDBChain &ChainQ, const vector<vector<byte> > &ProfileQ,
  vector<byte> &ComboLettersQ, vector<uint> &KmerBitsQ)
	{
	const DSSParams &Params = *DBS.m_Params;
	const uint MinU = uint(Params.m_MinU + 0.5);
	vector<uint> DBIdxs;
	vector<uint> Order;
	USort1(DBS, MinU, KmerBitsQ, DBIdxs, Order);
	const uint N = SIZE(Order);
	if (N == 0)
		return UINT_MAX;
	const uint MAXACCEPTS = Params.m_MaxAccepts;
	const uint MAXREJECTS = Params.m_MaxRejects;
	uint AcceptCount = 0;
	uint RejectCount = 0;
	uint TopHit = UINT_MAX;
	double BestEvalue = DBL_MAX;
	for (uint i = 0; i < N; ++i)
		{
		if (AcceptCount >= MAXACCEPTS)
			break;
		if (RejectCount >= MAXREJECTS)
			break;
		uint ChainIndexR = DBIdxs[Order[i]];
		const PDBChain &ChainR = *DBS.m_DBChains[ChainIndexR];
		const vector<byte> &ComboLettersR = *DBS.m_DBComboLettersVec[ChainIndexR];
		const vector<vector<byte> > &ProfileR = *DBS.m_DBProfiles[ChainIndexR];
		DA.Align_ComboFilter(ChainQ, ChainR, 
			ComboLettersQ, ComboLettersR, ProfileQ, ProfileR);
		if (DA.m_Path.empty())
			continue;
		if (DA.m_EvalueA > DBS.m_MaxEvalue)
			{
			++RejectCount;
			continue;
			}
		++AcceptCount;
		double Evalue = DA.GetEvalue(true);
		if (Evalue <= DBS.m_MaxEvalue)
			{
			if (TopHit == UINT_MAX || Evalue < BestEvalue)
				{
				TopHit = ChainIndexR;
				BestEvalue = Evalue;
				}
			}
		}
	if (DBS.m_fTsv != 0)
		{
		if (TopHit == UINT_MAX)
			{
			uint Idx = DBS.GetDBChainCount();
			const char *Label = ChainQ.m_Label.c_str();
			fprintf(DBS.m_fTsv, "C\t%u\t%s\n", Idx, Label);
			}
		else
			{
			const PDBChain &ChainT = *DBS.m_DBChains[TopHit];
			fprintf(DBS.m_fTsv, "M\t%u\t%s\t%s\t%.3g\n",
			  TopHit, ChainQ.m_Label.c_str(), ChainT.m_Label.c_str(), BestEvalue);
			}
		}
	return TopHit;
	}

struct ThreadUserData
	{
	atomic<uint> *ptrInputCount = 0;
	DBSearcher *ptrDBS = 0;
	ChainReader2 *ptrCR = 0;
	};

static void ThreadBody(uint ThreadIndex, void *ptrUserData)
	{
	static atomic<time_t> LastTime;
	static mutex Lock;

	ThreadUserData *ptrUD  = (ThreadUserData *) ptrUserData;
	DBSearcher &DBS = *ptrUD->ptrDBS;
	ChainReader2 &CR = *ptrUD->ptrCR;
	atomic<uint> &InputCount = *ptrUD->ptrInputCount;

	const DSSParams &Params = *DBS.m_Params;
	DSS D;
	D.SetParams(Params);

	DSSAligner DA;
	DA.m_Params = &Params;


	vector<byte> ComboLettersQ;
	vector<uint> KmersQ;
	vector<uint> KmerBitsQ;
	uint ClusterCount = 0;
	vector<vector<byte> > ProfileQ;
	for (;;)
		{
		PDBChain *ptrChainQ = CR.GetNext();
		if (ptrChainQ == 0)
			break;
		Lock.lock();
		time_t Now = time(0);
		if (Now - LastTime > 0)
			{
			ClusterCount = DBS.GetDBChainCount();
			Progress("%s chains, %u clusters\r",
			  IntToStr(InputCount), ClusterCount);
			LastTime = Now;
			}
		Lock.unlock();

		D.Init(*ptrChainQ);
		D.GetProfile(ProfileQ);
		D.GetComboLetters(ComboLettersQ);
		D.GetComboKmers(ComboLettersQ, KmersQ);
		D.GetComboKmerBits(KmersQ, KmerBitsQ);

		++InputCount;
		uint TopHit = 
		  ClusterQuery(DBS, DA, *ptrChainQ, ProfileQ, ComboLettersQ, KmerBitsQ);
		if (TopHit == UINT_MAX)
			{
			vector<vector<byte> > *ptrProfileQ = new vector<vector<byte> >(ProfileQ);
			vector<byte> *ptrComboLettersQ = new vector<byte>(ComboLettersQ);
			vector<uint> *ptrKmerBitsQ = new vector<uint>(KmerBitsQ);
			DBS.AddChain(ptrChainQ, ptrProfileQ, ptrComboLettersQ, ptrKmerBitsQ);
			}
		}
	}

void cmd_cluster()
	{
	const string &QFN = g_Arg1;

	DBSearcher DBS;
	DSSParams Params;
	Params.SetFromCmdLine(10000);

#define p(x, y, d)	{ if (optset_##x) Params.m_##y = opt_##x; else Params.m_##y = d; }
	p(maxaccepts, MaxAccepts, 1)
	p(maxrejects, MaxRejects, 128)
	p(minu, MinU, 0)
#undef p

	DBS.m_Params = &Params;

	DBS.InitEmpty();
	if (optset_evalue)
		DBS.m_MaxEvalue = (float) opt_evalue;
	else
		DBS.m_MaxEvalue = 10;

	DBS.m_fTsv = CreateStdioFile(opt_output);

	ChainReader2 CR;
	CR.Open(QFN);

	atomic<uint> InputCount;
	ThreadUserData UD;
	UD.ptrInputCount = &InputCount;
	UD.ptrDBS = &DBS;
	UD.ptrCR = &CR;
	RunThreads(ThreadBody, (void *) &UD);

	uint ClusterCount = DBS.GetDBChainCount();
	ProgressLog("%s chains, %u clusters\n", IntToStr(InputCount), ClusterCount);

	CloseStdioFile(DBS.m_fTsv);
	}
