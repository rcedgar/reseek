#include "myutils.h"
#include "sort.h"
#include "dbsearcher.h"
#include "chainreader2.h"
#include <thread>

uint GetUBits(const vector<uint> &KmerBitsQ, const vector<uint> &KmerBitsR);

void DBSearcher::USort(const vector<uint> &QueryKmerBits,
  vector<uint> &DBChainIndexes, vector<uint> &Order)
	{
	StartTimer(USort);
	DBChainIndexes.clear();

	const uint MinU = uint(round(m_Params->m_MinU));
	vector<uint> Us;
	uint DBSize = GetDBSize();
	for (uint Idx = 0; Idx < DBSize; ++Idx)
		{
		const vector<uint> &DBKmerBits = *m_DBKmerBitsVec[Idx];
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

void DBSearcher::StaticThreadUSort(uint ThreadIndex,
  DBSearcher *ptrDBS, ChainReader2 *ptrQCR)
	{
	ptrDBS->ThreadUSort(ThreadIndex, *ptrQCR);
	}

void DBSearcher::ThreadUSort(uint ThreadIndex, ChainReader2 &QCR)
	{
	asserta(ThreadIndex < SIZE(m_DAs));
	DSSAligner &DA = *m_DAs[ThreadIndex];
	DA.SetParams(*m_Params);
	DSS D;
	D.SetParams(*m_Params);
	vector<uint> Order;
	vector<byte> MuLettersQ;
	vector<uint> KmersQ;
	vector<uint> KmerBitsQ;
	vector<vector<byte> > ProfileQ;
	for (;;)
		{
		m_Lock.lock();
		const PDBChain *ptrChainQ = QCR.GetNext();
		if (ptrChainQ == 0)
			{
			m_Lock.unlock();
			break;
			}
		m_Lock.unlock();
		const PDBChain &ChainQ = *ptrChainQ;
		D.Init(ChainQ);
		D.GetMuLetters(MuLettersQ);
		D.GetMuKmers(MuLettersQ, KmersQ);
		D.GetMuKmerBits(KmersQ, KmerBitsQ);
		vector<uint> DBIdxs;
		USort(KmerBitsQ, DBIdxs, Order);
		const uint N = SIZE(Order);
		if (N == 0)
			continue;
		const uint MAXACCEPTS = 1;//@@TODO
		const uint MAXREJECTS = 64;//@@TODO
		uint AcceptCount = 0;
		uint RejectCount = 0;
		D.GetProfile(ProfileQ);
		for (uint i = 0; i < N; ++i)
			{
			if (AcceptCount >= MAXACCEPTS)
				break;
			if (RejectCount >= MAXREJECTS)
				break;
			uint ChainIndexR = DBIdxs[Order[i]];
			const PDBChain &ChainR = *m_DBChains[ChainIndexR];
			if (ChainQ.m_Label == ChainR.m_Label)
				continue;
			const vector<byte> &MuLettersR = *m_DBMuLettersVec[ChainIndexR];
			const vector<uint> &KmersR = *m_DBMuKmersVec[ChainIndexR];
			const vector<vector<byte> > &ProfileR = *m_DBProfiles[ChainIndexR];
			DA.Align_MuFilter(ChainQ, ChainR, 
			  MuLettersQ, KmersQ, MuLettersR, KmersR, ProfileQ, ProfileR);
			if (DA.m_Path.empty())
				continue;
			if (DA.m_EvalueA > m_MaxEvalue)
				{
				++RejectCount;
				continue;
				}
			++AcceptCount;
			if (DA.GetEvalue(true) <= m_MaxEvalue)
				{
				DA.ToTsv(m_fTsv, true);
				BaseOnAln(DA, true);
				}
			}
		}
	}
