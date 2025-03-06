#include "myutils.h"
#include "seqfeeder.h"

void SeqFeeder::Static_ThreadBody(SeqFeeder *SF)
	{
	SF->ThreadBody();
	}

void SeqFeeder::ThreadBody()
	{
	const auto FillWaitTime = std::chrono::milliseconds(m_FillWaitTime_ms);
	for (;;)
		{
		fprintf(stderr, "\nThreadBody_MuSeqSource\n");

		bool Any = false;
		for (;;)
			{
			for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
				{
				list<SeqInfo *> &SIList = m_SILists[ThreadIndex];
				size_t Size = SIList.size();
				if (Size >= m_MaxSIListSize)
					continue;

				mutex &ListLock = *m_ListLocks[ThreadIndex];
				mutex &OMLock = *m_OMLocks[ThreadIndex];
				ObjMgr &OM = *m_OMs[ThreadIndex];

				OMLock.lock();
				SeqInfo *SI = OM.GetSeqInfo();
				OMLock.unlock();

				bool Ok = m_SS->GetNext(SI);
				if (!Ok)
					{
					m_EOF = true;
					return;
					}
				Any = true;
				ListLock.lock();
				SIList.push_back(SI);
				ListLock.unlock();
				}
			}
		if (!Any)
			std::this_thread::sleep_for(FillWaitTime);
		}
	}

SeqFeeder::SeqFeeder(uint ThreadCount, SeqSource &SS)
	{
	m_ThreadCount = ThreadCount;
	m_SS = &SS;
	m_EOF = false;

	m_SILists.resize(ThreadCount);
	m_OMs.resize(ThreadCount);
	m_ListLocks.resize(ThreadCount);
	m_OMLocks.resize(ThreadCount);

	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		m_OMs[ThreadIndex] = new ObjMgr;
		m_ListLocks[ThreadIndex] = new mymutex;
		m_OMLocks[ThreadIndex] = new mymutex;
		}

	thread *t = new thread(Static_ThreadBody, this);
	}

void SeqFeeder::Down(uint ThreadIndex, SeqInfo *SI)
	{
	mutex &OMLock = *m_OMLocks[ThreadIndex];
	ObjMgr &OM = *m_OMs[ThreadIndex];
	OMLock.lock();
	OM.Down(SI);
	OMLock.unlock();
	}

SeqInfo *SeqFeeder::GetSI(uint ThreadIndex)
	{
	const auto GetWaitTime = std::chrono::milliseconds(m_GetWaitTime_ms);
	mutex &ListLock = *m_ListLocks[ThreadIndex];
	SeqInfo *SI = 0;
	for (;;)
		{
		list<SeqInfo *> &SIList = m_SILists[ThreadIndex];
		ListLock.lock();
		if (!SIList.empty())
			{
			SI = SIList.front();
			SIList.pop_front();
			}
		ListLock.unlock();
		if (SI != 0)
			return SI;

		if (m_EOF)
			return 0;

		std::this_thread::sleep_for(GetWaitTime);
		}
	}
