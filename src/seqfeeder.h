#pragma once

#include <list>
#include "mymutex.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"

#define SEQ_FEEDER_STATS	1

#if SEQ_FEEDER_STATS
#include "getticks.h"
#endif

class SeqFeeder
	{
public:
	SeqSource *m_SS = 0;
	vector<list<SeqInfo *> > m_SILists;
	vector<ObjMgr *> m_OMs;
	size_t m_MaxSIListSize = 100;
	vector<mymutex *> m_ListLocks;
	vector<mymutex *> m_OMLocks;
	uint m_ThreadCount = 0;
	bool m_EOF = false;
	uint m_FillWaitTime_ms = 100;
	uint m_GetWaitTime_ms = 10;
#if SEQ_FEEDER_STATS
	uint m_FillWaitCount = 0;
	uint m_GetWaitCount = 0;
	mutex m_StatsLock;
	TICKS m_FillWaitTicks = 0;
	TICKS m_GetWaitTicks = 0;
#endif

public:
	void Start(uint ThreadCount);

public:
	void ThreadBody();
	SeqInfo *GetSI(uint ThreadIndex);
	void Down(uint ThreadIndex, SeqInfo *SI);
	void Stats() const;

public:
	static void Static_ThreadBody(SeqFeeder *SF);
	};
