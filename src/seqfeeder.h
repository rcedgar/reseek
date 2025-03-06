#pragma once

#include "mymutex.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"

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
	uint m_FillWaitTime_ms = 1000;
	uint m_GetWaitTime_ms = 10;
private:
	SeqFeeder();

public:
	SeqFeeder(uint ThreadCount, SeqSource &SS);

public:
	void ThreadBody();
	SeqInfo *GetSI(uint ThreadIndex);
	void Down(uint ThreadIndex, SeqInfo *SI);

public:
	static void Static_ThreadBody(SeqFeeder *SF);
	};
