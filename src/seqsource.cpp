#include "myutils.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"

#define TIME_LOCKS	0

#if TIME_LOCKS
#include "getticks.h"
static TICKS g_tLocks;
static TICKS g_tUnLocks;
#endif

SeqSource::SeqSource()
	{
	m_SI = m_OM.GetSeqInfo();
	m_SeqCount = 0;
	m_DoGetLock = true;
	}

SeqSource::~SeqSource()
	{
	m_OM.Down(m_SI);
	}

bool SeqSource::GetNext(SeqInfo *SI)
	{
	m_Lock.lock();
	bool Ok = GetNextLo(SI);
	m_Lock.unlock();

	if (!Ok)
		{
		SI->m_Label = 0;
		SI->m_Seq = 0;
		return false;
		}

	++m_SeqCount;
	return true;
	}
