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
	m_Lock.mylock();
	bool Ok = GetNextLo(SI);
	m_Lock.myunlock();

	if (!Ok)
		{
		SI->m_Label = 0;
		SI->m_Seq = 0;
		return false;
		}

	++m_SeqCount;
	return true;
	}

void SeqSource::LogLockStats(const char *name)
	{
	m_Lock.logme(name);
	}

void SeqSource::ResetLockStats()
	{
	m_Lock.reset_stats();
	}
