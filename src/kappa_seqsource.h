#pragma once

#include <stdio.h>
#include <vector>
#include <deque>
#include <thread>
#include <mutex>
#include <condition_variable>
#include "fastaseqsource.h"
#include "flat_chain_reader.h"
#include "seqdb.h"
#include "bcadata.h"

enum KSS_SOURCE
	{
	KSSS_none,
	KSSS_fasta,
	KSSS_chains,
	KSSS_seqdb,
	KSSS_bcb
	};

struct KssBcbSlot
	{
	uint idx = UINT_MAX;
	uint L = 0;
	const string *label = 0;
	vector<uint8_t> kappa;
	};

struct KssBcbBatch
	{
	vector<KssBcbSlot> slots;
	uint count = 0;
	};

class kappa_seqsource : public SeqSource
	{
public:
	//bool m_IsFasta = false;
	flat_chain_reader m_CR;
	FASTASeqSource m_FSS;
	const flat_chain_t *m_chain = 0;
	const SeqDB *m_seqdb = 0;
	bool m_seqdb_codes = false;
	atomic<uint> m_seqdbidx = 0;
	atomic<uint> m_bcbidx = 0;
	const BCAData *m_bcb = 0;
	KSS_SOURCE m_KSSS = KSSS_none;

/////////////////////////////////////////////////////////////
// BCB async double-buffered prefetch.
//   A dedicated reader thread streams nu sequences from the
//   BCB file sequentially (one seek, then skip-and-read),
//   converts nu->kappa, and fills batches. kappa_filter workers claim
//   whole batches (claim_bcb_batch) and search lock-free; GetNext()
//   remains for other callers.
/////////////////////////////////////////////////////////////
	static const uint KSS_BCB_BATCH = 1024;
	static const uint KSS_BCB_NUM_BUFFERS = 16;

// Reader-thread-private scan cursor
	uint m_bcb_scan_next_idx = 0;

// Buffer pool and queues (guarded by m_bcb_qmutex)
	vector<KssBcbBatch *> m_bcb_all_batches;
	deque<KssBcbBatch *> m_bcb_free;
	deque<KssBcbBatch *> m_bcb_filled;
	mutex m_bcb_qmutex;
	condition_variable m_bcb_cv_filled;
	condition_variable m_bcb_cv_free;
	bool m_bcb_reader_eof = false;
	bool m_bcb_stop = false;
	bool m_bcb_reader_started = false;
	thread *m_bcb_reader = 0;

// Consumer-side current batch for GetNext() (guarded by SeqSource::m_Lock)
	KssBcbBatch *m_bcb_cur = 0;
	uint m_bcb_cur_pos = 0;

public:
	virtual bool GetIsNucleo() { return false; }

protected:
	virtual bool GetNextLo(SeqInfo *SI);
	bool get_next_bcb(SeqInfo *SI);
	void start_bcb_reader();
	void stop_bcb_reader();
	void bcb_reader_body();
	uint fill_bcb_batch(KssBcbBatch *batch);

public:
	kappa_seqsource() {}
	virtual ~kappa_seqsource() { stop_bcb_reader(); }

public:
	virtual unsigned GetPctDoneX10();

	virtual const char *GetFileNameC() const
		{ Die("kappa_seqsource::GetFileNameC()"); return 0; };
	virtual void Rewind()
		{ Die("kappa_seqsource::Rewind()"); };

public:
	void OpenFasta(const string &FileName);
	void OpenChains(const string &FileName);
	void OpenSeqDB(const SeqDB &DB, bool codes);
	void OpenBCB(const BCAData &bcb);
	void Close();

// Lock-free batch consumer for kappa_filter (uses m_bcb_qmutex only).
	KssBcbBatch *claim_bcb_batch();
	void release_bcb_batch(KssBcbBatch *batch);
	};
