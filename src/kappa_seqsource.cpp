#include "myutils.h"
#include "kappa_seqsource.h"
#include "seqinfo.h"
#include "alpha.h"
#include "flat_params.h"
#include "chaq.h"

bool FastaFileIsNucleo(FILE *f);
char GetFeatureChar(byte Letter, uint AlphaSize);

// Fills one batch by sequentially reading nu sequences and
// converting to kappa. Runs ONLY on the reader thread, which is
// the sole user of the BCB FILE* during the scan.
uint kappa_seqsource::fill_bcb_batch(KssBcbBatch *batch)
	{
	const uint ChainCount = m_bcb->GetChainCount();
	FILE *f = m_bcb->m_f;
	const uint maxL = flat_params::m_maxL;

	if (batch->slots.size() < KSS_BCB_BATCH)
		batch->slots.resize(KSS_BCB_BATCH);

	uint n = 0;
	for (; n < KSS_BCB_BATCH; ++n)
		{
		if (m_bcb_scan_next_idx >= ChainCount)
			break;

		const uint idx = m_bcb_scan_next_idx++;
		const uint origL = m_bcb->GetSeqLength(idx);
		uint L = origL;
		if (L > maxL)
			L = maxL;

#if DEBUG
		asserta(GetStdioFilePos64(f) == m_bcb->get_offset_nuseq(idx));
#endif

		KssBcbSlot &slot = batch->slots[n];
		slot.idx = idx;
		slot.L = L;
		slot.label = &m_bcb->GetLabel(idx);
		slot.kappa.resize(L);

		if (L > 0)
			{
			const uint64 nL = (uint64) fread(slot.kappa.data(), 1, L, f);
			if (nL != L)
				Die("kappa_seqsource::fill_bcb_batch() idx=%u L=%u", idx, L);
			chaq::codeseq_nu_to_kappa_inplace(slot.kappa.data(), L);
			}
		if (origL > L)
			SetStdioFilePos64(f, GetStdioFilePos64(f) + uint64(origL - L));
		}
	batch->count = n;
	return n;
	}

void kappa_seqsource::bcb_reader_body()
	{
	asserta(m_bcb != 0);
	FILE *f = m_bcb->m_f;
	asserta(f != 0);
	const uint ChainCount = m_bcb->GetChainCount();

	m_bcb_scan_next_idx = 0;
	if (ChainCount > 0)
		SetStdioFilePos64(f, m_bcb->get_offset_nuseq(0));

	for (;;)
		{
		KssBcbBatch *batch = 0;
			{
			std::unique_lock<mutex> lk(m_bcb_qmutex);
			m_bcb_cv_free.wait(lk,
				[this]{ return !m_bcb_free.empty() || m_bcb_stop; });
			if (m_bcb_stop)
				return;
			batch = m_bcb_free.front();
			m_bcb_free.pop_front();
			}

		const uint n = fill_bcb_batch(batch);
		const bool done = (m_bcb_scan_next_idx >= ChainCount);

			{
			std::lock_guard<mutex> lk(m_bcb_qmutex);
			if (n > 0)
				m_bcb_filled.push_back(batch);
			else
				m_bcb_free.push_back(batch);
			if (done)
				m_bcb_reader_eof = true;
			}
		m_bcb_cv_filled.notify_all();

		if (done)
			return;
		}
	}

void kappa_seqsource::start_bcb_reader()
	{
	if (m_bcb_reader_started)
		return;
	m_bcb_reader_started = true;
	m_bcb_reader_eof = false;
	m_bcb_stop = false;
	m_bcb_scan_next_idx = 0;
	m_bcb_cur = 0;
	m_bcb_cur_pos = 0;

	m_bcb_all_batches.clear();
	m_bcb_free.clear();
	m_bcb_filled.clear();
	uint nthreads = GetRequestedThreadCount();
	if (nthreads == 0)
		nthreads = 1;
	const uint num_buffers = KSS_BCB_BUFFERS_PER_THREAD * nthreads;
	for (uint i = 0; i < num_buffers; ++i)
		{
		KssBcbBatch *b = new KssBcbBatch;
		m_bcb_all_batches.push_back(b);
		m_bcb_free.push_back(b);
		}

	m_bcb_reader = new thread(&kappa_seqsource::bcb_reader_body, this);
	}

void kappa_seqsource::stop_bcb_reader()
	{
	if (!m_bcb_reader_started)
		return;
		{
		std::lock_guard<mutex> lk(m_bcb_qmutex);
		m_bcb_stop = true;
		}
	m_bcb_cv_free.notify_all();
	if (m_bcb_reader != 0)
		{
		if (m_bcb_reader->joinable())
			m_bcb_reader->join();
		delete m_bcb_reader;
		m_bcb_reader = 0;
		}
	for (size_t i = 0; i < m_bcb_all_batches.size(); ++i)
		delete m_bcb_all_batches[i];
	m_bcb_all_batches.clear();
	m_bcb_free.clear();
	m_bcb_filled.clear();
	m_bcb_cur = 0;
	m_bcb_cur_pos = 0;
	m_bcb_reader_eof = false;
	m_bcb_stop = false;
	m_bcb_reader_started = false;
	}

KssBcbBatch *kappa_seqsource::claim_bcb_batch()
	{
	if (!m_bcb_reader_started)
		start_bcb_reader();

	std::unique_lock<mutex> lk(m_bcb_qmutex);
	m_bcb_cv_filled.wait(lk,
		[this]{ return !m_bcb_filled.empty() || m_bcb_reader_eof; });
	if (m_bcb_filled.empty())
		return 0;
	KssBcbBatch *batch = m_bcb_filled.front();
	m_bcb_filled.pop_front();
	return batch;
	}

void kappa_seqsource::release_bcb_batch(KssBcbBatch *batch)
	{
	asserta(batch != 0);
		{
		std::lock_guard<mutex> lk(m_bcb_qmutex);
		m_bcb_free.push_back(batch);
		}
	m_bcb_cv_free.notify_one();
	}

bool kappa_seqsource::get_next_bcb(SeqInfo *SI)
	{
	if (!m_bcb_reader_started)
		start_bcb_reader();

	for (;;)
		{
		if (m_bcb_cur != 0 && m_bcb_cur_pos < m_bcb_cur->count)
			{
			const KssBcbSlot &slot = m_bcb_cur->slots[m_bcb_cur_pos++];
			SI->m_Index = slot.idx;
			SI->SetLabel(slot.label->c_str());
			SI->AllocL(slot.L);
			if (slot.L > 0)
				memcpy(SI->m_SeqBuffer, slot.kappa.data(), slot.L);
			SI->m_L = slot.L;
			++m_bcb_done_count;
			return true;
			}

		if (m_bcb_cur != 0)
			{
			release_bcb_batch(m_bcb_cur);
			m_bcb_cur = 0;
			m_bcb_cur_pos = 0;
			}

		std::unique_lock<mutex> lk(m_bcb_qmutex);
		m_bcb_cv_filled.wait(lk,
			[this]{ return !m_bcb_filled.empty() || m_bcb_reader_eof; });
		if (m_bcb_filled.empty())
			return false;
		m_bcb_cur = m_bcb_filled.front();
		m_bcb_filled.pop_front();
		m_bcb_cur_pos = 0;
		}
	}

// Caller must own memory because SeqSource may be shared
// between threads, so SeqInfo should be thread-private.
bool kappa_seqsource::GetNextLo(SeqInfo *SI)
	{
	switch (m_KSSS)
		{
	case KSSS_fasta:
		{
		bool Ok = m_FSS.GetNext(SI);
		if (!Ok)
			return false;
		byte *Seq = SI->m_SeqBuffer;
		for (uint i = 0; i < SI->m_L; ++i)
			Seq[i] = g_CharToLetterMu[Seq[i]];
		return true;
		}

	case KSSS_chains:
		{
		if (m_chain != 0)
			delete m_chain;
		m_chain = m_CR.GetNext();
		Die("kappa_seqsource::GetNextLo, from flat_chain not implemented");
		return false;
		}

	case KSSS_bcb:
		return get_next_bcb(SI);

	case KSSS_seqdb:
		{
		uint idx = m_seqdbidx++;
		SI->m_Index = idx;
		if (idx >= m_seqdb->GetSeqCount()) return false;

		const string &label = m_seqdb->GetLabel(idx);
		SI->SetLabel(label.c_str());

		uint L = m_seqdb->GetSeqLength(idx);
		const char *seq = m_seqdb->GetSeq(idx).c_str();
		SI->m_L = L;
		if (m_seqdb_codes)
			{
			SI->m_Seq = (byte *) seq;
#if DEBUG
			for (uint i = 0; i < L; ++i) 
				assert(seq[i] < KAPPA_AS);
#endif
			}
		else
			{
			SI->AllocL(L);
			for (uint i = 0; i < L; ++i)
				{
				uint8_t code = g_CharToLetterMu[seq[i]];
				if (code >= KAPPA_AS) code = 0;
				SI->m_SeqBuffer[i] = code;
				}
			SI->m_Seq = SI->m_SeqBuffer;
			}
		return true;
		}

	default: Die("m_KSSS=%d", int(m_KSSS));
		}

	return false;
	}

void kappa_seqsource::OpenFasta(const string &FileName)
	{
	m_KSSS = KSSS_fasta;
	m_seqdb = 0;
	m_FSS.Open(FileName);
	}

void kappa_seqsource::OpenChains(const string &FileName)
	{
	m_KSSS = KSSS_chains;
	m_seqdb = 0;
	m_CR.Open(FileName);
	}

void kappa_seqsource::OpenBCB(const BCAData &bcb)
	{
	m_KSSS = KSSS_bcb;
	m_bcb = &bcb;
	m_bcb_done_count = 0;
	}

void kappa_seqsource::OpenSeqDB(const SeqDB &DB, bool codes)
	{
	m_KSSS = KSSS_seqdb;
	m_seqdb = &DB;
	m_seqdbidx = 0;
	m_seqdb_codes = codes;
	}

uint kappa_seqsource::GetPctDoneX10()
	{
	switch (m_KSSS)
		{
	case KSSS_fasta:
		return m_FSS.GetPctDoneX10();

	case KSSS_chains:
		break;

	case KSSS_seqdb:
		{
		uint N = m_seqdb->GetSeqCount();
		uint pctx10 = uint((m_seqdbidx*1000.0)/N);
		if (pctx10 == 0) pctx10 = 1;
		if (pctx10 >= 999) pctx10 = 998;
		return pctx10;
		}

	case KSSS_bcb:
		{
		uint N = m_bcb->GetChainCount();
		uint pctx10 = uint((m_bcb_done_count*1000.0)/N);
		if (pctx10 == 0) pctx10 = 1;
		if (pctx10 >= 999) pctx10 = 998;
		return pctx10;
		}

	default: asserta(false);
		}

	Die("kappa_seqsource::GetPctDoneX10() not implemented");
	return 0;
	}

void kappa_seqsource::Close()
	{
	stop_bcb_reader();
	}
