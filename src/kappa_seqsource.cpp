#include "myutils.h"
#include "kappa_seqsource.h"
#include "seqinfo.h"
#include "alpha.h"
#include "flat_params.h"

bool FastaFileIsNucleo(FILE *f);
char GetFeatureChar(byte Letter, uint AlphaSize);

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
		{
		uint idx = m_bcbidx++;
		SI->m_Index = idx;
		if (idx >= m_bcb->GetChainCount()) return false;
		const string &label = m_bcb->GetLabel(idx);
		SI->SetLabel(label.c_str());
		uint L = m_bcb->GetSeqLength(idx);
		if (L > flat_params::m_maxL)
			L = flat_params::m_maxL;
		SI->AllocL(L);
		uint L2 = m_bcb->read_codeseq_nu(SI->m_SeqBuffer, idx, L);
		asserta(L2 == L);
		SI->m_L = L;
		chaq::codeseq_nu_to_kappa_inplace(SI->m_SeqBuffer, L);
		return true;
		}

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
		uint pctx10 = uint((m_bcbidx*1000.0)/N);
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
	}
