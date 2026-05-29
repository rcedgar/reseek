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

	case KSSS_seqdb:
		{
		uint idx = m_seqdbidx++;
		if (idx >= m_seqdb->GetSeqCount()) return false;

		const string &label = m_seqdb->GetLabel(idx);
		SI->SetLabel(label.c_str());

		uint L = m_seqdb->GetSeqLength(idx);
		SI->AllocL(L);
		const char *seq = m_seqdb->GetSeq(idx).c_str();
		for (uint i = 0; i < L; ++i)
			{
			uint8_t code = g_CharToLetterMu[seq[i]];
			if (code >= KAPPA_AS) code = 0;
			SI->m_SeqBuffer[i] = code;
			}
		SI->m_L = L;
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

void kappa_seqsource::OpenSeqDB(const SeqDB &DB)
	{
	m_KSSS = KSSS_chains;
	m_seqdb = &DB;
	m_seqdbidx = 0;
	}

void kappa_seqsource::Close()
	{
	}
