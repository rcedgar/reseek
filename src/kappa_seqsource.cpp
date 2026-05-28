#include "myutils.h"
#include "kappa_seqsource.h"
#include "seqinfo.h"
#include "alpha.h"

bool FastaFileIsNucleo(FILE *f);
char GetFeatureChar(byte Letter, uint AlphaSize);

// Caller must own memory because SeqSource may be shared
// between threads, so SeqInfo should be thread-private.
bool kappa_seqsource::GetNextLo(SeqInfo *SI)
	{
	if (m_IsFasta)
		{
		bool Ok = m_FSS.GetNext(SI);
		if (!Ok)
			return false;
		byte *Seq = SI->m_SeqBuffer;
		for (uint i = 0; i < SI->m_L; ++i)
			Seq[i] = g_CharToLetterMu[Seq[i]];
		return true;
		}

	if (m_chain != 0)
		delete m_chain;
	m_chain = m_CR.GetNext();
	Die("kappa_seqsource::GetNextLo, from flat_chain not implemented");
	return false;
	}

void kappa_seqsource::OpenFasta(const string &FileName)
	{
	m_IsFasta = true;
	m_FSS.Open(FileName);
	}

void kappa_seqsource::OpenChains(const string &FileName)
	{
	m_IsFasta = false;
	m_CR.Open(FileName);
	}

void kappa_seqsource::Close()
	{
	}
