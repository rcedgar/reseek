#include "myutils.h"
#include "museqsource.h"
#include "seqinfo.h"
#include "alpha.h"

bool FastaFileIsNucleo(FILE *f);
char GetFeatureChar(byte Letter, uint AlphaSize);

MuSeqSource::MuSeqSource()
	{
	}

MuSeqSource::~MuSeqSource()
	{
	}

// Caller must own memory because SeqSource may be shared
// between threads, so SeqInfo should be thread-private.
bool MuSeqSource::GetNextLo(SeqInfo *SI)
	{
	if (m_IsFasta)
		return m_FSS.GetNext(SI);

	if (m_Chain != 0)
		delete m_Chain;
	m_Chain = m_CR.GetNext();
	if (m_Chain == 0)
		return false;
	m_DSS.Init(*m_Chain);
	const uint L = m_Chain->GetSeqLength();
	const uint AlphaSize = m_DSS.GetAlphaSize(FEATURE_Mu);
	asserta(AlphaSize == 36);
	SI->AllocL(L);
	SI->m_L = L;
	SI->SetLabel(m_Chain->m_Label.c_str());
	if (m_ASCII)
		{
		string Seq;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint Letter = m_DSS.GetFeature(FEATURE_Mu, Pos);
			char c = GetFeatureChar(Letter, AlphaSize);
			Seq += c;
			}
		const byte *ByteSeq = (const byte *) Seq.c_str();
		SI->SetCopy(m_SeqCount, m_Chain->m_Label.c_str(), ByteSeq, L);
		return true;
		}
	else
		{
		byte *Seq = SI->m_SeqBuffer;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint Letter = m_DSS.GetFeature(FEATURE_Mu, Pos);
			Seq[Pos] = Letter;
			}
		return true;
		}
	}

void MuSeqSource::OpenFasta(const string &FileName)
	{
	m_IsFasta = true;
	m_FSS.Open(FileName);
	}

void MuSeqSource::Open(const string &FileName, const DSSParams &Params)
	{
	m_IsFasta = false;
	m_Params = &Params;
	m_DSS.SetParams(Params);
	m_CR.Open(FileName);
	}

void MuSeqSource::Close()
	{
	}
