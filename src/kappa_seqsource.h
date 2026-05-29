#pragma once

#include <stdio.h>
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

public:
	virtual bool GetIsNucleo() { return false; }

protected:
	virtual bool GetNextLo(SeqInfo *SI);

public:
	kappa_seqsource() {}
	virtual ~kappa_seqsource() {}

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
	};
