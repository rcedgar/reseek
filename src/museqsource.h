#pragma once

#include <stdio.h>
#include "seqsource.h"
#include "chainreader2.h"
#include "dss.h"

class MuSeqSource : public SeqSource
	{
public:
	ChainReader2 m_CR;
	const DSSParams *m_Params = 0;
	const PDBChain *m_Chain = 0;
	DSS m_DSS;

public:
	virtual bool GetIsNucleo() { return false; }

protected:
	virtual bool GetNextLo(SeqInfo *SI);

public:
	MuSeqSource();
	virtual ~MuSeqSource();

public:
	virtual unsigned GetPctDoneX10() { Die("MuSeqSource::GetPctDoneX10()"); return 0; };
	virtual const char *GetFileNameC() const { return "MuSeqSource::GetFileNameC()"; return 0; };
	virtual void Rewind() { Die("MuSeqSource::Rewind()"); };

public:
	void Open(const string &FileName, const DSSParams &Params);
	void Close();
	};
