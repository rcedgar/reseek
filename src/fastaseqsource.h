#ifndef fastaseqsource_h
#define fastaseqsource_h

#include <stdio.h>
#include "fileseqsource.h"

class FASTASeqSource : public FileSeqSource
	{
public:
	bool m_AllowDigits = false;

public:
	virtual bool GetIsNucleo();

protected:
	virtual bool GetNextLo(SeqInfo *SI);

public:
	FASTASeqSource();
	virtual ~FASTASeqSource();
	};

#endif // fastaseqsource_h
