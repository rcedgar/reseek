#pragma once

#include "pdbchain.h"
#include <list>

class ChainReader2
	{
public:
	enum STATE
		{
		STATE_Closed,
		STATE_PendingFile,
		STATE_ReadingCALFile,
		STATE_ReadingPDBFile,
		STATE_ReadingCIFFile,
		};

public:
	mutex m_Lock;
	STATE m_State = STATE_Closed;
	list<string> m_PendingFiles;
	list<string> m_PendingDirs;
	FILE *m_fCal = 0;
	string m_CalFilePendingLine;

// One vector of strings per chain
	list<vector<string> > m_LinesList;

public:
	void Open(const string &FileName);
	PDBChain *GetNext();
	PDBChain *PendingFile();
	void ReadNextDir();
	void Close();
	void PushFileOrDir(const string &Path);
	PDBChain *StartReadingCALFile(const string &FileName) {return 0;}
	PDBChain *StartReadingPDBFile(const string &FileName) {return 0;}
	PDBChain *StartReadingCIFFile(const string &FileName) {return 0;}
	PDBChain *GetNext_PDBFile() {return 0;}
	PDBChain *GetNext_CIFFile() {return 0;}
	PDBChain *GetNext_CALFile() {return 0;}
	};
