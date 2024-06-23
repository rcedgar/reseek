#pragma once

#include "pdbchain.h"
#include "linereader2.h"
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
	string m_RootFileName;
	string m_CurrentFileName;
	STATE m_State = STATE_Closed;
	list<string> m_PendingFiles;
	list<string> m_PendingDirs;
	FILE *m_fCal = 0;
	LineReader2 m_LR;
	string m_Line;
	vector<string> m_Lines;
	vector<PDBChain *> m_Chains_PDB;
	vector<PDBChain *> m_Chains_CIF;
	uint m_ChainIdx_PDB = 0;
	uint m_ChainIdx_CIF = 0;

public:
	void Open(const string &FileName);
	PDBChain *GetNext();
	PDBChain *PendingFile();
	void ReadNextDir();
	void Close();
	void PushFileOrDir(const string &Path);
	PDBChain *GetFirst_CAL(const string &FileName);
	PDBChain *GetNext_CAL();

	PDBChain *GetFirst_PDB(const string &FileName);
	PDBChain *GetNext_PDB();

	PDBChain *GetFirst_CIF(const string &FileName);
	PDBChain *GetNext_CIF();

public:
	PDBChain *ChainFromLines_CAL(const vector<string> &Lines);
	PDBChain *ChainsFromLines_PDB(const vector<string> &Lines,
	  vector<PDBChain *> &Chains);
	PDBChain *ChainsFromLines_CIF(const vector<string> &Lines,
	  vector<PDBChain *> &Chains);
	};
