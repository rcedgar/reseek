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
	string m_Label_PDB;
	time_t m_TimeLastProgressMsg = 0;
	uint m_ChainCount = 0;
	bool m_Trace = false;

public:
	void Open(const string &FN);
	PDBChain *GetNext();

private:
	PDBChain *GetNextLo();
	PDBChain *PendingFile();
	void ReadNextDir();
	void Close();
	bool PushFileOrDir(const string &Path);
	PDBChain *GetFirst_CAL(const string &FN);
	PDBChain *GetNext_CAL();

	PDBChain *GetFirst_PDB(const string &FN);
	PDBChain *GetNext_PDB();

	PDBChain *GetFirst_CIF(const string &FN);
	PDBChain *GetNext_CIF();
	bool FileNameHasStructureExt(const string &FN) const;
	bool IsStructureExt(const string &Ext) const;

public:
	static void GetFallbackLabelFromFN(const string &FN,
	  string &Label);
	};
