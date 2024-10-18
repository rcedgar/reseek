#pragma once

#include "pdbchain.h"
#include "linereader2.h"
#include "pdbfilescanner.h"
#include "bcadata.h"
#include <list>

class ChainReader2
	{
public:
	enum STATE
		{
		STATE_Closed,
		STATE_PendingFile,
		STATE_ReadingCALFile,
		STATE_ReadingBCAFile,
		STATE_ReadingPDBFile,
		STATE_ReadingCIFFile,
		STATE_ReadingVec,
		};

public:
	mutex m_Lock;
	STATE m_State = STATE_Closed;
	FILE *m_fCal = 0;
	LineReader2 m_LR;
	string m_Line;
	vector<string> m_Lines;
	vector<PDBChain *> m_Chains_PDB;
	vector<PDBChain *> m_Chains_CIF;
	vector<PDBChain *> *m_ptrChains = 0;
	uint m_ChainIdx_PDB = 0;
	uint m_ChainIdx_CIF = 0;
	uint m_ChainIdx_Vec = 0;
	string m_Label_PDB;
	BCAData m_BCA;
	uint64 m_ChainIdx_BCA = 0;
	uint m_ChainCount = 0;
	string m_CurrentFN;
	bool m_Trace = false;
	bool m_SaveLines = false;

// FS object shared with other threads
	PDBFileScanner *m_ptrFS = 0;

public:
	void Open(const string &FileName);
	void Open(PDBFileScanner &FS);
	void Open(vector<PDBChain *> &Chains);
	PDBChain *GetNext();

private:
	void Close();
	PDBChain *GetNextLo1();
	PDBChain *GetFirst(const string &FN);

	PDBChain *GetFirst_BCA(const string &FN);
	PDBChain *GetNext_BCA();

	PDBChain *GetFirst_CAL(const string &FN);
	PDBChain *GetNext_CAL();

	PDBChain *GetFirst_PDB(const string &FN);
	PDBChain *GetNext_PDB();

	PDBChain *GetFirst_CIF(const string &FN);
	PDBChain *GetNext_CIF();

	PDBChain *GetNext_Vec();

	PDBChain *ChainFromLines_CAL(const vector<string> &Lines) const;
	void ChainsFromLines_PDB(const vector<string> &Lines,
		vector<PDBChain *> &Chains, const string &Label) const;
	void ChainsFromLines_CIF(const vector<string> &Lines,
		vector<PDBChain *> &Chains, const string &FallbackLabel) const;
	bool IsATOMLine_PDB(const string &Line) const;
	bool IsChainEndLine_PDB(const string &Line) const;

public:
	static void GetFallbackLabelFromFN(const string &FN, string &Label);
	};
