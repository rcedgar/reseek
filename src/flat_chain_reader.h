#pragma once

#include "flat_chain.h"
#include "linereader2.h"
#include "pdbfilescanner.h"
#include "bcadata.h"
#include <list>

class flat_chain_reader
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
		};

public:
	mutex m_CRPerThreadLock;
	mutex m_CRGlobalLock;
	STATE m_State = STATE_Closed;
	FILE *m_fCal = 0;
	LineReader2 m_LR;
	string m_Line;
	vector<string> m_Lines;
	vector<flat_chain *> m_Chains_PDB;
	vector<flat_chain *> m_Chains_CIF;
	vector<flat_chain *> *m_ptrChains = 0;
	uint m_ChainIdx_PDB = 0;
	uint m_ChainIdx_CIF = 0;
	uint m_ChainIdx_Vec = 0;
	string m_Label_PDB;
	BCAData m_BCA;
	uint64 m_ChainIdx_BCA = 0;
	string m_CurrentFN;
	bool m_Trace = false;
	bool m_SaveLines = false;

// FS object shared with other threads
	PDBFileScanner *m_ptrFS = 0;

public:
	static uint m_CRGlobalChainCount;
	static uint m_CRGlobalFormatErrors;

public:
	void Open(const string &FileName);
	void Open(PDBFileScanner &FS);
	void Open(vector<flat_chain *> &Chains);
	flat_chain *GetNext();

private:
	void Close();
	flat_chain *GetNextLo1();
	flat_chain *GetFirst(const string &FN);

	flat_chain *GetFirst_BCA(const string &FN);
	flat_chain *GetNext_BCA();

	flat_chain *GetFirst_CAL(const string &FN);
	flat_chain *GetNext_CAL();

	flat_chain *GetFirst_PDB(const string &FN);
	flat_chain *GetNext_PDB();

	flat_chain *GetFirst_CIF(const string &FN);
	flat_chain *GetNext_CIF();

	flat_chain *GetNext_Vec();

	void ChainsFromLines_PDB(const vector<string> &Lines,
		vector<flat_chain *> &Chains, const string &Label) const;
	void ChainsFromLines_CIF(const vector<string> &Lines,
		vector<flat_chain *> &Chains, const string &FallbackLabel);
	bool IsATOMLine_PDB(const string &Line) const;
	bool IsChainEndLine_PDB(const string &Line) const;
	uint GetCIFFieldIdx(const map<string, uint> &FieldToIdx, const string &Name);
	void IncFormatErrors();

public:
	static void GetFallbackLabelFromFN(const string &FN, string &Label);
	};

void read_flat_chains(const string &fn, vector<flat_chain *> &chains);