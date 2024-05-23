#pragma once
#include <mutex>
#include "pdbchain.h"

class CalReader
	{
	mutex m_Lock;

public:
	FILE *m_f = 0;
	string m_PendingLine;
	vector<string> m_Lines;
	bool m_EOF = false;
	uint64 m_FileSize = UINT64_MAX;

public:
	void Clear()
		{
		m_f = 0;
		m_PendingLine.clear();
		m_Lines.clear();
		m_EOF = false;
		m_FileSize = UINT64_MAX;
		}

	void Close()
		{
		if (m_f != 0)
			CloseStdioFile(m_f);
		Clear();
		}

	void Open(const string &FileName);
	bool GetNext(PDBChain &Chain);
	void GetStrPctDone(string &s) const;
	double GetPctDone() const;
	};
