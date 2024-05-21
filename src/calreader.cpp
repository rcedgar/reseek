#include "myutils.h"
#include "calreader.h"
#include "omplock.h"

void CalReader::Open(const string &FileName)
	{
	Clear();
	m_f = OpenStdioFile(FileName);
	m_FileSize = GetStdioFileSize64(m_f);
	bool Ok = ReadLineStdioFile(m_f, m_PendingLine);
	if (!Ok)
		Close();
	if (m_PendingLine[0] != '>')
		Die("Invalid .cal file, does not start with '>': %s",
		  FileName.c_str());
	}

bool CalReader::GetNext(PDBChain &Chain)
	{
	if (m_f == 0)
		return false;

	Lock("CalReader::GetNext");
	if (m_EOF)
		{
		Unlock("CalReader::GetNext");
		return false;
		}
	asserta(!m_PendingLine.empty() && m_Lines.empty());
	asserta(m_PendingLine[0] == '>');

	m_Lines.push_back(m_PendingLine);
	m_PendingLine.clear();
	string Line;
	for (;;)
		{
		bool Ok = ReadLineStdioFile(m_f, Line);
		if (!Ok)
			{
			m_EOF = true;
			break;
			}
		if (Line[0] == '>')
			{
			m_PendingLine = Line;
			break;
			}
		m_Lines.push_back(Line);
		}

	Chain.FromCalLines(m_Lines);
	m_Lines.clear();
	Unlock("CalReader::GetNext");
	return true;
	}

double CalReader::GetPctDone() const
	{
	if (m_FileSize == UINT64_MAX || m_FileSize == 0)
		return 0;
	uint64 Pos = GetStdioFilePos64(m_f);
	double Pct = double(Pos)*100.0/m_FileSize;
	return Pct;
	}

void CalReader::GetStrPctDone(string &s) const
	{
	s.clear();
	if (m_FileSize == UINT64_MAX || m_FileSize == 0)
		return;
	uint64 Pos = GetStdioFilePos64(m_f);
	double Pct = double(Pos)*100.0/m_FileSize;
	if (Pct < 0.1)
		Ps(s, "%.3f", Pct);
	else if (Pct < 1)
		Ps(s, "%.2f", Pct);
	else
		Ps(s, "%.1f", Pct);
	}
