#include "myutils.h"
#include "chainreader.h"

void ChainReader::Open(const string &FileName)
	{
	Clear();

	m_FileName = FileName;

	if (EndsWith(m_FileName, ".files"))
		{
		ReadLinesFromFile(FileName, m_FileNames);
		ReadChainsFromFileNameVec(m_FileNames, m_FilesChains);
		m_FilesChainIndex = 0;
		m_Type = CR_Files;
		}
	else if (EndsWith(m_FileName, ".cal"))
		{
		m_Type = CR_CAL;
		m_CR.Open(m_FileName);
		}
	else
		{
		m_Type = CR_PDB;
		ReadLinesFromFile(FileName, m_Lines);
		m_LineIndex = 0;
		}

	m_EndOfInput = false;
	}

bool ChainReader::GetNext(PDBChain &Chain)
	{
	if (m_EndOfInput)
		return false;

	bool Ok = false;
	{
	switch (m_Type)
		{
	case CR_None:
		Die("ChainReader::GetNext() not open");
	
	case CR_CAL:
		Ok = GetNext_CAL(Chain);
		break;
	
	case CR_PDB:
		Ok = GetNext_PDB(Chain);
		break;

	case CR_Files:
		Ok = GetNext_Files(Chain);
		break;
	
	default:
		asserta(false);
		}
	}

	return Ok;
	}

bool ChainReader::GetNext_CAL(PDBChain &Chain)
	{
	bool Ok = m_CR.GetNext(Chain);
	if (!Ok)
		{
		m_EndOfInput = true;
		return false;
		}
	return true;
	}

bool ChainReader::GetNext_PDB(PDBChain &Chain)
	{
	string Label;
	GetLabelFromFileName(m_FileName, Label);

	if (m_LineIndex >= SIZE(m_Lines))
		{
		m_EndOfInput = true;
		return false;
		}

	char CurrentChainChar = 0;
	vector<string> Lines;
	for (;;)
		{
		if (m_LineIndex >= SIZE(m_Lines))
			break;
		const string &Line = m_Lines[m_LineIndex++];
		if (Line == "TER   ")
			break;
		if (!PDBChain::IsATOMLine(Line))
			continue;
		char ChainChar = PDBChain::GetChainCharFromATOMLine(Line);
		if (CurrentChainChar == 0)
			CurrentChainChar = ChainChar;
		else
			{
			if (ChainChar != CurrentChainChar)
				{
				--m_LineIndex;
				if (CurrentChainChar != 0)
					Label += CurrentChainChar;
				char ChainChar2 = Chain.FromPDBLines(Label, Lines);
				if (ChainChar2 == 0)
					{
					CurrentChainChar = ChainChar;
					Lines.clear();
					Lines.push_back(Line);
					continue;
					}
				asserta(ChainChar2 == CurrentChainChar);
				return true;
				}
			}
		Lines.push_back(Line);
		}
	if (Lines.empty())
		{
		m_EndOfInput = true;
		return false;
		}

	PDBChain::AppendChainToLabel(Label, CurrentChainChar);
	Chain.FromPDBLines(Label, Lines);
	return true;
	}

bool ChainReader::GetNext_Files(PDBChain &Chain)
	{
	if (m_FilesChainIndex >= SIZE(m_FilesChains))
		{
		m_EndOfInput = true;
		return false;
		}
	Chain = *m_FilesChains[m_FilesChainIndex++];
	return true;
	}

uint ChainReader::GetMilDone()
	{
	if (m_EndOfInput)
		return 1000;
	if (m_PctDone == 0)
		{
		m_PctDone = 0.0001;
		return 0;
		}
	m_PctDone = GetPctDone();
	uint Mil = uint(m_PctDone/1001);
	if (Mil == 0)
		Mil = 1;
	if (Mil >= 1001)
		Mil = 1000;
	return Mil;
	}

double ChainReader::GetPctDone() const
	{
	switch (m_Type)
		{
	case CR_None:	return 0;
	case CR_PDB:	return GetPctDone_PDB();
	case CR_CAL:	return GetPctDone_CAL();
	case CR_Files:	return GetPctDone_Files();
		}
	asserta(false);
	return 0;
	}

double ChainReader::GetPctDone_CAL() const
	{
	return m_CR.GetPctDone();
	}

double ChainReader::GetPctDone_PDB() const
	{
	return 50.0;
	}

double ChainReader::GetPctDone_Files() const
	{
	double Pct = GetPct(m_FilesChainIndex, SIZE(m_FilesChains) + 1);
	return Pct;
	}

void ChainReader::GetStrPctDone(string &s) const
	{
	s.clear();
	double Pct = GetPctDone();
	if (Pct < 0.1)
		Ps(s, "%.3f", Pct);
	else if (Pct < 1)
		Ps(s, "%.2f", Pct);
	else
		Ps(s, "%.1f", Pct);
	}
