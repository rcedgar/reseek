#include "myutils.h"
#include "chainreader2.h"

void ChainReader2::Close()
	{
	m_State = STATE_Closed;
	}

void ChainReader2::Open(const string &Path)
	{
	m_RootFileName = Path;
	m_CurrentFileName.clear();
	m_PendingFiles.clear();
	m_PendingDirs.clear();
	PushFileOrDir(Path);
	m_State = STATE_PendingFile;
	}

void ChainReader2::PushFileOrDir(const string &Path)
	{
	if (Path.empty())
		Die("Empty pathname");
	if (EndsWith(Path, "/"))
		m_PendingDirs.push_back(Path);
	else
		m_PendingFiles.push_back(Path);
	}

void ChainReader2::ReadNextDir()
	{
	asserta(!m_PendingDirs.empty());
	const string &Dir = m_PendingDirs.front();
	m_PendingDirs.pop_front();
	vector<string> FileNames;
	vector<bool> IsSubDirs;
	mylistdir(Dir, FileNames, IsSubDirs);
	for (uint i = 0; i < SIZE(FileNames); ++i)
		{
		const string &Path = Dir + string("/") + FileNames[i];
		if (IsSubDirs[i])
			m_PendingDirs.push_back(Path);
		else
			m_PendingFiles.push_back(Path);
		}
	}

// Files first, then directories to reduce queue
PDBChain *ChainReader2::PendingFile()
	{
	for (uint SanityCounter = 0; SanityCounter < 100; ++SanityCounter)
		{
		if (m_PendingFiles.empty())
			ReadNextDir();
		if (m_PendingFiles.empty())
			return 0;

		m_CurrentFileName = m_PendingFiles.front();
		m_PendingFiles.pop_front();

		string Ext;
		GetExtFromPathName(m_CurrentFileName, Ext);
		if (Ext == "cal")
			{
			PDBChain *Chain = GetFirst_CAL(m_CurrentFileName);
			if (Chain != 0)
				return Chain;
			}
		else if (Ext == "pdb" || Ext == "ent")
			{
			PDBChain *Chain = GetFirst_PDB(m_CurrentFileName);
			if (Chain != 0)
				return Chain;
			}
		else if (Ext == "cif")
			{
			PDBChain *Chain = GetFirst_CIF(m_CurrentFileName);
			if (Chain != 0)
				return Chain;
			}
		else if (Ext == "files")
			{
			vector<string> Paths;
			ReadLinesFromFile(m_CurrentFileName, Paths);
			for (uint i = 0; i < SIZE(Paths); ++i)
				PushFileOrDir(Paths[i]);
			}
		}

	Die("Excessive looping in ChainReader2::PendingFile()");
	return 0;
	}

PDBChain *ChainReader2::GetNext()
	{
	for (uint SanityCounter = 0; SanityCounter < 100; ++SanityCounter)
		{
		switch (m_State)
			{
		case STATE_Closed:
			return 0;

		case STATE_PendingFile:
			{
			PDBChain *Chain = PendingFile();
			if (Chain == 0)
				Close();
			return Chain;
			}

		case STATE_ReadingCALFile:
			{
			PDBChain *Chain = GetNext_CAL();
			if (Chain != 0)
				return Chain;
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingPDBFile:
			{
			PDBChain *Chain = GetNext_PDB();
			if (Chain != 0)
				return Chain;
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingCIFFile:
			{
			PDBChain *Chain = GetNext_CIF();
			if (Chain != 0)
				return Chain;
			m_State = STATE_PendingFile;
			continue;
			}

		default:
			asserta(false);
			}
		}
	Die("Excessive looping in ChainReader2::GetNext()");
	return 0;
	}

void ChainReader2::GetFallbackLabelFromFN(const string &FileName, string &Label)
	{
	GetBaseName(FileName, Label);
	string Ext;
	GetExtFromPathName(FileName, Ext);
	ToLower(Ext);

// Special-case for downloaded PDB files e.g. pdb1iv1.ent
	if (Ext == "pdb" || Ext == "ent")
		{
		if (Label.size() == 7 && Label[0] == 'p' && Label[1] == 'd' && Label[2] == 'b')
			{
			Label = Label.substr(3, string::npos);
			ToUpper(Label);
			}
		}
	}

PDBChain *ChainFromCALLines(const vector<string> &Lines)
	{
	if (Lines.empty())
		return 0;
	PDBChain *Chain = new PDBChain;
	Chain->FromCalLines(Lines);
	return Chain;
	}

PDBChain *ChainReader2::GetFirst_CAL(const string &FileName)
	{
	m_LR.Open(FileName);
	bool Ok = m_LR.ReadLine(m_Line);
	if (!Ok)
		Die("Failed to read first line of CAL file '%s'",
		  FileName.c_str());
	return GetNext_CAL();
	}

PDBChain *ChainReader2::GetNext_CAL()
	{
	if (m_LR.m_EOF)
		{
		m_LR.Close();
		return 0;
		}
	if (m_Line.empty() || m_Line[0] != '>')
		Die("Expected '>' in CAL file '%s'",
		  m_CurrentFileName.c_str());
	m_Lines.clear();
	while (m_LR.ReadLine(m_Line))
		{
		if (m_Line.c_str()[0] == '>')
			break;
		m_Lines.push_back(m_Line);
		}
	return ChainFromCALLines(m_Lines);
	}

PDBChain *ChainReader2::GetFirst_PDB(const string &FileName)
	{
	ReadLinesFromFile(FileName, m_Lines);
	string FallbackLabel;
	GetFallbackLabelFromFN(FileName, FallbackLabel);
	PDBChain::ChainsFromLines_PDB(m_Lines, m_Chains_PDB, FallbackLabel);
	m_ChainIdx_PDB = 0;
	return GetNext_PDB();
	}

PDBChain *ChainReader2::GetFirst_CIF(const string &FileName)
	{
	ReadLinesFromFile(FileName, m_Lines);
	string FallbackLabel;
	GetFallbackLabelFromFN(FileName, FallbackLabel);
	PDBChain::ChainsFromLines_CIF(m_Lines, m_Chains_CIF, FallbackLabel);
	m_ChainIdx_CIF = 0;
	return GetNext_CIF();
	}

PDBChain *ChainReader2::GetNext_PDB()
	{
	const uint N = SIZE(m_Chains_PDB);
	if (m_ChainIdx_PDB == N)
		return 0;
	asserta(m_ChainIdx_PDB < N);
	return m_Chains_PDB[m_ChainIdx_PDB++];
	}

PDBChain *ChainReader2::GetNext_CIF()
	{
	const uint N = SIZE(m_Chains_CIF);
	if (m_ChainIdx_CIF == N)
		return 0;
	asserta(m_ChainIdx_CIF < N);
	return m_Chains_CIF[m_ChainIdx_CIF++];
	}
