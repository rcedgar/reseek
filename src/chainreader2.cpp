#include "myutils.h"
#include "chainreader2.h"

void ChainReader2::Close()
	{
	m_State = STATE_Closed;
	}

void ChainReader2::Open(const string &Path)
	{
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

		const string FileName = m_PendingFiles.front();
		m_PendingFiles.pop_front();

		string Ext;
		GetExtFromPathName(FileName, Ext);
		if (Ext == "cal")
			{
			PDBChain *Chain = StartReadingCALFile(FileName);
			if (Chain != 0)
				return Chain;
			}
		else if (Ext == "pdb" || Ext == "ent")
			{
			PDBChain *Chain = StartReadingPDBFile(FileName);
			if (Chain != 0)
				return Chain;
			}
		else if (Ext == "cif")
			{
			PDBChain *Chain = StartReadingCIFFile(FileName);
			if (Chain != 0)
				return Chain;
			}
		else if (Ext == "files")
			{
			vector<string> Paths;
			ReadLinesFromFile(FileName, Paths);
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
			PDBChain *Chain = GetNext_CALFile();
			if (Chain != 0)
				return Chain;
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingPDBFile:
			{
			PDBChain *Chain = GetNext_PDBFile();
			if (Chain != 0)
				return Chain;
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingCIFFile:
			{
			PDBChain *Chain = GetNext_CIFFile();
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
