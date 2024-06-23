#include "myutils.h"
#include "chainreader2.h"

void ChainReader2::Close()
	{
	if (m_Trace) Log("ChainReader2::Close()\n");
	m_State = STATE_Closed;
	}

void ChainReader2::Open(const string &Path)
	{
	m_Trace = opt_trace_chainreader2;
	if (m_Trace) Log("ChainReader2::Open(%s)\n", Path.c_str());
	m_RootFileName = Path;
	m_CurrentFileName.clear();
	m_PendingFiles.clear();
	m_PendingDirs.clear();
	PushFileOrDir(Path);
	m_State = STATE_PendingFile;
	m_TimeLastProgressMsg = time(0);
	m_ChainCount = 0;
	}

void ChainReader2::PushFileOrDir(const string &Path)
	{
	if (m_Trace) Log("ChainReader2::PushFileOrDir(%s)\n", Path.c_str());
	if (Path.empty())
		Die("Empty pathname");
	if (EndsWith(Path, "/"))
		m_PendingDirs.push_back(Path);
	else
		m_PendingFiles.push_back(Path);
	}

void ChainReader2::ReadNextDir()
	{
	if (m_PendingDirs.empty())
		return;
	const string Dir = m_PendingDirs.front();
	if (m_Trace) Log("ChainReader2::ReadNextDir(%s)\n", Dir.c_str());
	m_PendingDirs.pop_front();
	vector<string> FileNames;
	vector<bool> IsSubDirs;
	mylistdir(Dir, FileNames, IsSubDirs);
	for (uint i = 0; i < SIZE(FileNames); ++i)
		{
		const string &FN = FileNames[i];
		if (FN[0] == '.')
			continue;
		const string &Path = Dir + string("/") + FN;
		if (IsSubDirs[i])
			{
			if (m_Trace) Log("  m_PendingDirs.push_back(%s)\n", Path.c_str());
			m_PendingDirs.push_back(Path);
			}
		else
			{
			if (m_Trace) Log("  m_PendingFiles.push_back(%s)\n", Path.c_str());
			m_PendingFiles.push_back(Path);
			}
		}
	}

// Files first, then directories to reduce queue
PDBChain *ChainReader2::PendingFile()
	{
	for (uint SanityCounter = 0; SanityCounter < 100; ++SanityCounter)
		{
		if (m_Trace) Log("ChainReader2::PendingFile() SanityCounter=%u\n", SanityCounter);
		if (m_PendingFiles.empty())
			ReadNextDir();
		if (m_PendingFiles.empty())
			return 0;

		m_CurrentFileName = m_PendingFiles.front();
		m_PendingFiles.pop_front();

		string Ext;
		GetExtFromPathName(m_CurrentFileName, Ext);
		ToLower(Ext);
		if (m_Trace) Log("  m_CurrentFileName=%s Ext=%s\n", m_CurrentFileName.c_str(), Ext.c_str());
		if (Ext == "cal")
			{
			m_State = STATE_ReadingCALFile;
			PDBChain *Chain = GetFirst_CAL(m_CurrentFileName);
			if (Chain != 0)
				return Chain;
			}
		else if (Ext == "pdb" || Ext == "ent")
			{
			m_State = STATE_ReadingPDBFile;
			PDBChain *Chain = GetFirst_PDB(m_CurrentFileName);
			if (Chain != 0)
				return Chain;
			}
		else if (Ext == "cif")
			{
			m_State = STATE_ReadingCIFFile;
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
	for (;;)
		{
		time_t Now = time(0);
		PDBChain *Chain = GetNextLo();
		if (Chain == 0)
			{ if (m_Trace) Log("ChainReader2::GetNextLo() Chain=0\n"); }
		else
			{ if (m_Trace) Log("ChainReader2::GetNextLo() Chain=%s\n", Chain->m_Label.c_str()); }

		if (Chain == 0)
			{
			Progress("%u chains read               \n", m_ChainCount);
			return 0;
			}
		brk(Chain->m_Label.find(".pdb") != string::npos);
		++m_ChainCount;
		if (Now > m_TimeLastProgressMsg)
			{
			m_TimeLastProgressMsg = Now;
			Progress("%u chains read >%s\r",
			  m_ChainCount, Chain->m_Label.c_str());
			}
		if (Chain->GetSeqLength() > 0)
			return Chain;
		delete Chain;
		}
	return 0;
	}

PDBChain *ChainReader2::GetNextLo()
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
		Die("%s: Expected '>' in CAL file",
		  m_CurrentFileName.c_str());

	const string Label = m_Line.substr(1);
	if (m_Trace) Log("ChainReader2::GetNext_CAL() Label=%s\n", Label.c_str());
	m_Lines.clear();
	while (m_LR.ReadLine(m_Line))
		{
		if (m_Line.c_str()[0] == '>')
			break;
		m_Lines.push_back(m_Line);
		}

	PDBChain *Chain = new PDBChain;
	Chain->m_Label = Label;

/***
>102l
M       43.619  -1.924  8.869
N       40.445  -0.876  10.670
I       38.254  2.240   11.220
F       40.340  3.621   14.036
***/
	const uint N = SIZE(m_Lines);
	vector<string> Fields;
	for (uint LineNr = 1; LineNr < N; ++LineNr)
		{
		const string &Line = m_Lines[LineNr];
		Split(Line, Fields, '\t');
		if (Fields.size() != 4 || Fields[0].size() != 1)
			Die("%s: Invalid CAL record '%s'",
			  m_CurrentFileName.c_str(), Line.c_str());

		char aa = Fields[0][0];
		double X = StrToFloat(Fields[1]);
		double Y = StrToFloat(Fields[2]);
		double Z = StrToFloat(Fields[3]);

		Chain->m_Seq.push_back(aa);
		Chain->m_Xs.push_back(X);
		Chain->m_Ys.push_back(Y);
		Chain->m_Zs.push_back(Z);
		}
	return Chain;
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
	PDBChain *Chain = m_Chains_PDB[m_ChainIdx_PDB++];
	if (m_Trace) Log("ChainReader2::GetNext_PDB() %u/%u Label=%s\n", m_ChainIdx_PDB, N, Chain->m_Label.c_str());
	return Chain;
	}

PDBChain *ChainReader2::GetNext_CIF()
	{
	const uint N = SIZE(m_Chains_CIF);
	if (m_ChainIdx_CIF == N)
		return 0;
	asserta(m_ChainIdx_CIF < N);
	PDBChain *Chain = m_Chains_CIF[m_ChainIdx_CIF++];
	if (m_Trace) Log("ChainReader2::GetNext_CIF() %u/%u Label=%s\n", m_ChainIdx_CIF, N, Chain->m_Label.c_str());
	return Chain;
	}
