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
	bool Ok = PushFileOrDir(Path);
	if (!Ok)
		Die("Not found '%s'", Path.c_str());
	m_State = STATE_PendingFile;
	m_TimeLastProgressMsg = time(0);
	m_ChainCount = 0;
	}

bool ChainReader2::IsStructureExt(const string &Ext) const
	{
	string LowerExt = Ext;
	ToLower(LowerExt);

#define x(s)	if (Ext == #s || Ext == string(#s) + string(".gz")) return true;
	x(pdb)
	x(ent)
	x(cif)
	x(mmcif)
	x(cal)
#undef x
	return false;
	}

bool ChainReader2::FileNameHasStructureExt(const string &FN) const
	{
	string Ext;
	GetExtFromPathName(FN, Ext);
	return IsStructureExt(Ext);
	}

bool ChainReader2::PushFileOrDir(const string &Path)
	{
	if (m_Trace) Log("ChainReader2::PushFileOrDir(%s)\n", Path.c_str());
	if (Path.empty())
		return false;
	if (IsDirectory(Path))
		{
		m_PendingDirs.push_back(Path);
		return true;
		}
	else if (IsRegularFile(Path))
		{
		m_PendingFiles.push_back(Path);
		return true;
		}
	return false;
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
		if (FN == "." || FN == "..")
			continue; 
		const string &Path = Dir + string("/") + FN;
		if (IsDirectory(Path))
			m_PendingDirs.push_back(Path);
		else if (FileNameHasStructureExt(FN))
			m_PendingFiles.push_back(Path);
		}
	}

// Files first, then directories to reduce queue
PDBChain *ChainReader2::PendingFile()
	{
	for (uint SanityCounter = 0; SanityCounter < 100; ++SanityCounter)
		{
		if (m_Trace) Log("ChainReader2::PendingFile() SanityCounter=%u\n", SanityCounter);
		while (m_PendingFiles.empty() && !m_PendingDirs.empty())
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
		else if (Ext == "pdb" || Ext == "pdb.gz" || Ext == "ent" || Ext == "ent.gz")
			{
			m_State = STATE_ReadingPDBFile;
			PDBChain *Chain = GetFirst_PDB(m_CurrentFileName);
			if (Chain != 0)
				return Chain;
			}
		else if (Ext == "cif" || Ext == "cif.gz" || Ext == "mmcif" || Ext == "mmcif.gz")
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
	m_Lock.lock();
	PDBChain *Chain = GetNextLo1();
	m_Lock.unlock();
	return Chain;
	}

PDBChain *ChainReader2::GetNextLo1()
	{
	uint SanityCounter = 0;
	for (;;)
		{
		if (m_Trace) Log("GetNextLo1() state=%d\n", m_State);
		time_t Now = time(0);
		PDBChain *Chain = GetNextLo2();
		if (Chain == 0)
			{
			++SanityCounter;
			if (SanityCounter >= 100)
				Die("Excessive looping in ChainReader2::GetNextLo1()");
			}
		if (Chain == 0)
			{ if (m_Trace) Log("ChainReader2::GetNextLo1() Chain=0\n"); }
		else
			{ if (m_Trace) Log("ChainReader2::GetNextLo1() Chain=%s\n", Chain->m_Label.c_str()); }

		if (Chain == 0)
			{
			Progress("%s chains                \n", IntToStr(m_ChainCount));
			return 0;
			}
		++m_ChainCount;
		if (Now > m_TimeLastProgressMsg)
			{
			m_TimeLastProgressMsg = Now;
			Progress("%s chains >%s\r",
			  IntToStr(m_ChainCount), Chain->m_Label.c_str());
			}
		if (Chain->GetSeqLength() > 0)
			return Chain;
		delete Chain;
		}
	return 0;
	}

PDBChain *ChainReader2::GetNextLo2()
	{
	for (uint SanityCounter = 0; SanityCounter < 100; ++SanityCounter)
		{
		if (m_Trace) Log("GetNextLo2() state=%d\n", m_State);
		switch (m_State)
			{
		case STATE_Closed:
			return 0;

		case STATE_PendingFile:
			{
			PDBChain *Chain = PendingFile();
			if (m_Trace) Log("PendingFile()=0, Close() and stop.\n");
			if (Chain == 0)
				Close();
			return Chain;
			}

		case STATE_ReadingCALFile:
			{
			PDBChain *Chain = GetNext_CAL();
			if (Chain != 0)
				return Chain;
			if (m_Trace) Log("GetNext_CAL()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingPDBFile:
			{
			PDBChain *Chain = GetNext_PDB();
			if (Chain != 0)
				return Chain;
			if (m_Trace) Log("GetNext_PDB()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingCIFFile:
			{
			PDBChain *Chain = GetNext_CIF();
			if (Chain != 0)
				return Chain;
			if (m_Trace) Log("GetNext_CIF()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		default:
			asserta(false);
			}
		}
	Die("Excessive looping in ChainReader2::GetNextLo2()");
	return 0;
	}

void ChainReader2::GetFallbackLabelFromFN(const string &FN, string &Label)
	{
	GetStemName(FN, Label);
	string Ext;
	GetExtFromPathName(FN, Ext);
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

PDBChain *ChainReader2::GetFirst_CAL(const string &FN)
	{
	m_LR.Open(FN);
	bool Ok = m_LR.ReadLine(m_Line);
	if (!Ok)
		Die("Failed to read first line of CAL file '%s'",
		  FN.c_str());
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

PDBChain *ChainReader2::GetFirst_PDB(const string &FN)
	{
	ReadLinesFromFile(FN, m_Lines);
	string FallbackLabel;
	GetFallbackLabelFromFN(FN, FallbackLabel);
	PDBChain::ChainsFromLines_PDB(m_Lines, m_Chains_PDB, FallbackLabel);
	m_ChainIdx_PDB = 0;
	return GetNext_PDB();
	}

PDBChain *ChainReader2::GetFirst_CIF(const string &FN)
	{
	ReadLinesFromFile(FN, m_Lines);
	string FallbackLabel;
	GetFallbackLabelFromFN(FN, FallbackLabel);
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
