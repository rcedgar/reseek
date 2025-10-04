#include "myutils.h"
#include "chainreader2.h"
#include "timing.h"

uint ChainReader2::m_CRGlobalChainCount;
uint ChainReader2::m_CRGlobalFormatErrors;

void ChainReader2::Close()
	{
	m_CRGlobalLock.lock();
	if (m_Trace) Log("ChainReader2::Close()\n");
	if (m_State != STATE_Closed)
		{
		m_State = STATE_Closed;
		if (m_ptrFS != 0)
			delete m_ptrFS;
		m_ptrFS = 0;
		}
	m_CRGlobalLock.unlock();
	}

void ChainReader2::Open(const string &FileName)
	{
	asserta(m_State == STATE_Closed);
	asserta(m_ptrFS == 0);
	PDBFileScanner *FS = new PDBFileScanner;
	FS->Open(FileName);
	Open(*FS);
	}

void ChainReader2::Open(PDBFileScanner &FS)
	{
	asserta(m_State == STATE_Closed);
	m_ptrFS = &FS;
	m_Trace = opt(trace_chainreader2);
	m_ptrFS->m_Trace = opt(trace_chainreader2);
	if (m_Trace) Log("ChainReader2::Open()\n");
	m_State = STATE_PendingFile;
	m_CRGlobalChainCount = 0;
	}

void ChainReader2::Open(vector<PDBChain *> &Chains)
	{
	asserta(m_State == STATE_Closed);
	m_ptrChains = &Chains;
	m_ChainIdx_Vec = 0;
	}

// Files first, then directories to reduce queue
PDBChain *ChainReader2::GetFirst(const string &FN)
	{
	m_CurrentFN = FN;

	string Ext;
	GetExtFromPathName(FN, Ext);
	ToLower(Ext);
	if (m_Trace) Log("  FN=%s Ext=%s\n", FN.c_str(), Ext.c_str());
	if (Ext == "cal")
		{
		m_State = STATE_ReadingCALFile;
		PDBChain *Chain = GetFirst_CAL(FN);
		if (Chain != 0)
			return Chain;
		}
	else if (Ext == "bca")
		{
		m_State = STATE_ReadingBCAFile;
		PDBChain *Chain = GetFirst_BCA(FN);
		if (Chain != 0)
			return Chain;
		}
	else if (Ext == "pdb" || Ext == "pdb.gz" || Ext == "ent" || Ext == "ent.gz")
		{
		m_State = STATE_ReadingPDBFile;
		PDBChain *Chain = GetFirst_PDB(FN);
		if (Chain != 0)
			return Chain;
		}
	else if (Ext == "cif" || Ext == "cif.gz" || Ext == "mmcif" || Ext == "mmcif.gz")
		{
		m_State = STATE_ReadingCIFFile;
		PDBChain *Chain = GetFirst_CIF(FN);
		if (Chain != 0)
			return Chain;
		}
	else
		Die("ChainReader2::GetNext(%s), unknown extension", FN.c_str());
	return 0;
	}

PDBChain *ChainReader2::GetNext()
	{
	StartTimer(ChainReader2_GetNext);
	for (uint SanityCounter = 0; ; ++SanityCounter)
		{
		if (SanityCounter > 100)
			Warning("Excessive looping in ChainReader2::GetNext()");

		m_CRPerThreadLock.lock();
		EndTimer(ChainReader2_GetNext);
		StartTimer(ChainReader2_GetNextLo1);
		PDBChain *Chain = GetNextLo1();
		EndTimer(ChainReader2_GetNextLo1);
		StartTimer(ChainReader2_GetNext);
		m_CRPerThreadLock.unlock();

		if (Chain == 0)
			{
			EndTimer(ChainReader2_GetNext);
			return 0;
			}

		if (Chain->GetSeqLength() == 0)
			{
			delete Chain;
			continue;
			}
		m_CRGlobalLock.lock();
		Chain->m_Idx = m_CRGlobalChainCount++;
		m_CRGlobalLock.unlock();
		EndTimer(ChainReader2_GetNext);
		return Chain;
		}
	}

PDBChain *ChainReader2::GetNextLo1()
	{
	for (uint SanityCounter = 0; SanityCounter < 100; ++SanityCounter)
		{
		if (m_Trace) Log("GetNextLo1() state=%d\n", m_State);
		switch (m_State)
			{
		case STATE_Closed:
			return 0;

		case STATE_PendingFile:
			{
			string FN;
			bool Ok = m_ptrFS->GetNext(FN);
			if (!Ok)
				return 0;
			PDBChain *Chain = GetFirst(FN);
			if (Chain != 0)
				return Chain;
			continue;
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

		case STATE_ReadingBCAFile:
			{
			PDBChain *Chain = GetNext_BCA();
			if (Chain != 0)
				return Chain;
			if (m_Trace) Log("GetNext_BCA()=0, state->PendingFile\n");
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
	Die("Excessive looping in ChainReader2::GetNextLo1()");
	return 0;
	}

void ChainReader2::GetFallbackLabelFromFN(const string &FN, string &Label)
	{
	GetStemName(FN, Label);
	string Ext;
	GetExtFromPathName(FN, Ext);
	ToLower(Ext);

// Special-case for downloaded PDB files e.g. pdb1iv1.ent
	if (Ext == "pdb" || Ext == "ent" || Ext == "pdb.gz" || Ext == "ent.gz")
		{
		if (Label.size() == 7 && Label[0] == 'p' && Label[1] == 'd' && Label[2] == 'b')
			{
			Label = Label.substr(3, string::npos);
			ToUpper(Label);
			}
		}
	}

PDBChain *ChainReader2::GetFirst_BCA(const string &FN)
	{
	m_BCA.Open(FN);
	m_ChainIdx_BCA = 0;
	return GetNext_BCA();
	}

PDBChain *ChainReader2::GetNext_BCA()
	{
	uint64 ChainCount = m_BCA.GetChainCount();
	if (m_ChainIdx_BCA >= ChainCount)
		{
		m_BCA.Close();
		return 0;
		}
	PDBChain *Chain = new PDBChain;
	m_BCA.ReadChain(m_ChainIdx_BCA++, *Chain);
	return Chain;
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
		  m_CurrentFN.c_str());

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
	for (uint LineNr = 0; LineNr < N; ++LineNr)
		{
		const string &Line = m_Lines[LineNr];
		if (Line.empty())
			continue;
		Split(Line, Fields, '\t');
		if (Fields.size() != 4 || Fields[0].size() != 1)
			Die("%s: Invalid CAL record '%s'",
			  m_CurrentFN.c_str(), Line.c_str());

		char aa = Fields[0][0];
		float X = StrToFloatf(Fields[1]);
		float Y = StrToFloatf(Fields[2]);
		float Z = StrToFloatf(Fields[3]);

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
	string Label;
	GetFallbackLabelFromFN(FN, Label);
	ChainsFromLines_PDB(m_Lines, m_Chains_PDB, Label);
	m_ChainIdx_PDB = 0;
	return GetNext_PDB();
	}

PDBChain *ChainReader2::GetFirst_CIF(const string &FN)
	{
	ReadLinesFromFile(FN, m_Lines);

	string FallbackLabel;
	GetFallbackLabelFromFN(FN, FallbackLabel);
	ChainsFromLines_CIF(m_Lines, m_Chains_CIF, FallbackLabel);
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

PDBChain *ChainReader2::GetNext_Vec()
	{
	asserta(m_ptrChains != 0);
	const uint N = SIZE(*m_ptrChains);
	if (m_ChainIdx_Vec == N)
		return 0;
	asserta(m_ChainIdx_CIF < N);
	PDBChain *Chain = (*m_ptrChains)[m_ChainIdx_Vec++];
	if (m_Trace) Log("ChainReader2::GetNext_Vec() %u/%u Label=%s\n", m_ChainIdx_Vec, N, Chain->m_Label.c_str());
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

PDBChain *ChainReader2::ChainFromLines_CAL(const vector<string> &Lines) const
	{
	PDBChain *Chain = new PDBChain;
	Chain->FromCalLines(Lines);
	return Chain;
	}

bool ChainReader2::IsATOMLine_PDB(const string &Line) const
	{
	if (SIZE(Line) < 27)
		return false;
	if (strncmp(Line.c_str(), "ATOM  ", 6) == 0)
		return true;
	return false;
	}

bool ChainReader2::IsChainEndLine_PDB(const string &Line) const
	{
	if (StartsWith(Line, "TER ") || StartsWith(Line, "ENDMDL"))
		return true;
	return false;
	}
