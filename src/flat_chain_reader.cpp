#include "myutils.h"
#include "flat_chain.h"
#include "flat_chain_reader.h"

uint flat_chain_reader::m_CRGlobalChainCount;
uint flat_chain_reader::m_CRGlobalFormatErrors;

void flat_chain_reader::Close()
	{
	m_CRGlobalLock.lock();
	if (m_Trace) Log("flat_chain_reader::Close()\n");
	if (m_State != STATE_Closed)
		{
		m_State = STATE_Closed;
		if (m_ptrFS != 0)
			delete m_ptrFS;
		m_ptrFS = 0;
		}
	m_CRGlobalLock.unlock();
	}

void flat_chain_reader::Open(const string &FileName)
	{
	asserta(m_State == STATE_Closed);
	asserta(m_ptrFS == 0);
	PDBFileScanner *FS = new PDBFileScanner;
	FS->Open(FileName);
	Open(*FS);
	}

void flat_chain_reader::Open(PDBFileScanner &FS)
	{
	asserta(m_State == STATE_Closed);
	m_ptrFS = &FS;
	m_Trace = opt(trace_chainreader2);
	m_ptrFS->m_Trace = opt(trace_chainreader2);
	if (m_Trace) Log("flat_chain_reader::Open()\n");
	m_State = STATE_PendingFile;
	m_CRGlobalChainCount = 0;
	}

void flat_chain_reader::Open(vector<flat_chain_t *> &Chains)
	{
	asserta(m_State == STATE_Closed);
	m_ptrChains = &Chains;
	m_ChainIdx_Vec = 0;
	}

// Files first, then directories to reduce queue
flat_chain_t* flat_chain_reader::GetFirst(const string &FN)
	{
	m_CurrentFN = FN;

	string Ext;
	GetExtFromPathName(FN, Ext);
	ToLower(Ext);
	if (m_Trace) Log("  FN=%s Ext=%s\n", FN.c_str(), Ext.c_str());
	if (Ext == "cal")
		{
		m_State = STATE_ReadingCALFile;
		flat_chain_t* Chain = GetFirst_CAL(FN);
		if (Chain)
			return Chain;
		}
	else if (Ext == "bca")
		{
		m_State = STATE_ReadingBCAFile;
		flat_chain_t* Chain = GetFirst_BCA(FN);
		if (Chain)
			return Chain;
		}
	else if (Ext == "pdb" || Ext == "pdb.gz" || Ext == "ent" || Ext == "ent.gz")
		{
		m_State = STATE_ReadingPDBFile;
		flat_chain_t* Chain = GetFirst_PDB(FN);
		if (Chain)
			return Chain;
		}
	else if (Ext == "cif" || Ext == "cif.gz" || Ext == "mmcif" || Ext == "mmcif.gz")
		{
		m_State = STATE_ReadingCIFFile;
		flat_chain_t* Chain = GetFirst_CIF(FN);
		if (Chain)
			return Chain;
		}
	else
		Die("flat_chain_reader::GetNext(%s), unknown extension", FN.c_str());
	return 0;
	}

flat_chain_t* flat_chain_reader::GetNext()
	{
	for (uint SanityCounter = 0; ; ++SanityCounter)
		{
		if (SanityCounter > 100)
			Warning("Excessive looping in flat_chain_reader::GetNext()");

		m_CRPerThreadLock.lock();
		flat_chain_t* Chain = GetNextLo1();
		m_CRPerThreadLock.unlock();

		if (!Chain)
			return 0;

		if (Chain->get_length() == 0)
			continue;
		return Chain;
		}
	}

flat_chain_t* flat_chain_reader::GetNextLo1()
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
			flat_chain_t* Chain = GetFirst(FN);
			if (Chain)
				return Chain;
			continue;
			}

		case STATE_ReadingCALFile:
			{
			flat_chain_t* Chain = GetNext_CAL();
			if (Chain)
				return Chain;
			if (m_Trace) Log("GetNext_CAL()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingBCAFile:
			{
			flat_chain_t* Chain = GetNext_BCA();
			if (Chain)
				return Chain;
			if (m_Trace) Log("GetNext_BCA()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingPDBFile:
			{
			flat_chain_t* Chain = GetNext_PDB();
			if (Chain)
				return Chain;
			if (m_Trace) Log("GetNext_PDB()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		case STATE_ReadingCIFFile:
			{
			flat_chain_t* Chain = GetNext_CIF();
			if (Chain)
				return Chain;
			if (m_Trace) Log("GetNext_CIF()=0, state->PendingFile\n");
			m_State = STATE_PendingFile;
			continue;
			}

		default:
			asserta(false);
			}
		}
	Die("Excessive looping in flat_chain_reader::GetNextLo1()");
	return 0;
	}

flat_chain_t* flat_chain_reader::GetFirst_BCA(const string &FN)
	{
	m_BCA.Open(FN);
	m_ChainIdx_BCA = 0;
	return GetNext_BCA();
	}

flat_chain_t* flat_chain_reader::GetNext_BCA()
	{
	uint64 ChainCount = m_BCA.GetChainCount();
	if (m_ChainIdx_BCA >= ChainCount)
		{
		m_BCA.Close();
		return 0;
		}
	flat_chain_t* chain = m_BCA.read_flat_chain(m_ChainIdx_BCA++);
	return chain;
	}

flat_chain_t* flat_chain_reader::GetFirst_CAL(const string &FN)
	{
	m_LR.Open(FN);
	bool Ok = m_LR.ReadLine(m_Line);
	if (!Ok)
		Die("Failed to read first line of CAL file '%s'",
		  FN.c_str());
	return GetNext_CAL();
	}

flat_chain_t* flat_chain_reader::GetNext_CAL()
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
	if (m_Trace) Log("flat_chain_reader::GetNext_CAL() Label=%s\n", Label.c_str());
	m_Lines.clear();
	while (m_LR.ReadLine(m_Line))
		{
		if (m_Line.c_str()[0] == '>')
			break;
		m_Lines.push_back(m_Line);
		}

/***
>102l
M       43.619  -1.924  8.869
N       40.445  -0.876  10.670
I       38.254  2.240   11.220
F       40.340  3.621   14.036
***/
	const uint N = SIZE(m_Lines);
	vector<string> Fields;
	vector<char> aas;
	vector<float> Xs, Ys, Zs;
	aas.reserve(RESERVE_CHAIN_LENGTH);
	Xs.reserve(RESERVE_CHAIN_LENGTH);
	Ys.reserve(RESERVE_CHAIN_LENGTH);
	Zs.reserve(RESERVE_CHAIN_LENGTH);
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

		aas.push_back(aa);
		Xs.push_back(X);
		Ys.push_back(Y);
		Zs.push_back(Z);
		}
	auto chain = flat_chain_t::newflat(Label, aas, Xs, Ys, Zs);
	return chain;
	}

flat_chain_t* flat_chain_reader::GetFirst_PDB(const string &FN)
	{
	ReadLinesFromFile(FN, m_Lines);
	string Label;
	GetFallbackLabelFromFN(FN, Label);
	ChainsFromLines_PDB(m_Lines, m_Chains_PDB, Label);
	m_ChainIdx_PDB = 0;
	return GetNext_PDB();
	}

flat_chain_t* flat_chain_reader::GetFirst_CIF(const string &FN)
	{
	ReadLinesFromFile(FN, m_Lines);

	string FallbackLabel;
	GetFallbackLabelFromFN(FN, FallbackLabel);
	ChainsFromLines_CIF(m_Lines, m_Chains_CIF, FallbackLabel);
	m_ChainIdx_CIF = 0;
	return GetNext_CIF();
	}

flat_chain_t* flat_chain_reader::GetNext_PDB()
	{
	const uint N = SIZE(m_Chains_PDB);
	if (m_ChainIdx_PDB == N)
		return 0;
	asserta(m_ChainIdx_PDB < N);
	flat_chain_t* Chain = m_Chains_PDB[m_ChainIdx_PDB++];
	if (m_Trace) Log("flat_chain_reader::GetNext_PDB() %u/%u Label=%s\n", m_ChainIdx_PDB, N, Chain->m_label.c_str());
	return Chain;
	}

flat_chain_t* flat_chain_reader::GetNext_Vec()
	{
	asserta(m_ptrChains != 0);
	const uint N = SIZE(*m_ptrChains);
	if (m_ChainIdx_Vec == N)
		return 0;
	asserta(m_ChainIdx_CIF < N);
	flat_chain_t* Chain = (*m_ptrChains)[m_ChainIdx_Vec++];
	if (m_Trace) Log("flat_chain_reader::GetNext_Vec() %u/%u Label=%s\n", m_ChainIdx_Vec, N, Chain->m_label.c_str());
	return Chain;
	}

flat_chain_t* flat_chain_reader::GetNext_CIF()
	{
	const uint N = SIZE(m_Chains_CIF);
	if (m_ChainIdx_CIF == N)
		return 0;
	asserta(m_ChainIdx_CIF < N);
	flat_chain_t* Chain = m_Chains_CIF[m_ChainIdx_CIF++];
	if (m_Trace) Log("flat_chain_reader::GetNext_CIF() %u/%u Label=%s\n", m_ChainIdx_CIF, N, Chain->m_label.c_str());
	return Chain;
	}

bool flat_chain_reader::IsATOMLine_PDB(const string &Line) const
	{
	if (SIZE(Line) < 27)
		return false;
	if (strncmp(Line.c_str(), "ATOM  ", 6) == 0)
		return true;
	return false;
	}

bool flat_chain_reader::IsChainEndLine_PDB(const string &Line) const
	{
	if (StartsWith(Line, "TER ") || StartsWith(Line, "ENDMDL"))
		return true;
	return false;
	}

void flat_chain_reader::ChainsFromLines_PDB(const vector<string> &Lines,
  vector<flat_chain_t *> &Chains, const string &Label) const
	{
	Chains.clear();
	const uint N = SIZE(Lines);
	vector<string> ChainLines;
	char CurrChainChar = 0;
	bool AnyAtoms = false;
	bool EndOfChainFound = false;
	for (uint i = 0; i < N; ++i)
		{
		const string &Line = Lines[i];
		if (IsChainEndLine_PDB(Line))
			EndOfChainFound = true;
		if (IsATOMLine_PDB(Line))
			{
			if (Line.size() < 54)
				continue;
			char ChainChar = Line[21];
			if (ChainChar != CurrChainChar)
				{
				if (AnyAtoms && !ChainLines.empty())
					{
					flat_chain_t* Chain = flat_chain_t::newflat(0);
					string ChainStr;
					bool Ok = Chain->from_pdb_lines(Label, ChainLines, m_SaveLines);
					if (Ok)
						Chains.push_back(Chain);
					ChainLines.clear();
					EndOfChainFound = false;
					AnyAtoms = false;
					}
				CurrChainChar = ChainChar;
				}
			if (!EndOfChainFound)
				ChainLines.push_back(Line);
			AnyAtoms = true;
			}
		}

	if (!ChainLines.empty() && AnyAtoms)
		{
		flat_chain_t* Chain = flat_chain_t::newflat(0);
		bool Ok = Chain->from_pdb_lines(Label, ChainLines, m_SaveLines);
		ChainLines.clear();
		Chains.push_back(Chain);
		}
	}
