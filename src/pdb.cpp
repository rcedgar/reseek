#include "myutils.h"
#include "chainreader2.h"

static bool LabelAlreadyHasChain(const string &Label, 
  const string &ChainStr)
	{
	uint chn = SIZE(ChainStr);
	if (chn != 1)
		return false;
	uint labn = SIZE(Label);
	if (labn < 6)
		return false;
	if (tolower(Label[labn-1]) != tolower(ChainStr[chn-1]))
		return false;
	char c = Label[labn-2];
	if (c == '_' || c == ':' || c == '.')
		return true;
	return false;
	}

void ChainizeLabel(string &Label, const string &_ChainStr)
	{
	if (opt_nochainchar)
		return;
	string ChainStr = _ChainStr;
	if (ChainStr == "" || ChainStr == " ")
		ChainStr = '_';
	if (LabelAlreadyHasChain(Label, _ChainStr))
		return;
	Label += (optset_chainsep ? string(opt_chainsep) : "_");
	Label += ChainStr;
	}

bool PDBChain::FromPDBLines(const string &Label,
  const vector<string> &Lines, bool SaveLines)
	{
	Clear();
	if (SaveLines)
		m_Lines = Lines;
	m_Label = Label;
	const uint N = SIZE(Lines);
	uint ResidueCount = 0;
	int CurrentResidueNumber = INT_MAX;
	vector<string> ATOMLines;
	string ChainStr;
	for (uint LineNr = 0; LineNr < N; ++LineNr)
		{
		const string &Line = Lines[LineNr];
	// Can be multiple models for same chain, use first only
		if (StartsWith(Line, "TER ") || StartsWith(Line, "ENDMDL"))
			break;
		const size_t L = Line.size();

		char LineChainChar = Line[21];
		string LineChainStr;
		LineChainStr.push_back(LineChainChar);
		if (ChainStr == "")
			ChainStr = LineChainStr;
		else if (ChainStr != LineChainStr)
			Die("PDBChain::FromPDBLines() two chains %s, %s",
			  ChainStr.c_str(), LineChainStr.c_str());

		char aa;
		float X, Y, Z;
		bool IsCA = GetFieldsFromATOMLine(Line, X, Y, Z, aa);
		if (!IsCA)
			continue;

		m_Seq.push_back(aa);
		m_Xs.push_back(X);
		m_Ys.push_back(Y);
		m_Zs.push_back(Z);
		}

	ChainizeLabel(m_Label, ChainStr);
	bool Ok = (SIZE(m_Xs) > 0);
	return Ok;
	}

void ChainReader2::ChainsFromLines_PDB(const vector<string> &Lines,
  vector<PDBChain *> &Chains, const string &Label) const
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
			if (Line.size() < 57)
				continue;
			char ChainChar = Line[21];
			if (ChainChar != CurrChainChar)
				{
				if (AnyAtoms && !ChainLines.empty())
					{
					PDBChain *Chain = new PDBChain;
					string ChainStr;
					bool Ok = Chain->FromPDBLines(Label, ChainLines, m_SaveLines);
					if (Ok)
						Chains.push_back(Chain);
					else
						delete Chain;
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
		PDBChain *Chain = new PDBChain;
		bool Ok = Chain->FromPDBLines(Label, ChainLines, m_SaveLines);
		ChainLines.clear();
		Chains.push_back(Chain);
		}
	}
