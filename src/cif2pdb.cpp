#include "myutils.h"
#include "pdbchain.h"
#include "chainreader2.h"
#include "cif.h"
#include <map>

static uint GetFieldIdx(const map<string, uint> &FieldToIdx,
						const string &Name)
	{
	map<string, uint>::const_iterator iter = FieldToIdx.find(Name);
	if (iter == FieldToIdx.end())
		Die("CIF field not found '%s'", Name.c_str());
	uint Idx = iter->second;
	return Idx;
	}

static void MakePDBAtomLine(const string &ATOM, const string &ChainStr,
							uint AtomNr, uint ResNr, const string &atomname_Fld,
							const string &aa_Fld, float X, float Y, float Z, 
							string &Line)
	{
	if (ATOM == "HETATM")
		Line = ATOM;
	else
		Line = "ATOM  ";
	assert(SIZE(Line) == 6);

	assert(AtomNr <= 999999);
	Psa(Line, "%-6u", AtomNr);
	assert(SIZE(Line) == 12);

	Psa(Line, "%-4.4s", atomname_Fld.c_str());
	assert(SIZE(Line) == 16);

	Line += ' ';	// altLoc
	assert(SIZE(Line) == 17);

	Line += aa_Fld;
	assert(SIZE(Line) == 20);

	if (SIZE(ChainStr) == 0)
		{
		Line += ' ';
		Line += 'A';
		}
	else if (SIZE(ChainStr) == 1)
		{
		Line += ' ';
		Line += ChainStr[0];
		}
	else if (SIZE(ChainStr) == 2)
		Line += ChainStr;
	else
		{
		Line += ChainStr[0];
		Line += ChainStr[1];
		}
	assert(SIZE(Line) == 22);

	asserta(ResNr <= 9999);
	Psa(Line, "%4d", ResNr);
	assert(SIZE(Line) == 26);

	Line += ' '; // insertion code
	assert(SIZE(Line) == 27);
	
	Line += "   ";
	assert(SIZE(Line) == 30);

	Psa(Line, "%8.3f", X);
	Psa(Line, "%8.3f", Y);
	Psa(Line, "%8.3f", Z);
	assert(SIZE(Line) == 54);
	}

void ReadCIF(const vector<string> &Lines,
			 vector<vector<string> > &PDBAtomLinesVec)
	{
	PDBAtomLinesVec.clear();
	string CurrentChainStr;
	PDBChain *Chain = 0;
	const uint N = SIZE(Lines);
	CIF_PARSER_STATE PS = PS_WaitingForLoop;
	vector<string> FieldList;
	vector<string> ATOMLines;
	for (uint i = 0; i < N; ++i)
		{
		const string &Line = Lines[i];
		if (!Line.empty() && Line[0] == '#')
			continue;
		if (PS == PS_Finished)
			break;

		switch (PS)
			{
		case PS_WaitingForLoop:
			{
			if (Line == "loop_")
				{
				PS = PS_AtLoop;
				continue;
				}
			break;
			}

		case PS_AtLoop:
			{
			if (StartsWith(Line, "_atom_site."))
				{
				PS = PS_InFieldList;
				string Name = Line;
				StripWhiteSpace(Name);
				FieldList.push_back(Name);
				}
			else
				PS = PS_WaitingForLoop;
			break;
			}

		case PS_InFieldList:
			{
			if (StartsWith(Line, "_atom_site."))
				{
				string Name = Line;
				StripWhiteSpace(Name);
				FieldList.push_back(Name);
				}
			else if (Line == "loop_")
				PS = PS_AtLoop;
			else if (StartsWith(Line, "ATOM ") ||
					 StartsWith(Line, "HETATM"))
				{
				PS = PS_InATOMs;
				ATOMLines.push_back(Line);
				}
			break;
			}

		case PS_InATOMs:
			{
			if (StartsWith(Line, "ATOM ") ||
				StartsWith(Line, "HETATM"))
				ATOMLines.push_back(Line);
			else
				PS = PS_Finished;
			break;
			}

		default:
			asserta(false);
			}
		}

	const uint FieldCount = SIZE(FieldList);
	map<string, uint> FieldToIdx;
	for (uint Idx = 0; Idx < FieldCount ; ++Idx)
		FieldToIdx[FieldList[Idx]] = Idx;

	uint Chain_FldIdx = GetFieldIdx(FieldToIdx, "_atom_site.auth_asym_id");
	uint CA_FldIdx = GetFieldIdx(FieldToIdx, "_atom_site.label_atom_id");
	uint ResNr_FldIdx = GetFieldIdx(FieldToIdx, "_atom_site.label_seq_id");
	uint AtomNr_FldIdx = GetFieldIdx(FieldToIdx, "_atom_site.id");
	uint X_FldIdx = GetFieldIdx(FieldToIdx, "_atom_site.Cartn_x");
	uint Y_FldIdx = GetFieldIdx(FieldToIdx, "_atom_site.Cartn_y");
	uint Z_FldIdx = GetFieldIdx(FieldToIdx, "_atom_site.Cartn_z");
	uint aa_FldIdx = GetFieldIdx(FieldToIdx, "_atom_site.label_comp_id");
	uint ModelNr_FldIdx = GetFieldIdx(FieldToIdx, "_atom_site.pdbx_PDB_model_num");

	if (Chain_FldIdx == UINT_MAX) return;
	if (AtomNr_FldIdx == UINT_MAX) return;
	if (ResNr_FldIdx == UINT_MAX) return;
	if (CA_FldIdx == UINT_MAX) return;
	if (X_FldIdx == UINT_MAX) return;
	if (Y_FldIdx == UINT_MAX) return;
	if (Z_FldIdx == UINT_MAX) return;
	if (aa_FldIdx == UINT_MAX) return;
	uint CurrentModelNr = UINT_MAX;

	vector<string> PDBAtomLines;
	const uint ATOMLineCount = SIZE(ATOMLines);
	for (uint i = 0; i < ATOMLineCount; ++i)
		{
		const string &Line = ATOMLines[i];
		vector<string> Fields;
		SplitWhite(Line, Fields);
		uint n = SIZE(Fields);
		if (n != FieldCount)
			Die("Expected %u fields got %u in '%s'\n",
			  FieldCount, n, Line.c_str());

		const string &CA_Fld = Fields[CA_FldIdx];

		if (ModelNr_FldIdx != UINT_MAX)
			{
			uint ModelNr = (uint) atoi(Fields[ModelNr_FldIdx].c_str());
			if (CurrentModelNr != UINT_MAX && ModelNr != CurrentModelNr)
				break;
			CurrentModelNr = ModelNr;
			}

		const string &Chain_Fld = Fields[Chain_FldIdx];
		string ChainStr = Chain_Fld;
		if (ChainStr == "")
			ChainStr = "__";
		if (ChainStr != CurrentChainStr)
			{
			if (!PDBAtomLines.empty())
				PDBAtomLinesVec.push_back(PDBAtomLines);
			CurrentChainStr = ChainStr;
			}

		const string &aa_Fld = Fields[aa_FldIdx];
		if (SIZE(aa_Fld) != 3)
			continue;
		uint AtomNr = StrToUint(Fields[AtomNr_FldIdx]);
		uint ResNr = StrToUint(Fields[ResNr_FldIdx]);
		float X = StrToFloatf(Fields[X_FldIdx]);
		float Y = StrToFloatf(Fields[Y_FldIdx]);
		float Z = StrToFloatf(Fields[Z_FldIdx]);
	//	char aa = GetOneFromThree(aa_Fld);

		string ATOM = "ATOM";
		if (StartsWith(Line, "HETATM"))
			ATOM = "HETATM";
		string PDBLine;
		MakePDBAtomLine(ATOM, ChainStr, AtomNr, ResNr, CA_Fld, aa_Fld, 
						X, Y, Z, PDBLine);
		PDBAtomLines.push_back(PDBLine);
		}
	if (!PDBAtomLines.empty())
		PDBAtomLinesVec.push_back(PDBAtomLines);
	}

void cmd_cif2pdb()
	{
	vector<string> CIFLines;
	ReadLinesFromFile(g_Arg1, CIFLines);

	vector<vector<string> > PDBAtomLinesVec;
	ReadCIF(CIFLines, PDBAtomLinesVec);

	if (!optset_output)
		return;

	FILE *f = CreateStdioFile(opt(output));

	const uint ChainCount = SIZE(PDBAtomLinesVec);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		vector<string> &Lines = PDBAtomLinesVec[ChainIdx];
		for (uint i = 0; i < SIZE(Lines); ++i)
			fprintf(f, "%s\n", Lines[i].c_str());
		}
	CloseStdioFile(f);
	}
