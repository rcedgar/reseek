#include "myutils.h"
#include "pdbchain.h"
#include "calreader.h"

void ReadLinesFromFile(const string &FileName, vector<string> &Lines);

// $d/data/pdb/a0/pdb1a04.ent
void GetLabelFromFileNamePDBEnt(const string &FileName, string &Label)
	{
	vector<string> Fields;
	Split(FileName, Fields, '/');
	uint FieldCount = SIZE(Fields);
	asserta(FieldCount > 0);
	Label = Fields[FieldCount-1];

	size_t n = Label.size();
	if (n < 4)
		return;

	if (EndsWith(FileName, ".pdb") || EndsWith(FileName, ".ent"))
		Label.resize(n-4);
	if (Label[0] == 'p' && Label[1] == 'd' && Label[2] == 'b')
		Label = Label.substr(3, string::npos);

	for (uint i = 0; i < SIZE(Label); ++i)
		Label[i] = toupper(Label[i]);
	}

void GetLabelFromFileName(const string &FileName, string &Label)
	{
	if (EndsWith(FileName, ".ent"))
		{
		GetLabelFromFileNamePDBEnt(FileName, Label);
		return;
		}

	vector<string> Fields;
	Split(FileName, Fields, '/');
	uint FieldCount = SIZE(Fields);
	asserta(FieldCount > 0);
	Label = Fields[FieldCount-1];
	const uint n = SIZE(Label);

	if (n > 4 && EndsWith(Label, ".pdb"))
		Label = Label.substr(0, n-4);
	else if (n == 6 && Label[4] == '_')
		Label = Label.substr(0, 4);
	}

void ReadChains(const vector<string> &FileNames,
  vector<PDBChain *> &Structures)
	{
	Structures.clear();
	const uint N = SIZE(FileNames);
	for (uint i = 0; i < N; ++i)
		{
		const string &FileName = FileNames[i];
		if (EndsWith(FileName, ".gz"))
			Die("Compressed files not supported '%s'", FileName.c_str());

		string Label;
		GetLabelFromFileName(FileName, Label);

		vector<string> Lines;
		ReadLinesFromFile(FileName, Lines);

		ProgressStep(i, N, "Reading chains %s", Label.c_str());

		vector<PDBChain *> Chains;
		PDBChain::ChainsFromLines(Label, Lines, Chains);
		const uint NC = SIZE(Chains);
		for (uint j = 0; j < NC; ++j)
			Structures.push_back(Chains[j]);
		}
	}

void ReadChainsCal(const string &FileName, vector<PDBChain *> &Structures)
	{
	CalReader CR;
	CR.Open(FileName);
	uint N = 0;
	for (;;)
		{
		PDBChain &Q = *new PDBChain;
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;
		Structures.push_back(&Q);
		++N;
#if DEBUG
		if (N%1000 == 0)
#else
		if (N%100000 == 0)
#endif
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("Reading %s %s%% done\r",
			  FileName.c_str(), sPct.c_str());
			}
		}
	Progress("Reading %s 100.0%% done\n", FileName.c_str());
	}

void ReadChains(const string &FileName, vector<PDBChain *> &Structures)
	{
	if (FileName.empty())
		Die("Missing chains filename");

	if (EndsWith(FileName, ".cal") || EndsWith(FileName, ".ppc"))
		{
		ReadChainsCal(FileName, Structures);
		return;
		}
	else if (EndsWith(FileName, ".files"))
		{
		vector<string> FileNames;
		ReadLinesFromFile(FileName, FileNames);
		ReadChains(FileNames, Structures);
		}

	string Label;
	GetLabelFromFileName(FileName, Label);

	vector<string> Lines;
	ReadLinesFromFile(FileName, Lines);

	vector<PDBChain *> Chains;
	PDBChain::ChainsFromLines(Label, Lines, Chains);
	const uint NC = SIZE(Chains);
	for (uint j = 0; j < NC; ++j)
		Structures.push_back(Chains[j]);
	}
