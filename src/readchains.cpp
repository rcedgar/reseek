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
	else if (n > 4 && EndsWith(Label, ".cif"))
		Label = Label.substr(0, n-4);
	else if (n == 6 && Label[4] == '_')
		Label = Label.substr(0, 4);
	}

void ReadChainsFromFileNameVec(const vector<string> &FileNames,
  vector<PDBChain *> &Chains)
	{
	Chains.clear();
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
			Chains.push_back(Chains[j]);
		}
	}

static void ReadChainsFromCalFile(const string &FileName, vector<PDBChain *> &Chains)
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
		Chains.push_back(&Q);
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

static void ReadChainsFromCIFFile(const string &FileName, 
  vector<PDBChain *> &Chains)
	{
	Die("TODO");
	}

static void ReadChainsFromPDBFile(const string &FileName, 
  vector<PDBChain *> &Chains)
	{
	string Label;
	GetLabelFromFileName(FileName, Label);

	vector<string> Lines;
	ReadLinesFromFile(FileName, Lines);

	PDBChain::ChainsFromLines(Label, Lines, Chains);
	const uint NC = SIZE(Chains);
	for (uint j = 0; j < NC; ++j)
		Chains.push_back(Chains[j]);
	}

void ReadChainsFromDirectory(const string &DirName,
  vector<PDBChain *> &Chains, bool Recursive)
	{
	vector<string> FileNames;
	vector<bool> IsSubDirs;
	mylistdir(DirName, FileNames, IsSubDirs);
	const uint N = SIZE(FileNames);
	asserta(SIZE(IsSubDirs) == N);
	vector<string> SubDirs;
	for (uint i = 0; i < N; ++i)
		{
		const string &FileName = FileNames[i];
		bool IsSubDir = IsSubDirs[i];
		if (IsSubDir)
			{
			SubDirs.push_back(FileName);
			continue;
			}
		ReadChains(FileName, Chains);
		}

	if (Recursive)
		{
		const uint M = SIZE(SubDirs);
		for (uint i = 0; i < M; ++i)
			{
			string Dir = DirName + string("/") + SubDirs[i] + string("/");
			ReadChains(Dir, Chains);
			}
		}
	}

void ReadChains(const string &FileName, vector<PDBChain *> &Chains)
	{
	if (FileName.empty())
		Die("Missing chains filename");

	if (EndsWith(FileName, "/"))
		{
		ReadChainsFromDirectory(FileName, Chains, false);
		return;
		}

	string Ext;
	GetExtFromPathName(FileName, Ext);
	ToLower(Ext);

	if (EndsWith(FileName, ".files"))
		{
		vector<string> FileNames;
		ReadLinesFromFile(FileName, FileNames);
		ReadChainsFromFileNameVec(FileNames, Chains);
		return;
		}

	if (Ext == ".cal")
		ReadChainsFromCalFile(FileName, Chains);
	else if (Ext == ".pdb" || Ext == ".ent")
		ReadChainsFromPDBFile(FileName, Chains);
	else if (Ext == ".cif" || Ext == ".mmcif")
		ReadChainsFromCIFFile(FileName, Chains);
	}
