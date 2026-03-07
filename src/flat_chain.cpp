#include "myutils.h"
#include "pdbchain.h"
#include "pdbfilescanner.h"
#include "flat_chain.h"
#include "flat_chain_reader.h"

void ChainizeLabel(string &Label, const string &_ChainStr);

void flat_chain::set_xyz(const vector<float> &Xs,
	const vector<float> &Ys, const vector<float> &Zs)
	{
	const uint32_t L = SIZE(Xs);
	down0(m_xyz);
	m_xyz = create_chainxyz(L);
	for (uint32_t i = 0; i < L; ++i)
		{
		uint16_t ic_x = PDBChain::CoordToIC(Xs[i]);
		uint16_t ic_y = PDBChain::CoordToIC(Ys[i]);
		uint16_t ic_z = PDBChain::CoordToIC(Zs[i]);
		m_xyz->set(0, i, ic_x);
		m_xyz->set(1, i, ic_y);
		m_xyz->set(2, i, ic_z);
		}
	}

void flat_chain::set_aa(const vector<char> &aas)
	{
	const uint32_t L = SIZE(aas);
	down0(m_aa);
	m_aa = create_chainaa(L);
	memcpy(m_aa->m_data, aas.data(), L);
	}

bool flat_chain::from_pdb_lines(const string &label,
	const vector<string> &lines, bool save_lines)
	{
	clear();
	if (save_lines)
		m_lines = lines;
	m_label = label;
	const uint N = SIZE(lines);
	uint ResidueCount = 0;
	int CurrentResidueNumber = INT_MAX;
	string ChainStr;
	vector<char> aas;
	vector<float> Xs, Ys, Zs;
	aas.reserve(RESERVE_CHAIN_LENGTH);
	Xs.reserve(RESERVE_CHAIN_LENGTH);
	Ys.reserve(RESERVE_CHAIN_LENGTH);
	Zs.reserve(RESERVE_CHAIN_LENGTH);
	for (uint lineNr = 0; lineNr < N; ++lineNr)
		{
		const string &line = lines[lineNr];
	// Can be multiple models for same chain, use first only
		if (StartsWith(line, "TER ") || StartsWith(line, "ENDMDL"))
			break;
		const size_t L = line.size();

		char lineChainChar = line[21];
		string lineChainStr;
		lineChainStr.push_back(lineChainChar);
		if (ChainStr == "")
			ChainStr = lineChainStr;
		else if (ChainStr != lineChainStr)
			Die("flat_chain::from_pdb_lines() two chains %s, %s",
			  ChainStr.c_str(), lineChainStr.c_str());

		char aa;
		float X, Y, Z;
		bool IsCA = PDBChain::GetFieldsFromATOMLine(line, X, Y, Z, aa);
		if (!IsCA)
			continue;

		aas.push_back(aa);
		Xs.push_back(X);
		Ys.push_back(Y);
		Zs.push_back(Z);
		}
	set_xyz(Xs, Ys, Zs);
	set_aa(aas);

	ChainizeLabel(m_label, ChainStr);
	bool Ok = (SIZE(Xs) > 0);
	return Ok;
	}

void read_flat_chains(const string &fn, vector<flat_chain *> &chains)
	{
	PDBFileScanner FS;
	FS.Open(fn);

	flat_chain_reader CR;
	CR.Open(FS);
	for (;;)
		{
		flat_chain *chain = CR.GetNext();
		if (chain == 0)
			break;
		chains.push_back(chain);
		}
	}


void flat_chain::to_fasta(const string &fn) const
	{
	FILE *f = CreateStdioFile(fn);
	to_fasta(f);
	CloseStdioFile(f);
	}

void flat_chain::to_cal(const string &fn) const
	{
	FILE *f = CreateStdioFile(fn);
	to_cal(f);
	CloseStdioFile(f);
	}

void flat_chain::to_fasta(FILE *f) const
	{
	if (f == 0)
		return;
	string seq(m_aa->m_data, m_aa->m_size);
	SeqToFasta(f, m_label, seq);
	}

void flat_chain::to_cal(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, ">%s\n", m_label.c_str());
	uint L = get_length();
	for (uint i = 0; i < L; ++i)
		{
		char aa = get_aa(i);
		float x, y, z;
		get_coords(i, x, y, z);
		fprintf(f, "%c\t%.1f\t%.1f\t%.1f\n", aa, x, y, z);
		}
	}

void cmd_test()
	{
	vector<flat_chain *> chains;
	read_flat_chains(g_Arg1, chains);
	FILE *f = CreateStdioFile(opt(output));
	uint n = SIZE(chains);
	ProgressLog("%u chains\n", n);
	for (uint i = 0; i < n; ++i)
		chains[i]->to_cal(f);
	log_flat_stats();
	CloseStdioFile(f);
	}
