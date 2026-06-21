#include "myutils.h"
#include "flat_chain.h"
#include "flat_chain_reader.h"
#include "flat_params.h"
#include "flat_aligner.h"
#include "chain_data.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "kabsch.h"
#include "pdbfilescanner.h"
#include "flat_helpers.h"

static const uint32_t s_bits = bit_distmx | bit_mega_prof;

static void read_flat_chains_save_lines(const string &fn,
	vector<flat_chain_t *> &chains)
	{
	PDBFileScanner FS;
	FS.Open(fn);

	flat_chain_reader CR;
	CR.m_SaveLines = true;
	CR.Open(FS);
	for (;;)
		{
		flat_chain_t *chain = CR.GetNext();
		if (chain == 0)
			break;
		chains.push_back(chain);
		}
	}

static void check_chain_lengths(const vector<flat_chain_t *> &chains,
	const string &fn)
	{
	const uint n = uint(chains.size());
	for (uint i = 0; i < n; ++i)
		{
		const uint L = chains[i]->get_length();
		if (L == 0)
			Die("Empty chain in %s", fn.c_str());
		if (L > flat_params::m_maxL)
			Die("Chain %s length %u exceeds max %u",
				chains[i]->m_label.c_str(), L, flat_params::m_maxL);
		}
	}

static void build_chain_data_vec(const flat_params &params,
	const vector<flat_chain_t *> &chains,
	vector<chain_data *> &cdvec)
	{
	const uint n = uint(chains.size());
	cdvec.clear();
	cdvec.resize(n, 0);
	if (n == 0)
		return;

	chain_data **cd = myalloc(chain_data *, n);
	chain_data::fill_chain_data_vec(params, chains, s_bits, cd);
	for (uint i = 0; i < n; ++i)
		cdvec[i] = cd[i];
	myfree(cd);
	}

static void XformLine(const double t[3],
	const double u[3][3], string &Line)
	{
	float x, y, z;
	GetXYZFromATOMLine(Line, x, y, z);

	double Pt[3];
	double XPt[3];

	Pt[0] = x;
	Pt[1] = y;
	Pt[2] = z;
	transform(t, u, Pt, XPt);

	x = (float) XPt[0];
	y = (float) XPt[1];
	z = (float) XPt[2];
	SetXYZInATOMLine(Line, x, y, z, Line);
	}

static void XformLines(const double t[3],
	const double u[3][3], vector<string> &Lines)
	{
	const uint N = SIZE(Lines);
	for (uint i = 0; i < N; ++i)
		{
		string &Line = Lines[i];
		if (IsATOMLine(Line))
			XformLine(t, u, Line);
		}
	}

static void FlatKabsch(const flat_chain_t &chainQ,
	const flat_chain_t &chainT,
	uint loQ, uint loT,
	const char *path, uint ncol,
	double t[3], double u[3][3])
	{
	uint M = 0;
	for (uint col = 0; col < ncol; ++col)
		if (path[col] == 'M')
			++M;
	if (M == 0)
		Die("No aligned positions for superposition");

	double **x = myalloc(double *, M);
	double **y = myalloc(double *, M);
	uint posQ = loQ;
	uint posT = loT;
	uint m = 0;
	for (uint col = 0; col < ncol; ++col)
		{
		char c = path[col];
		if (c == 'M')
			{
			x[m] = myalloc(double, 3);
			y[m] = myalloc(double, 3);
			float xq, yq, zq;
			float xt, yt, zt;
			chainQ.get_coords(posQ, xq, yq, zq);
			chainT.get_coords(posT, xt, yt, zt);
			x[m][0] = xq;
			x[m][1] = yq;
			x[m][2] = zq;
			y[m][0] = xt;
			y[m][1] = yt;
			y[m][2] = zt;
			++m;
			++posQ;
			++posT;
			}
		else if (c == 'D')
			++posQ;
		else if (c == 'I')
			++posT;
		else
			asserta(false);
		}
	Kabsch(x, y, int(M), t, u);
	for (uint i = 0; i < M; ++i)
		{
		myfree(x[i]);
		myfree(y[i]);
		}
	myfree(x);
	myfree(y);
	}

static float AlignPairFlat(flat_aligner &fa,
	const chain_data &cdQ, const chain_data &cdT,
	bool do_output)
	{
	const flat_chain_t *chainQ = cdQ.m_chain;
	const flat_chain_t *chainT = cdT.m_chain;

	fa.cacheT(cdT.m_label, cdT.m_mega_prof, 0, cdT.m_L);
	fa.alignQ(cdQ.m_label, cdQ.m_mega_prof, 0, cdQ.m_L);
	const float score = fa.m_score;

	if (!do_output)
		return score;

	if (optset_aln)
		{
		FILE *f = CreateStdioFile(opt(aln));
		fa.write_aln(f);
		CloseStdioFile(f);
		}

	if (optset_output || optset_output2)
		{
		if (chainQ->m_lines.empty() || chainT->m_lines.empty())
			Die("-output/-output2 require PDB/CIF input with saved ATOM lines");

		double t[3];
		double u[3][3];
		FlatKabsch(*chainQ, *chainT,
			fa.m_loQ, fa.m_loT,
			fa.m_path_buffer, fa.m_ncol, t, u);

		vector<string> linesQ = chainQ->m_lines;
		XformLines(t, u, linesQ);

		if (optset_output)
			{
			FILE *f = CreateStdioFile(opt(output));
			for (uint i = 0; i < SIZE(linesQ); ++i)
				fprintf(f, "%s\n", linesQ[i].c_str());
			CloseStdioFile(f);
			}

		if (optset_output2)
			{
			FILE *f = CreateStdioFile(opt(output2));
			for (uint i = 0; i < SIZE(linesQ); ++i)
				{
				string line = linesQ[i];
				if (line.size() > 21)
					line[21] = '1';
				fprintf(f, "%s\n", line.c_str());
				}

			const vector<string> &linesT = chainT->m_lines;
			for (uint i = 0; i < SIZE(linesT); ++i)
				{
				string line = linesT[i];
				if (line.size() > 21)
					line[21] = '2';
				fprintf(f, "%s\n", line.c_str());
				}
			CloseStdioFile(f);
			}
		}

	return score;
	}

void cmd_flat_alignpair()
	{
	if (!optset_input2)
		Die("Must specify -input2");
	if (!optset_alphadir)
		Die("Must specify -alphadir");
	if (!optset_varstr)
		Die("Must specify -varstr");

	const string &qfn = g_Arg1;
	const string &tfn = opt(input2);

	vector<flat_chain_t *> chainsQ;
	vector<flat_chain_t *> chainsT;
	read_flat_chains_save_lines(qfn, chainsQ);
	read_flat_chains_save_lines(tfn, chainsT);

	const uint chain_count_q = uint(chainsQ.size());
	const uint chain_count_t = uint(chainsT.size());
	if (chain_count_q == 0)
		Die("No chains found in %s", qfn.c_str());
	if (chain_count_t == 0)
		Die("No chains found in %s", tfn.c_str());

	check_chain_lengths(chainsQ, qfn);
	check_chain_lengths(chainsT, tfn);

	flat_params params;
	params.init_from_cmdline();

	vector<chain_data *> cdvecQ;
	vector<chain_data *> cdvecT;
	build_chain_data_vec(params, chainsQ, cdvecQ);
	build_chain_data_vec(params, chainsT, cdvecT);

	flat_aligner fa;
	fa.m_params = &params;
	fa.alloc();

	float best_score = -FLT_MAX;
	uint best_chain_index_q = UINT_MAX;
	uint best_chain_index_t = UINT_MAX;
	for (uint chain_index_q = 0; chain_index_q < chain_count_q; ++chain_index_q)
		{
		for (uint chain_index_t = 0; chain_index_t < chain_count_t; ++chain_index_t)
			{
			float score = AlignPairFlat(fa,
				*cdvecQ[chain_index_q],
				*cdvecT[chain_index_t],
				false);
			if (score > best_score)
				{
				best_score = score;
				best_chain_index_q = chain_index_q;
				best_chain_index_t = chain_index_t;
				}
			}
		}

	if (best_chain_index_q == UINT_MAX)
		Die("No alignment found");

	AlignPairFlat(fa,
		*cdvecQ[best_chain_index_q],
		*cdvecT[best_chain_index_t],
		true);

	fa.freemem();
	}
