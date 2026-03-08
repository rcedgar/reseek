#include "myutils.h"
#include "flat_dist_types.h"
#include "flat_base.h"
#include "flat_chain.h"
#include "flat_distmx.h"
#include "chaq.h"

static uint M = 64;
static uint m = 8;

static void upd_nen_di(const chaindistmx_t *dm,
	uint L, uint i, vector<uint> &counts)
	{
	const uint16_t *distmx = dm->m_data;
	uint n = SIZE(counts);
	uint nendi = 0;
	sid_t minsid = UINT16_MAX;
	for (uint j = 0; j < L; ++j)
		{
		uint di = uint(abs(int(i) - int(j)));
		if (di >= M)
			continue;
		if (di < int(m))
			continue;
		uint k = banded_ij_to_k(M, i, j);
		sid_t sid = distmx[k];
		if (sid < minsid)
			{
			nendi = di;
			minsid = sid;
			}
		}
	if (nendi < n)
		++counts[nendi];
	}

static void upd_fen_di(const chaindistmx_t *dm,
	uint L, uint i, vector<uint> &counts)
	{
	const uint16_t *distmx = dm->m_data;
	uint n = SIZE(counts);
	uint fendi = 0;
	sid_t maxsid = 0;
	for (uint j = 0; j < L; ++j)
		{
		uint di = uint(abs(int(i) - int(j)));
		if (di >= M)
			continue;
		if (di < int(m))
			continue;
		uint k = banded_ij_to_k(M, i, j);
		sid_t sid = distmx[k];
		if (sid > maxsid)
			{
			fendi = di;
			maxsid = sid;
			}
		}
	if (fendi < n)
		++counts[fendi];
	}

static void upd_nen_dx(const chaindistmx_t *dm,
	uint L, uint i, vector<uint> &counts)
	{
	const uint16_t *distmx = dm->m_data;
	uint n = SIZE(counts);
	uint nendi = 0;
	sid_t minsid = UINT16_MAX;
	for (uint j = 0; j < L; ++j)
		{
		uint di = uint(abs(int(i) - int(j)));
		if (di >= M)
			continue;
		if (di < int(m))
			continue;
		uint k = banded_ij_to_k(M, i, j);
		sid_t sid = distmx[k];
		if (sid < minsid)
			{
			nendi = di;
			minsid = sid;
			}
		}
	uint dist_Angstroms = sid2ic[minsid]/10;
	if (dist_Angstroms < n)
		++counts[dist_Angstroms];
	}

static void upd_fen_dx(const chaindistmx_t *dm,
	uint L, uint i, vector<uint> &counts)
	{
	const uint16_t *distmx = dm->m_data;
	uint n = SIZE(counts);
	uint fendi = 0;
	sid_t maxsid = 0;
	for (uint j = 0; j < L; ++j)
		{
		uint di = uint(abs(int(i) - int(j)));
		if (di >= M)
			continue;
		if (di < int(m))
			continue;
		uint k = banded_ij_to_k(M, i, j);
		sid_t sid = distmx[k];
		if (sid > maxsid)
			{
			fendi = di;
			maxsid = sid;
			}
		}
	uint dist_Angstroms = sid2ic[maxsid]/10;
	if (dist_Angstroms < n)
		++counts[dist_Angstroms];
	}

static void write_dist(FILE *f,
	const vector<flat_chain *> &chains,
	const vector<const chaindistmx_t *> &dms,
	const string &name, uint N)
	{
	if (f == 0)
		return;
	vector<uint> counts(N);

	Progress("%s ...", name.c_str());
	const uint nchains = SIZE(dms);
	for (uint chidx = 0; chidx < nchains; ++chidx)
		{
		const flat_chain *chain = chains[chidx];
		const uint L = chain->get_length();
		const chaindistmx_t *dm = dms[chidx];
		for (uint i = 0; i < L; ++i)
			{
			if (name == "nen_di")
				upd_nen_di(dm, L, i, counts);
			else if (name == "nen_dx")
				upd_nen_dx(dm, L, i, counts);
			if (name == "fen_di")
				upd_fen_di(dm, L, i, counts);
			else if (name == "fen_dx")
				upd_fen_dx(dm, L, i, counts);
			}
		}
	Progress(" done\n");

	fprintf(f, "counts\t%s(%u)\t%u\n", name.c_str(), M, N);
	for (uint i = 0; i < N; ++i)
		fprintf(f, "%u\t%i\n", i, counts[i]);
	}

static void get_dms(vector<flat_chain *> &chains,
	vector<const chaindistmx_t *> &dms)
	{
	dms.clear();

	chaq c;
	const uint nchains = SIZE(chains);
	dms.reserve(nchains);
	for (uint chidx = 0; chidx < nchains; ++chidx)
		{
		ProgressStep(chidx, nchains, "Reading chains");
		const flat_chain *chain = chains[chidx];

		const uint L = chain->get_length();
		if (L < 8)
			continue;

		c.init(chain);
		const chaindistmx_t *dm = c.get_distmx(M);
		dms.push_back(dm);
		}
	}

void cmd_nbr_stats()
	{
	vector<flat_chain *> chains;
	read_flat_chains(g_Arg1, chains);
	FILE *fOut = CreateStdioFile(opt(output));

	M = 64;
	vector<const chaindistmx_t *> dms;
	get_dms(chains, dms);

	write_dist(fOut, chains, dms, "nen_di", M);
	write_dist(fOut, chains, dms, "nen_dx", 75);

	write_dist(fOut, chains, dms, "fen_di", M);
	write_dist(fOut, chains, dms, "fen_dx", 100);

	M = 128;
	get_dms(chains, dms);
	write_dist(fOut, chains, dms, "fen_di", M);
	write_dist(fOut, chains, dms, "fen_dx", 100);

	M = 256;
	get_dms(chains, dms);
	write_dist(fOut, chains, dms, "fen_di", M);
	write_dist(fOut, chains, dms, "fen_dx", 100);

	M = 512;
	get_dms(chains, dms);
	write_dist(fOut, chains, dms, "fen_di", M);
	write_dist(fOut, chains, dms, "fen_dx", 100);

	CloseStdioFile(fOut);

	log_flat_stats();
	}
