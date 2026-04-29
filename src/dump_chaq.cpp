#include "myutils.h"
#include "chaq.h"
#include "flat_chain.h"
#include "flat_distmx.h"
#include "flat_params.h"
#include "alpha.h"
#include "quantize.h"

static void dump_int16(FILE *f, uint16_t i)
	{
	if (i == UINT16_MAX)
		fprintf(f, "\t.");
	else if (i == UINT16_MAX-1)
		fprintf(f, "\t!");
	else
		fprintf(f, "\t%u", i);
	}

static void log_nen(flat_chain_t *chain, uint pos)
	{
	// int DSSParams::m_NEN_w = 12;
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;

	const uint L = chain->get_length();
	sid_t *nensids = myalloc(sid_t, L);
	sid_t *rensids = myalloc(sid_t, L);
	sid_t *pensids = myalloc(sid_t, L);
	sid_t *mensids = myalloc(sid_t, L);
	uint16_t *distmx = myalloc(sid_t, L*M);
	uint16_t *nens = myalloc(uint16_t, L);
	uint16_t *rens = myalloc(uint16_t, L);
	uint16_t *pens = myalloc(uint16_t, L);
	uint16_t *mens = myalloc(uint16_t, L);

	memset(nens, 0xff, L*sizeof(nens[0]));
	memset(rens, 0xff, L*sizeof(nens[0]));
	memset(pens, 0xff, L*sizeof(nens[0]));
	memset(mens, 0xff, L*sizeof(nens[0]));

	memset(nensids, 0xff, L*sizeof(nensids[0]));
	memset(rensids, 0xff, L*sizeof(nensids[0]));
	memset(pensids, 0xff, L*sizeof(nensids[0]));
	memset(mensids, 0xff, L*sizeof(nensids[0]));
	memset(distmx, 0xff, L*M*sizeof(distmx[0]));

	chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);

	chaq::fill_pen_men_vecs(
		distmx, L, pens, pensids, mens, mensids);

	chaq::fill_nen_ren_vecs(pens, mens, pensids, mensids, L,
		nens, rens, nensids, rensids);

	uint nen = nens[pos];
	uint nen2 = UINT_MAX;
	uint minsid = UINT16_MAX;
	for (uint i = 0; i < L; ++i)
		{
		if (abs(int(pos) - int(i)) < int(m))
			{
			Log("[%3u]  (too close)\n", i);
			continue;
			}
		uint k = banded_ij_to_k(i, pos);
		uint sid = distmx[k];
		float d = sid2dist(sid);
		Log("[%3u]  %5u  %7.1f\n", i, sid, d);
		if (sid < minsid)
			{
			nen2 = i;
			minsid = sid;
			}
		}
	Log("\nnen=%u, nen2=%u, minsid=%u, mind=%.1f\n",
		nen, nen2, minsid, sid2dist(minsid));
	}

static void dump_chaq(FILE *f, flat_chain_t *chain)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint m = flat_params::m_nn_min_offset;
	const uint L = chain->get_length();
	sid_t *nensids = myalloc(sid_t, L);
	sid_t *rensids = myalloc(sid_t, L);
	sid_t *pensids = myalloc(sid_t, L);
	sid_t *mensids = myalloc(sid_t, L);
	uint16_t *distmx = myalloc(sid_t, L*M);
	uint16_t *nens = myalloc(uint16_t, L);
	uint16_t *rens = myalloc(uint16_t, L);
	uint16_t *pens = myalloc(uint16_t, L);
	uint16_t *mens = myalloc(uint16_t, L);

	memset(nens, 0xff, L*sizeof(nens[0]));
	memset(rens, 0xff, L*sizeof(nens[0]));
	memset(pens, 0xff, L*sizeof(nens[0]));
	memset(mens, 0xff, L*sizeof(nens[0]));

	memset(nensids, 0xff, L*sizeof(nensids[0]));
	memset(rensids, 0xff, L*sizeof(nensids[0]));
	memset(pensids, 0xff, L*sizeof(nensids[0]));
	memset(mensids, 0xff, L*sizeof(nensids[0]));
	memset(distmx, 0xff, L*M*sizeof(distmx[0]));

	chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);

	chaq::fill_pen_men_vecs(
		distmx, L, pens, pensids, mens, mensids);

	chaq::fill_nen_ren_vecs(pens, mens, pensids, mensids, L,
		nens, rens, nensids, rensids);

	fprintf(f, "pos");
	fprintf(f, "\taa");

	fprintf(f, "\tnen");
	fprintf(f, "\tren");
	fprintf(f, "\tpen");
	fprintf(f, "\tmen");

	fprintf(f, "\tnensid");
	fprintf(f, "\trensid");
	fprintf(f, "\tpensid");
	fprintf(f, "\tmensid");

	fprintf(f, "\tnendist");
	fprintf(f, "\trendist");
	fprintf(f, "\tpendist");
	fprintf(f, "\tmendist");
	fprintf(f, "\n");

	for (uint pos = 0; pos < L; ++pos)
		{
		fprintf(f, "%u", pos);
		fprintf(f, "\t%c", chain->m_aa->m_data[pos]);

		uint16_t nen = nens[pos];
		uint16_t ren = rens[pos];
		uint16_t pen = pens[pos];
		uint16_t men = mens[pos];

		uint16_t nensid = nensids[pos];
		uint16_t rensid = rensids[pos];
		uint16_t pensid = pensids[pos];
		uint16_t mensid = mensids[pos];

		float nendist = (nensid == UINT16_MAX ? -1 : sid2dist(nensid));
		float rendist = (rensid == UINT16_MAX ? -1 : sid2dist(rensid));
		float pendist = (pensid == UINT16_MAX ? -1 : sid2dist(pensid));
		float mendist = (mensid == UINT16_MAX ? -1 : sid2dist(mensid));

		dump_int16(f, nen);
		dump_int16(f, ren);
		dump_int16(f, pen);
		dump_int16(f, men);

		dump_int16(f, nensid);
		dump_int16(f, rensid);
		dump_int16(f, pensid);
		dump_int16(f, mensid);

		fprintf(f, "\t%.3g", nendist);
		fprintf(f, "\t%.3g", rendist);
		fprintf(f, "\t%.3g", pendist);
		fprintf(f, "\t%.3g", mendist);

		fprintf(f, "\n");
		}
	}

void cmd_dump_chaq()
	{
	const string &chainfn = g_Arg1;
	vector<vector<uint8_t> > codeseqs;
	if (optset_bandwidth) flat_params::m_distmx_bandwidth = opt(bandwidth);
	FILE *f = CreateStdioFile(opt(output));

	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	uint nchain = SIZE(chains);

	vector<uint16_t> counts(UINT16_MAX+1);
	uint nbad = 0;
	for (uint i = 0; i < nchain; ++i)
		dump_chaq(f, chains[i]);

	CloseStdioFile(f);
	}
