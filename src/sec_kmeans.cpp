#include "myutils.h"
#include "flat_chain.h"
#include "sec_kmeans.h"
#include "chaq.h"
#include "quarts.h"
#include "alpha.h"

static bool check_backbone(const flat_chain_t* chain)
	{
	const uint L = chain->get_length();
	for (uint i = 1; i < L; ++i)
		{
		float d = chain->slow_float_dist(i-1, i);
		if (d < 3.7 || d > 3.9)
			return false;
		}
	return true;
	}

void cmd_sec_variance()
	{
	asserta(optset_output);
	FILE *fOut = CreateStdioFile(opt(output));
	vector<flat_chain_t *>chains;
	read_flat_chains(g_Arg1, chains);
	const uint nrchains = SIZE(chains);
	const uint M = 32;
	const int w = 16;
	vector<vector<float> > dists(w+1);
	for (uint chain_idx = 0; chain_idx < nrchains; ++chain_idx)
		{
		const flat_chain_t* chain = chains[chain_idx];
		if (!check_backbone(chain))
			continue;
		const uint L = chain->get_length();
		for (uint pos = 0; pos < L; ++pos)
		for (uint dij = 1; dij < w; ++dij)
			{
			uint pos2 = pos + dij;
			if (pos2 >= L)
				break;
			float d = chain->slow_float_dist(pos, pos2);
			dists[dij].push_back(d);
			}
		}

	fprintf(fOut, "dij");
	fprintf(fOut, "\tN");
	fprintf(fOut, "\tMin");
	fprintf(fOut, "\tLoQ");
	fprintf(fOut, "\tMed");
	fprintf(fOut, "\tHiQ");
	fprintf(fOut, "\tMax");
	fprintf(fOut, "\tAvg");
	fprintf(fOut, "\tStdDev");
	fprintf(fOut, "\n");
	for (uint dij = 1; dij < w; ++dij)
		{
		const vector<float> &ds = dists[dij];
		QuartsFloat QF;
		GetQuartsFloat(ds, QF);
		ProgressLog("dij=%2u  ", dij);
		QF.ProgressLogMe();

		fprintf(fOut, "%u\t", dij);
		QF.ToTsv(fOut);
		}
	CloseStdioFile(fOut);
	}

void intseq2fasta(FILE *f, const string &label, const uint8_t *intseq, uint L)
	{
	if (f == 0)
		return;
	fprintf(f, ">%s\n", label.c_str());
	for (uint i = 0; i < L; ++i)
		{
		if (i > 0 && i%80 == 0)
			fputc('\n', f);
		char c = g_LetterToCharMu[intseq[i]];
		fputc(c, f);
		}
	fputc('\n', f);
	}

void cmd_sec_fasta()
	{
	const uint M = 32;

	vector<flat_chain_t *>chains;
	read_flat_chains(g_Arg1, chains);
	const uint nrchains = SIZE(chains);

	sec_kmeans SK;
	SK.from_tsv(opt(input));
	SK.m_M = M;

	FILE *ffa = CreateStdioFile(opt(output));

	for (uint chainidx = 0; chainidx < nrchains; ++chainidx)
		{
		const flat_chain_t* chain = chains[chainidx];
		const uint L = chain->get_length();
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain->m_xyz->m_data, L, M, distmx);
		uint8_t *intseq = myalloc(uint8_t, L);
		SK.get_intseq(distmx, L, intseq);
		intseq2fasta(ffa, chain->m_label, intseq, L);
		}
	CloseStdioFile(ffa);
	}

static void validate_offs(
	const vector<int> &off1s,
	const vector<int> &off2s)
	{
	size_t n = off1s.size();
	asserta(off2s.size() == n);
	for (uint i = 0; i < n; ++i)
		{
		int off1 = off1s[i];
		int off2 = off2s[i];
		asserta(off1 < off2);
		int dij = max(off1, off2) - min(off1, off2);
		asserta(dij > 1);
		for (uint j = 0; j < i; ++i)
			{
			if (off1s[j] == off1 && off2s[j] == off2)
				Die("dupe %d,%d", off1, off2);
			}
		}
	}

static void get_rand_off12(uint M, int &off1, int &off2)
	{
	for (uint iter = 0; iter < 100; ++iter)
		{
		uint u = randu32()%(M+1);
		bool plus = (randu32()%2 == 0);
		off1 = (plus ? u : -int(u));

		u = randu32()%(M+1);
		plus = (randu32()%2 == 0);
		off2 = (plus ? u : -int(u));
		int dij = max(off1, off2) - min(off1, off2);
		if (dij > 1)
			return;
		}
	asserta(false);
	}

static bool ok_to_append(
	const vector<int> &off1s,
	const vector<int> &off2s,
	int off1,
	int off2)
	{
	size_t n = off1s.size();
	asserta(off2s.size() == n);
	for (size_t i = 0; i < n; ++i)
		{
		if (off1s[i] == off1 && off2s[i] == off2)
			return false;
		if (off2s[i] == off1 && off1s[i] == off2)
			return false;
		}
	return true;
	}

// r_D_M where D=dimension M=max offset
static void get_random_offs(const string &spec,
	vector<int> &off1s, vector<int> &off2s)
	{
	vector<string> flds;
	Split(spec, flds, '_');
	asserta(flds.size() == 3);
	assert(flds[0] == "r");
	uint D = StrToUint(flds[1]);
	uint M = StrToUint(flds[2]);
	asserta(D > 1 && D < 32);
	asserta(M > 2 && M < 32);
	off1s.clear();
	off2s.clear();

	for (uint d = 0; d < D; ++d)
		{
		bool ok = false;
		for (uint iter = 0; iter < 100; ++iter)
			{
			int off1, off2;
			get_rand_off12(M, off1, off2);
			if (ok_to_append(off1s, off2s, off1, off2))
				{
				off1s.push_back(min(off1,off2));
				off2s.push_back(max(off1,off2));
				ok = true;
				break;
				}
			}
		asserta(ok);
		}
	asserta(off1s.size() == D);
	asserta(off2s.size() == D);
	}

void cmd_sec_kmeans()
	{
	asserta(optset_alpha_size);
	vector<flat_chain_t *>chains;
	read_flat_chains(g_Arg1, chains);
	const uint nrchains = SIZE(chains);

/***
     _________________________________________
	|     |     |     | [7] |     | o-o |  o  |
	|-3,+3|-2,+3|-1,+3| 0,+3|+1,+3|+2,+3|+3,+3|    o = zero
	|_____|_____|_____|_____|_____|_____|_____|    o-o = 3.81 A (adjacent)
	|     | [2] | [4] | [5] | o-o |  o  | o-o |
	|-3,+2|-2,+2|-1,+2| 0,+2|+1,+2|+2,+2|+3,+2|
	|_____|_____|_____|_____|_____|_____|_____|
	|     | [1] | [3] | o-o |  o  | o-o |     |
	|-3,+1|-2,+1|-1,+1| 0,+1|+1,+1|+2,+1|+3,+1|
	|_____|_____|_____|_____|_____|_____|_____|
	| [8] | [0] | o-o |  o  | o-o |     |     |
	|-3,0 |-2,0 |-1,0 | 0,0 |+1,0 |+2,0 |+3,0 |
	|_____|_____|_____|_____|_____|_____|_____|
	|     | o-o |  o  | o-o |     |     |     |
	|-3,-1|-2,-1|-1,-1| 0,-1|+1,-1|+2,-1|+3,-1|
	|_____|_____|_____|_____|_____|_____|_____|
	| o-o |  o  | o-o |     |     |     |     |
	|-3,-2|-2,-2|-1,-2| 0,-2|+1,-2|+2,-2|+3,-2|
	|_____|_____|_____|_____|_____|_____|_____|
	|  o  | o-o |     |     |     |     | [6] |
	|-3,-3|-2,-3|-1,-3| 0,-3|+1,-3|+2,-3|+3,-3|
	|_____|_____|_____|_____|_____|_____|_____|

***/
	//                                   0   1   2   3   4   5   6   7   8
	const vector<int> default_off1s = { -2, -2, -2, -1, -1,  0, -3,  0, -3 };
	const vector<int> default_off2s = {  0,  1,  2,  1,  2,  2,  3,  3,  0 };
	//                                   3   3   5   3   3   3   7   4   4  dij

	vector<int> off1s;
	vector<int> off2s;
	uint D = UINT_MAX;
	if (optset_spec)
		{
		const string &spec = opt(spec);
		if (spec[0] == 'r')
			get_random_offs(spec, off1s, off2s);
		else
			{
			vector<string> flds;
			Split(spec, flds, ',');
			size_t n = flds.size();
			asserta(n%2 == 0);
			n /= 2;
			for (uint i = 0; i < n; ++i)
				{
				int off1 = StrToInt(flds[2*i]);
				int off2 = StrToInt(flds[2*i+1]);
				int minoff = min(off1,off2);
				int maxoff = min(off1,off2);
				off1s.push_back(minoff);
				off2s.push_back(maxoff);
				}
			}
		}
	else
		{
		off1s = default_off1s;
		off2s = default_off2s;
		}
	validate_offs(off1s, off2s);

	Log("off1s "); for (uint i = 0; i < off1s.size(); ++i) Log(" %3d", off1s[i]); Log("\n");
	Log("off2s "); for (uint i = 0; i < off2s.size(); ++i) Log(" %3d", off2s[i]); Log("\n");
	Log("-spec ");
	for (uint i = 0; i < off1s.size(); ++i)
		{
		if (i > 0) Log(",");
		Log("%d,%d", off1s[i], off2s[i]);
		}
	Log("\n");
	
	const uint K = opt(alpha_size);
	const uint32_t M = 32; // dist mx band width

	uint niter = 1000;
	if (optset_iters) niter = opt(iters);

	sec_kmeans SK;
	SK.init(K, M, off1s, off2s);

	sid_t *vs = myalloc(sid_t, SK.m_N*SK.m_D);
	SK.m_cluster_idxs = myalloc(uint, SK.m_N);
	SK.set_vs(chains);
	SK.train(niter);
	SK.logme();
	SK.ss4stats();
	SK.to_tsv(opt(output));
	if (optset_fasta)
		{
		FILE *ffa = CreateStdioFile(opt(fasta));

		for (uint chainidx = 0; chainidx < nrchains; ++chainidx)
			{
			const flat_chain_t* chain = chains[chainidx];
			const uint L = chain->get_length();
			sid_t *distmx = myalloc(sid_t, L*M);
			chaq::fill_distmx(chain->m_xyz->m_data, L, M, distmx);
			uint8_t *intseq = myalloc(uint8_t, L);
			SK.get_intseq(distmx, L, intseq);
			intseq2fasta(ffa, chain->m_label, intseq, L);
			}
		CloseStdioFile(ffa);
		}

#if DEBUG
	sec_kmeans SK2;
	SK2.from_tsv(opt(output));
	SK2.logme();
#endif
	log_flat_stats();
	}
