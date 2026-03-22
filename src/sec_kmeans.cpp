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
		int dij = max(off1, off2) - min(off1, off2);
		asserta(dij > 1);
		}
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
	if (optset_spec)
		{
		vector<string> flds;
		Split(opt(spec), flds, ',');
		size_t n = flds.size();
		asserta(n%2 == 0);
		n /= 2;
		for (uint i = 0; i < n; ++i)
			{
			off1s.push_back(StrToInt(flds[2*i]));
			off2s.push_back(StrToInt(flds[2*i+1]));
			}
		}
	else
		{
		off1s = default_off1s;
		off2s = default_off2s;
		}
	validate_offs(off1s, off2s);

	const uint K = opt(alpha_size);
	const uint32_t M = 32; // dist mx band width

	sec_kmeans SK;
	SK.init(K, M, off1s, off2s);

	sid_t *vs = myalloc(sid_t, SK.m_N*SK.m_D);
	SK.m_cluster_idxs = myalloc(uint, SK.m_N);
	SK.set_vs(chains);
	SK.train();
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
