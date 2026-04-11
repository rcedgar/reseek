#include "myutils.h"
#include "flat_chain.h"
#include "cigar.h"
#include "chaq.h"

void trunc_label(string &Label);

void read_flat_chains_idx_trunclabel(
	const string &fn,
	vector<flat_chain_t *> &chains,
	unordered_map<string, uint> &label2idx);

float flat_getlddt_muscle_some_floats(
	const uint32_t *posQs,
	const uint32_t LQ,
	const uint32_t *posTs,
	const uint32_t LT,
	const uint ncol,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	const uint M,
	uint32_t *nr_considered_vec,
	uint32_t *nr_preserved_vec);

void path2posvecs(
	const string &path,
	uint loQ, uint LQ,
	uint loT, uint LT,
	vector<uint> &posQs,
	vector<uint> &posTs)
	{
	posQs.clear();
	posTs.clear();
	const uint colcount = uint(path.size());
	posQs.reserve(colcount);
	posTs.reserve(colcount);
	uint posQ = loQ;
	uint posT = loT;
	for (uint col = 0; col < colcount; ++col)
		{
		char c = path[col];
		if (c == 'M')
			{
			assert(posQ < LQ);
			assert(posT < LT);
			posQs.push_back(posQ);
			posTs.push_back(posT);
			}
		if (c == 'M' || c == 'D')
			posQ++;
		if (c == 'M' || c == 'I')
			posT++;
		}
	}

static const uint M = 64;

static void structure_features(
	FILE *f,
	flat_chain_t *chainQ,
	flat_chain_t *chainT,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	uint loQ, uint LQ, uint loT, uint LT,
	const string &CIGAR)
	{
	asserta(f != 0);
	string path;
	CIGARToPath(CIGAR, path, true);

	const string labelQ = chainQ->m_label;
	const string labelT = chainT->m_label;
	const uint LQ2 = chainQ->get_length();
	const uint LT2 = chainT->get_length();
	asserta(LQ2 == LQ);
	asserta(LT2 == LT);

	vector<uint> posQs, posTs;
	path2posvecs(path, loQ, LQ, loT, LT, posQs, posTs);
	const uint ncol = uint(posQs.size());
	assert(posTs.size() == ncol);

	float lddt = 0;
	if (ncol > 8)
		{
		uint32_t *nr_considered_vec = myalloc(uint32_t, ncol);
		uint32_t *nr_preserved_vec = myalloc(uint32_t, ncol);
		lddt = flat_getlddt_muscle_some_floats(
			posQs.data(), LQ,
			posTs.data(), LT,
			ncol,
			distmxQ,
			distmxT,
			M,
			nr_considered_vec,
			nr_preserved_vec);
		myfree(nr_considered_vec);
		myfree(nr_preserved_vec);
		}
	fprintf(f, "%s", labelQ.c_str());
	fprintf(f, "\t%s", labelT.c_str());
	fprintf(f, "\t%.4g", lddt);
	fprintf(f, "\n");
	}

void cmd_structure_features()
	{
	asserta(optset_input); // tsv with q,t,score,CIGAR
	asserta(optset_output);
	vector<flat_chain_t *> chains;
	unordered_map<string, uint> label2chainidx;
	read_flat_chains_idx_trunclabel(g_Arg1, chains, label2chainidx);
	const uint nchain = uint(chains.size());

	vector<sid_t *> distmxs;
	vector<uint32_t> Ls;
	distmxs.reserve(nchain);
	Ls.reserve(nchain);
	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		ProgressStep(chainidx, nchain, "Distance mxs");
		const flat_chain_t *chain = chains[chainidx];
		const uint L = chain->get_length();
		Ls.push_back(L);
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain->m_xyz->m_data, L, M, distmx);
		distmxs.push_back(distmx);
		}

	FILE *fout = CreateStdioFile(opt(output));

	string line;
	vector<string> flds;
	FILE *fin = OpenStdioFile(opt(input));
	uint64_t file_size = GetStdioFileSize64(fin);
	uint64_t last_file_pos = 0;
	Progress("Reading hits 0%%\r");
	for (;;)
		{
		bool ok = ReadLineStdioFile(fin, line);
		if (!ok)
			break;
		uint64_t file_pos = GetStdioFilePos64(fin);
		if (file_pos - last_file_pos > 10e6)
			{
			double pct = file_pos*100.0/(file_size+1);
			Progress("Reading hits %.2f%%\r", pct);
			last_file_pos = file_pos;
			}
	
	//       0      1      2     3       4     5         6     7
	//       Q      T  score   loQ      LQ   loT        LT     CIGAR
	// d3mkbb_ d2gdma_ 44.94   1       133     1       153     14M2D24M2D10M9D12M3D7M2D14M1I49M
		Split(line, flds, '\t');
		asserta(flds.size() == 8);
		string labelQ = flds[0];
		string labelT = flds[1];
		uint loQ = StrToUint(flds[3]);
		uint LQ = StrToUint(flds[4]);
		uint loT = StrToUint(flds[5]);
		uint LT = StrToUint(flds[6]);
		const string &CIGAR = flds[7];

		trunc_label(labelQ);
		trunc_label(labelT);

		unordered_map<string, uint>::iterator iterQ =
			label2chainidx.find(labelQ);
		unordered_map<string, uint>::iterator iterT =
			label2chainidx.find(labelT);

		if (iterQ == label2chainidx.end()) Die("Chain not found >%s", labelQ.c_str());
		if (iterT == label2chainidx.end()) Die("Chain not found >%s", labelT.c_str());

		uint chainidxQ = iterQ->second;
		uint chainidxT = iterT->second;

		flat_chain_t *chainQ = chains[chainidxQ];
		flat_chain_t *chainT = chains[chainidxT];

		const uint LQ2 = chainQ->get_length();
		const uint LT2 = chainT->get_length();

		asserta(LQ2 == Ls[chainidxQ]);
		asserta(LQ2 == LQ);

		asserta(LT == Ls[chainidxT]);
		asserta(LT2 == LT);

		const sid_t *distmxQ = distmxs[chainidxQ];
		const sid_t *distmxT = distmxs[chainidxT];

		structure_features(fout, chainQ, chainT, distmxQ, distmxT,
			loQ, LQ, loT, LT, CIGAR);
		}
	Progress("Reading hits 100%%   \n");
	CloseStdioFile(fout);
	}