#include "myutils.h"
#include "landmark.h"
#include "flat_chain.h"
#include "chaq.h"
#include "seqdb.h"

static const uint M = 32;
static const uint MINL = 80;
static const uint MAXL = 1000;

static void load_msas(
	const string &msafilesfn,
	vector<SeqDB *> &MSAs,
	vector<string> &msastemnames)
	{
	MSAs.clear();
	msastemnames.clear();
	vector<string> msafns;
	ReadLinesFromFile(msafilesfn, msafns);
	const uint n = uint(msafns.size());
	MSAs.reserve(n);
	for (uint i = 0; i < n; ++i)
		{
		ProgressStep(i, n, "Loading MSAs");
		SeqDB *MSA = new SeqDB;
		MSA->FromFasta(msafns[i], true);
		asserta(MSA->IsAligned());
		MSAs.push_back(MSA);

		vector<string> flds;
		Split(msafns[i], flds, '/');
		msastemnames.push_back(flds[flds.size()-1]);
		}
	}

void cmd_landmark()
	{
	asserta(optset_input);
	asserta(optset_output);
	asserta(optset_output2);
	const string &chainsfn = g_Arg1;
	const string &msafilesfn = opt(input);
	FILE *ffa = CreateStdioFile(opt(fasta));

	vector<SeqDB *> MSAs;
	vector<string> msastemnames;
	load_msas(msafilesfn, MSAs, msastemnames);
	const uint nmsa = uint(MSAs.size());
	ProgressLog("%u msas\n", nmsa);

	vector<flat_chain_t *> chains;
	read_flat_chains(chainsfn, chains);
	const uint nchain = uint(chains.size());
	uint16_t *distmx = myalloc(uint16_t, MAXL*M);

	unordered_map<string, uint> label2chainidx;
	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		const string &label = chains[chainidx]->m_label;
		label2chainidx[label] = chainidx;
		}

	uint used = 0;
	uint notfound = 0;
	uint tooshort = 0;
	vector<uint> pos2col;
	uint total_residues = 0;
	uint total_landmarks = 0;
	uint total_overlaps = 0;
	uint n_coil = 0;
	uint n_helix = 0;
	uint n_strand = 0;
	for (uint msaidx = 0; msaidx < nmsa; ++msaidx)
		{
		ProgressStep(msaidx, nmsa, "Landmark");

		const string &msastemname = msastemnames[msaidx];
		string output_msafn = opt(output2) + msastemname;
		FILE *foutmsa = CreateStdioFile(output_msafn);

		SeqDB &MSA = *MSAs[msaidx];
		const uint nrow = MSA.GetSeqCount();
		const uint ncol = MSA.GetColCount();

		for (uint rowidx = 0; rowidx < nrow; ++rowidx)
			{
			const string &label = MSA.GetLabel(rowidx);
			const string &row = MSA.GetSeq(rowidx);
			unordered_map<string, uint>::const_iterator iter =
				label2chainidx.find(label);
			if (iter == label2chainidx.end())
				{
				++notfound;
				continue;
				}

			const uint chainidx = iter->second;

			flat_chain_t *chain = chains[chainidx];
			const uint L = chain->get_length();
			if (L < MINL || L > MAXL)
				{
				++tooshort;
				continue;
				}
			++used;

			total_residues += L;
			const ic_t *xyz = chain->m_xyz->m_data;
			chaq::fill_distmx(xyz, L, M, distmx);

			string ss_str;
			chaq::get_ss4_str(distmx, M, L, ss_str);
			asserta(ss_str.size() == L);

			uint pos = 0;
			string landmark_row;
			landmark_row.resize(ncol, '!');
			for (uint colidx = 0; colidx < ncol; ++colidx)
				{
				char c = row[colidx];
				if (isgap(c))
					landmark_row[colidx] = '-';
				else
					{
					landmark_row[colidx] = 'A';
					pos2col.push_back(colidx);
					}
				}

			string landmark_seq;
			landmark_seq.resize(L, 'A');
			SeqToFasta(ffa, label, landmark_seq);
			SeqToFasta(foutmsa, label, landmark_row);
			}
		CloseStdioFile(foutmsa);
		}
	ProgressLog("total_residues %u, landmarks %u, overlaps %u\n",
		total_residues, total_landmarks, total_overlaps);
	ProgressLog("helix %u, strand %u, coil %u\n",
		n_helix, n_strand, n_coil);
	ProgressLog("%u notound, %u tooshort, %u used\n",
		notfound, tooshort, used);

	CloseStdioFile(ffa);
	}
