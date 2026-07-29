#include "myutils.h"
#include "flat_params.h"
#include "bcadata.h"
#include "chaq.h"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>

void cmd_getchains()
	{
	if (!optset_label && !optset_labels)
		Die("Must set -label or -labels");

	vector<string> Labels;
	if (optset_labels)
		{
		ReadLinesFromFile(opt(labels), Labels);
		if (Labels.empty())
			Die("No labels found in '%s'", opt(labels));
		}
	if (optset_label)
		Labels.push_back(opt(label));

	unordered_set<string> requested;
	for (uint i = 0; i < SIZE(Labels); ++i)
		{
		ToUpper(Labels[i]);
		requested.insert(Labels[i]);
		}
	const uint n_requested = SIZE(requested);

	const bool want_cal = optset_cal;
	const bool want_bca = optset_bca;
	const bool want_bcb = optset_bcb;
	const bool want_fasta = optset_fasta;
	const bool want_any = want_cal || want_bca || want_bcb || want_fasta;

	BCAData src;
	src.Open(g_Arg1);
	const uint chain_count = src.GetChainCount();

	unordered_map<string, uint> label2idx;
	label2idx.reserve(chain_count);
	for (uint i = 0; i < chain_count; ++i)
		{
		string lab = src.GetLabel(i);
		ToUpper(lab);
		label2idx[lab] = i;
		}

	vector<uint> found_idxs;
	vector<string> missing;
	found_idxs.reserve(n_requested);

	for (auto iter = requested.begin(); iter != requested.end(); ++iter)
		{
		const string &lab = *iter;
		auto bcx_iter = label2idx.find(lab);
		if (bcx_iter == label2idx.end())
			{
			missing.push_back(lab);
			continue;
			}
		found_idxs.push_back(bcx_iter->second);
		}

	sort(found_idxs.begin(), found_idxs.end());

	const uint n_found = SIZE(found_idxs);
	const uint n_missing = SIZE(missing);

	if (want_any)
		{
		FILE *fCal = CreateStdioFile(opt(cal));
		FILE *fFasta = CreateStdioFile(opt(fasta));

		BCAData out_bca;
		BCAData out_bcb;
		BCAData *ptr_bca = 0;
		BCAData *ptr_bcb = 0;
		if (want_bca)
			{
			out_bca.Create(opt(bca), false);
			ptr_bca = &out_bca;
			}
		if (want_bcb)
			{
			out_bcb.Create(opt(bcb), true);
			ptr_bcb = &out_bcb;
			}

		chaq_vecs2 cv;
		bool cv_inited = false;
		if (ptr_bcb != 0)
			{
			chaq::alloc_chaq_vecs2(cv, flat_params::m_maxL);
			cv_inited = true;
			}

		uint8_t *nu_read = 0;
		if (src.m_HasNuSequences)
			nu_read = myalloc(uint8_t, flat_params::m_maxL);

		for (uint i = 0; i < n_found; ++i)
			{
			ProgressStep(i, n_found, "Writing output");
			const uint idx = found_idxs[i];
			flat_chain_t *chain = src.read_flat_chain(idx);
			const uint L = chain->get_length();

			if (nu_read != 0)
				{
				const uint nL = src.read_codeseq_nu(nu_read, idx,
					flat_params::m_maxL);
				asserta(nL == L);
				chain->set_nu_codes(nu_read, L);
				}

			if (fCal != 0)
				chain->to_cal(fCal);
			if (fFasta != 0)
				chain->to_fasta(fFasta);
			if (ptr_bca != 0)
				ptr_bca->write_flat_chain(chain, &cv);
			if (ptr_bcb != 0)
				ptr_bcb->write_flat_chain(chain, &cv);

			delete chain;
			}

		myfree(nu_read);
		if (cv_inited)
			chaq::free_chaq_vecs2(cv);

		if (ptr_bca != 0)
			out_bca.Close();
		if (ptr_bcb != 0)
			out_bcb.Close();
		CloseStdioFile(fCal);
		CloseStdioFile(fFasta);
		}

	src.Close();

	ProgressLog("%u / %u labels found, %u not found\n",
	  n_found, n_requested, n_missing);
	if (n_missing > 0)
		{
		Progress("%u not found", n_missing);
		Log("%u not found\n", n_missing);
		uint Counter = 0;
		for (uint i = 0; i < n_missing; ++i)
			{
			if (Counter++ < 3)
				Progress(" %s", missing[i].c_str());
			Log(">%s\n", missing[i].c_str());
			}
		Progress("\n");
		}
	}
