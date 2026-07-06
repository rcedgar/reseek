#include "myutils.h"
#include "flat_params.h"
#include "bcadata.h"
#include "chaq.h"

void cmd_bcx_subsample()
	{
	if (!optset_n)
		Die("Must specify -n (number of chains)");

	const bool want_cal = optset_cal;
	const bool want_bca = optset_bca;
	const bool want_bcb = optset_bcb;
	const bool want_fasta = optset_fasta;

	if (!want_cal && !want_bca && !want_bcb && !want_fasta)
		Die("Must specify one or more output options: "
		  "-cal, -bca, -bcb, -fasta");

	const uint min_chain_length =
		(optset_minchainlength ? opt(minchainlength) : 80);
	const uint n_want = opt(n);

	BCAData src;
	src.Open(g_Arg1);
	const uint chain_count = src.GetChainCount();

	vector<uint> eligible;
	eligible.reserve(chain_count);
	for (uint i = 0; i < chain_count; ++i)
		{
		if (src.GetSeqLength(i) >= min_chain_length)
			eligible.push_back(i);
		}

	const uint n_eligible = SIZE(eligible);
	if (n_eligible == 0)
		Die("No chains with length >= %u", min_chain_length);

	const uint n_out = min(n_want, n_eligible);
	if (n_want > n_eligible)
		ProgressLog("Warning: only %u eligible chains (requested %u)\n",
		  n_eligible, n_want);

	Shuffle(eligible);
	eligible.resize(n_out);

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

	const uint too_short = chain_count - n_eligible;

	for (uint i = 0; i < n_out; ++i)
		{
		ProgressStep(i, n_out, "Subsample");
		const uint idx = eligible[i];
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

	src.Close();
	if (ptr_bca != 0)
		out_bca.Close();
	if (ptr_bcb != 0)
		out_bcb.Close();

	ProgressLog("%u chains written (%u eligible, %u total, %u too short)\n",
	  n_out, n_eligible, chain_count, too_short);
	}
