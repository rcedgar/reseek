#include "myutils.h"
#include "seqdb.h"
#include "flat_bench_kmer.h"

#define	SHOW_PROGRESS	1

void flat_bench_kmer::align_pair_kmer(uint domidxq, uint domidxt)
	{
	assert(domidxq < m_kmerseqvec.size());
	assert(domidxt < m_kmerseqvec.size());

	const vector<uint32_t> &kmersq = m_kmerseqvec[domidxq];
	const vector<uint32_t> &kmerst = m_kmerseqvec[domidxt];

	set<uint32_t> kmersetq;
	for (auto kmer : kmersq)
		kmersetq.insert(kmer);
	uint n = 0;
	for (auto kmer: kmerst)
		{
		if (kmersetq.find(kmer) != kmersetq.end())
			++n;
		}

	AppendHit(domidxq, domidxt, float(n));
	}

void flat_bench_kmer::read_kmers(const string &fastafn, 
	uint alpha_size, uint k)
	{
	asserta(alpha_size <= 36);
	asserta(alpha_size != 20);
	asserta(k >= 2 && k <= 6);

	SeqDB DB;
	DB.FromFasta(fastafn);
	const uint nseq = DB.GetSeqCount();
	DB.ToLetters(g_CharToLetterMu);
	m_kmerseqvec.clear();
	m_kmerseqvec.resize(nseq);
	for (uint seqidx = 0; seqidx < nseq; ++seqidx)
		{
		const string &label = DB.GetLabel(seqidx);
		const uint L = DB.GetSeqLength(seqidx);
		const uint8_t *byteseq = DB.GetByteSeq(seqidx);
		if (L < 2*k)
			continue;
		const uint nk = L - k + 1;

		uint domidx = m_look->get_domidx(label);
		vector<uint32_t> &kmers = m_kmerseqvec[domidx];
		assert(domidx < m_kmerseqvec.size());

		kmers.reserve(nk);
		for (uint pos = 0; pos < nk; ++pos)
			{
			uint32_t kmer = 0;
			for (uint j = 0; j < k; ++j)
				{
				uint8_t code = byteseq[pos+j];
				asserta(code < alpha_size);
				kmer = kmer*alpha_size + code;
				}
			kmers.push_back(kmer);
			}
		}
	}

void flat_bench_kmer::ThreadBody_All(uint ThreadIdx)
	{
	assert(m_look);
	const uint ndom = m_look->get_ndom();
	const uint npair = m_look->get_pair_count_upper_triangle_with_diagonal();
	uint CurrentDomIdxT = UINT_MAX;
	for (;;)
		{
		uint DomIdxT = m_NextQueryIdx++;
		if (DomIdxT >= ndom)
			return;

		for (uint DomIdxQ = DomIdxT+1; DomIdxQ < ndom; ++DomIdxQ)
			{
			uint PairIdx = triangle_ij_to_k(DomIdxT, DomIdxQ, ndom);
			uint progress_count = m_progress_counter++;
#if SHOW_PROGRESS
			if (progress_count%1000 == 0)
				ProgressStep(progress_count, npair, "Aligning");
#endif
			align_pair_kmer(DomIdxT, DomIdxQ);
			}
		}
	}

void cmd_flat_bench_kmer()
	{
	asserta(optset_lookup);
	asserta(optset_alpha_size);
	asserta(optset_k);

	asserta(!optset_dope);
	asserta(!optset_fapattern);
	asserta(!optset_mxpattern);
	asserta(!optset_spec);
	asserta(!optset_varstr);

	const string &fastafn = g_Arg1;

	flat_bench_kmer FB;
	FB.ReadLookup(opt(lookup));
	FB.read_kmers(fastafn, opt(alpha_size), opt(k));

	FB.Alloc();
	FB.Search("all");
	FB.SetScoreOrder();
	FB.Bench();
	FB.WriteHits(opt(output), true);
	}
