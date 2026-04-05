#pragma once

#include "flat_bench.h"

class flat_bench_kmer : public flat_bench
	{
public:
	vector<vector<uint32_t> > m_kmerseqvec;

public:
	void read_kmers(const string &fastafn,
		uint alpha_size, uint k);
	void align_pair_kmer(uint DomIdxT, uint DomIdxQ);

public:
	virtual void ThreadBody_All(uint ThreadIdx);
	virtual void ThreadBody_Dope(uint ThreadIdx) { Die("dope"); }
	};
