#pragma once

#include "flat_bench.h"

class flat_aligner;

class flat_bench_struct_feature : public flat_bench
	{
public:
	vector<flat_chain_t *> m_chains;
	vector<sid_t *> m_distmxs;

public:
	void read_chains(const string &fn);
	void set_distmxs(uint M);

	float get_feature_value(uint idxQ, uint idxT,
		const flat_aligner &fa) const;

	float get_lddt(uint idxQ, uint idxT,
		const flat_aligner &fa) const;

	float get_dali(uint idxQ, uint idxT,
		const flat_aligner &fa) const;

	float get_dalix(uint idxQ, uint idxT,
		const flat_aligner &fa) const;

	float get_entropy(uint idxQ, uint idxT,
		const flat_aligner &fa) const;

public:
	virtual void ThreadBody_All(uint ThreadIdx);
	virtual void ThreadBody_Dope(uint ThreadIdx);
	};
