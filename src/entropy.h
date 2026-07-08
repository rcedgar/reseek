#pragma once

#include "seqdb.h"

class entropy
	{
public:
	uint m_window = 12;
	uint m_step = 12;
	uint m_nseq = 0;
	double m_max_possible_H = 0;

	vector<string> m_fafns;
	vector<string> m_feature_names;
	vector<uint> m_alpha_sizes;
	vector<vector<vector<uint8_t> > > m_profiles;
	vector<string> m_labels;
	vector<uint> m_seqlengths;
	map<string, uint> m_label2idx;
	float **m_unweighted_logoddsvec = 0;
	float **m_weighted_logoddsvec = 0;
	float *m_weights = 0;
	vector<uint> m_feature_subset;
	uint m_total_col_count = 0;

	map<vector<uint>, double> m_fis2H;

	SeqDB m_TPDB;
	SeqDB m_FPDB;

	vector<uint> m_profidxqs;
	vector<uint> m_profidxts;
	vector<vector<uint> > m_posvecq;
	vector<vector<uint> > m_posvect;
	vector<bool> m_pair_is_tp_vec;

	vector<float> m_scores;
	vector<bool> m_is_tps;

public:
	void set_max_possible_H()
		{
		uint w = m_window;
		double P = 1.0/w;
		m_max_possible_H = w*(-P*log(P));
		}

	double get_entropy(
		const vector<vector<uint8_t> > &profile,
		const vector<uint> &fis,
		uint start_pos) const;

	const char *get_namestr(const vector<uint> &fis, string &s) const;
	void load_profiles(const vector<string> &fafns);
	double get_mean_entropy(const vector<uint> &fis);
	void get_mean_entropy_vec(
		const vector<vector<uint> > &fivec,
		vector<double> &Hs,
		vector<uint> &order);
	void sort_fis(vector<uint> &fis) const;
	void get_next_fis(
		const vector<vector<uint> > &fivec,
		const vector<uint> &order,
		uint maxn,
		vector<vector<uint> > &next_fis) const;

	void load_fa2s(
		const string &tpfa2fn,
		const string &fpfa2fn);

	void parse_fa2(SeqDB &DB, const bool is_tp);

	void read_logoddsvec(const vector<string> &fns);

	float calc_col_score(
		const vector<vector<uint8_t> > &profq, uint posq,
		const vector<vector<uint8_t> > &proft, uint post) const;

	void set_logodds_subset(
		const vector<uint> &fis,
		const vector<float> &weights);

	void calc_col_scores();

	float roc_auc(const vector<float>& scores,
                      const vector<bool>& is_tp) const;

	};
