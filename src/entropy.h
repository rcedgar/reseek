#pragma once

class entropy
	{
public:
	uint m_window = 12;
	uint m_step = 12;
	uint m_nseq = 0;
	double m_max_possible_H = 0;

	vector<string> m_fafns;
	vector<string> m_feature_names;
	vector<vector<vector<uint8_t> > > m_profiles;

	map<vector<uint>, double> m_fis2H;

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
	};