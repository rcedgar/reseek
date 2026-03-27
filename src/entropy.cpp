#include "myutils.h"
#include "flat_chain.h"
#include "flat_helpers.h"
#include "alpha.h"
#include "chaq.h"
#include "seqdb.h"
#include "sort.h"
#include "entropy.h"

double entropy::get_entropy(
	const vector<vector<uint8_t> > &profile,
	const vector<uint> &fis,
	uint start_pos) const
	{
	const size_t nfeat = profile.size();
	const size_t nfi = fis.size();
	map<vector<uint8_t>, uint> unique2count;
	vector<uint8_t> v(nfi);
	for (uint pos = start_pos; pos < start_pos + m_window; ++pos)
		{
		for (uint fidx = 0; fidx < nfi; ++fidx)
			{
			uint fi = fis[fidx];
			v[fidx] = profile[fi][pos];
			}
		map<vector<uint8_t>, uint>::const_iterator iter =
			unique2count.find(v);
		if (iter == unique2count.end())
			unique2count[v] = 1;
		else
			unique2count[v] += 1;
		}

	double H = 0;
	for (map<vector<uint8_t>, uint>::const_iterator iter = unique2count.begin();
		iter != unique2count.end(); ++iter)
		{
		uint count = iter->second;
		double P = double(count)/m_window;
		H += -P*log(P);
		}
	return H;
	}

void entropy::load_profiles(const vector<string> &fafns)
	{
	m_fafns = fafns;
	m_profiles.clear();
	m_feature_names.clear();
	const size_t nfeat = m_fafns.size();

	vector<string> names;
	string namestr;
	for (uint i = 0; i < nfeat; ++i)
		{
		string name;
		GetStemName(m_fafns[i], name);
		ProgressStep(i, nfeat, "Loading profiles %s", name.c_str());
		m_feature_names.push_back(name);
		}

	vector<SeqDB *> DBs;
	const map<string, uint> *label2idx = 0;
	m_nseq = 0;
	for (size_t fi = 0; fi < nfeat; ++fi)
		{
		const string &fafn = m_fafns[fi];
		string feature_name;
		GetStemName(fafn, feature_name);
		uint alpha_size = get_alpha_size_from_feature_name(feature_name);
		const uint8_t *char2code = chaq::get_char2letter(alpha_size);
		SeqDB &DB = *new SeqDB;
		DB.FromFasta(fafn);
		DB.SetLabelToIndex();
		if (fi == 0)
			{
			m_nseq = DB.GetSeqCount();
			label2idx = &DB.m_LabelToIndex;
			m_profiles.resize(m_nseq);
			for (uint seqidx = 0; seqidx < m_nseq; ++seqidx)
				m_profiles[seqidx].resize(nfeat);
			}
		DB.ToLetters(char2code);
		uint labelidx = 0;
		uint undef = 0;
		for (map<string, uint>::const_iterator iter = label2idx->begin();
			iter != label2idx->end(); ++iter)
			{
			const string &label = iter->first;
			uint seqidx = DB.GetSeqIndex(label);
			const uint L = DB.GetSeqLength(seqidx);
			const uint8_t *seq = DB.GetByteSeq(seqidx);
			for (uint pos = 0; pos < L; ++pos)
				{
				uint8_t code = seq[pos];
				if (code >= alpha_size)
					{
					++undef;
					code = randu32()%alpha_size;
					}
				m_profiles[labelidx][fi].push_back(code);
				}
			++labelidx;
			}
		if (undef > 0)
			ProgressLog("%s undef=%u\n", fafn.c_str(), undef);
		}
	}

void entropy::sort_fis(vector<uint> &fis) const
	{
	sort(fis.begin(), fis.end());
	}

double entropy::get_mean_entropy(const vector<uint> &fis)
	{
	vector<uint> sorted_fis(fis);
	sort_fis(sorted_fis);
	string namestr;
	get_namestr(sorted_fis, namestr);
	map<vector<uint>, double>::const_iterator iter =
		m_fis2H.find(sorted_fis);
	if (iter != m_fis2H.end())
		{
		double meanH = iter->second;
		ProgressLog("n=%u;%s=%.4f (cached)\n",
			uint(fis.size()), namestr.c_str(), meanH);
		return meanH;
		}

	uint N = 0;
	double sumH = 0;
	for (uint seqidx = 0; seqidx < m_nseq; ++seqidx)
		{
		const vector<vector<uint8_t> > &profile = m_profiles[seqidx];
		const size_t L = profile[0].size();
		for (uint pos = 0; pos + m_window <= L; pos += m_step)
			{
			double H = get_entropy(profile, fis, pos);
			sumH += H;
			++N;
			}
		}
	double meanH = sumH/N;
	m_fis2H[sorted_fis] = meanH;

	ProgressLog("n=%u;%s=%.4f\n",
		uint(fis.size()), namestr.c_str(), meanH);
	return meanH;
	}

void cmd_entropy()
	{
	const string &fafnstr = g_Arg1;
	vector<string> m_fafns;
	Split(fafnstr, m_fafns, ',');
	const size_t nfeat = m_fafns.size();

	entropy E;

	E.m_window = 12;
	if (optset_window)
		E.m_window = opt(window);
	E.m_step = 12;
	if (optset_step)
		E.m_step = opt(step);

	E.load_profiles(m_fafns);

	vector<uint> fis;
	vector<string> names;
	string namestr;
	for (uint i = 0; i < nfeat; ++i)
		{
		fis.push_back(i);
		string name;
		GetStemName(m_fafns[i], name);
		names.push_back(name);
		if (i > 0)
			namestr += ',';
		namestr += name;
		}

	double meanH = E.get_mean_entropy(fis);
	ProgressLog("names=%s H=%.4f\n", namestr.c_str(), meanH);
	}

const char *entropy::get_namestr(const vector<uint> &fis, string &s) const
	{
	s.clear();
	for (int i = 0; i < fis.size(); ++i)
		{
		if (i > 0)
			s += ",";
		uint fi = fis[i];
		asserta(fi < m_feature_names.size());
		s += m_feature_names[fi];
		}
	return s.c_str();
	}

void entropy::get_mean_entropy_vec(
	const vector<vector<uint> > &fivec,
	vector<double> &Hs,
	vector<uint> &order)
	{
	Hs.clear();
	order.clear();
	const uint n = uint(fivec.size());
	for (uint i = 0; i < n; ++i)
		{
		const vector<uint> &fis = fivec[i];
		double H = get_mean_entropy(fis);
		Hs.push_back(H);
		}
	order.resize(n);
	QuickSortOrderDesc(Hs.data(), n, order.data());
	}

void entropy::get_next_fis(
	const vector<vector<uint> > &fivec,
	const vector<uint> &order,
	uint maxn,
	vector<vector<uint> > &next_fis) const
	{
	next_fis.clear();
	const uint nfeat = uint(m_feature_names.size());
	set<vector<uint> > newset;
	for (size_t i = 0; i < uint(order.size()); ++i)
		{
		const vector<uint> &fis = fivec[order[i]];
		for (uint fi = 0; fi < nfeat; ++fi)
			{
			bool found = false;
			for (uint j = 0; j < uint(fis.size()); ++j)
				{
				if (fis[j] == fi)
					{
					found = true;
					break;
					}
				}
			if (!found)
				{
				vector<uint> next = fis;
				next.push_back(fi);
				sort_fis(next);
				map<vector<uint>, double>::const_iterator iter =
					m_fis2H.find(next);
				if (iter == m_fis2H.end())
					{
					if (newset.find(next) == newset.end())
						{
						next_fis.push_back(next);
						newset.insert(next);
						if (newset.size() >= maxn)
							return;
						}
					}
				}
			}
		}
	}

void cmd_entropy_greedy()
	{
	const string &filesfn = g_Arg1;
	vector<string> fafns;
	ReadLinesFromFile(filesfn, fafns);
	const uint nfeat = uint(fafns.size());

	entropy E;
	E.m_window = 12;
	if (optset_window)
		E.m_window = opt(window);
	E.m_step = 12;
	if (optset_step)
		E.m_step = opt(step);
	ProgressLog("window=%u, step=%u\n",
		E.m_window, E.m_step);

	E.load_profiles(fafns);

	vector<uint> fis;
	fis.push_back(0);

	double H0 = E.get_mean_entropy(fis);
	vector<double> Hs;
	vector<vector<uint> > fivec;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		fis.clear();
		fis.push_back(fi);
		fivec.push_back(fis);
		}
	vector<uint> order;
	E.get_mean_entropy_vec(fivec, Hs, order);

	vector<vector<uint> > last_fivec(fivec);
	vector<vector<uint> > next_fivec(fivec);
	for (uint iter = 0; iter < 100; ++iter)
		{
		E.get_next_fis(last_fivec, order, 1024, next_fivec);
		if (next_fivec.empty())
			break;
		Progress("\niter %u (%u)\n", iter+1, uint(next_fivec.size()));
		E.get_mean_entropy_vec(next_fivec, Hs, order);
		last_fivec = next_fivec;
		}
	Progress("converged.\n");
	}
