#include "myutils.h"
#include "lookup.h"
#include "flat_chain.h"
#include "flat_helpers.h"

void lookup::from_tsv(const string &fn)
	{
	vector<string> labels;

	string line;
	vector<string> flds;
	FILE *f = OpenStdioFile(fn);
	while (ReadLineStdioFile(f, line))
		{
		Split(line, flds, '\t');
		asserta(flds.size() == 2);
		labels.push_back(flds[0] + "/" + flds[1]);
		}
	CloseStdioFile(f);
	from_labels(labels);
	}

void lookup::from_labels(const vector<string> &labels)
	{
	clear();
	reserve();

	string line;
	vector<string> flds;
	vector<string> flds2;
	for (auto label : labels)
		{
		uint domidx = uint(m_doms.size());
		Split(label, flds, '/');
		asserta(flds.size() == 2);
		const string &dom = flds[0];
		asserta(m_dom2idx.find(dom) == m_dom2idx.end());

		const string &fam = flds[1];
		Split(fam, flds2, '.');
		asserta(flds2.size() >= 3);
		const string fold = flds2[0] + string(".") + flds2[1];
		const string sf = fold + string(".") + flds2[2];
		uint famidx = UINT_MAX;
		uint sfidx = UINT_MAX;
		uint foldidx = UINT_MAX;

		unordered_map<string, uint>::const_iterator iterfam =
			m_fam2idx.find(fam);
		if (iterfam == m_fam2idx.end())
			{
			famidx = uint(m_fams.size());
			m_fams.push_back(fam);
			m_fam2idx[fam] = famidx;
			}
		else
			famidx = iterfam->second;
		asserta(famidx < m_fams.size());

		unordered_map<string, uint>::const_iterator itersf =
			m_sf2idx.find(sf);
		if (itersf == m_sf2idx.end())
			{
			sfidx = uint(m_sfs.size());
			m_sfs.push_back(sf);
			m_sf2idx[sf] = sfidx;
			}
		else
			sfidx = itersf->second;
		asserta(sfidx < m_sfs.size());

		unordered_map<string, uint>::const_iterator iterfold =
			m_fold2idx.find(fold);
		if (iterfold == m_fold2idx.end())
			{
			foldidx = uint(m_folds.size());
			m_folds.push_back(fold);
			m_fold2idx[fold] = foldidx;
			}
		else
			foldidx = iterfold->second;
		asserta(foldidx < m_folds.size());

		m_doms.push_back(dom);
		m_dom2idx[dom] = domidx;
		m_domidx2famidx.push_back(famidx);
		m_domidx2sfidx.push_back(sfidx);
		m_domidx2foldidx.push_back(foldidx);
		}
	size_t ndom = m_doms.size();
	asserta(m_dom2idx.size() == ndom);
	asserta(m_sf2idx.size() == m_sfs.size());
	asserta(m_fold2idx.size() == m_folds.size());

	fill();
	}

void lookup::to_tsv(const string &fn)
	{
	if (fn == "") return;
	FILE *f = CreateStdioFile(fn);
	for (size_t domidx = 0; domidx < m_doms.size(); ++domidx)
		{
		asserta(domidx < m_domidx2sfidx.size());
		uint sfidx = m_domidx2sfidx[domidx];
		asserta(sfidx < m_sfs.size());
		fprintf(f, "%s\t%s\n", m_doms[domidx].c_str(), m_sfs[sfidx].c_str());
		}
	CloseStdioFile(f);
	}

void lookup::fill_sf()
	{
	const uint ndom = uint(m_doms.size());
	const uint nsf = uint(m_sfs.size());
	m_sfidx2ndom.clear();
	m_sfidx2ndom.resize(nsf);
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		uint sfidx = m_domidx2sfidx[domidx];
		asserta(sfidx < nsf);
		++m_sfidx2ndom[sfidx];
		}

	m_NT = 0;
	m_NF = 0;
	for (uint sfidx = 0; sfidx < nsf ; ++sfidx)
		{
		uint sfndom = m_sfidx2ndom[sfidx];
		asserta(sfndom > 0);
		m_NT += (sfndom*(sfndom - 1))/2;
		}
	m_NT *= 2;
	m_NF = ndom*(ndom-1) - m_NT;
	m_NI = 0;

	m_pair_count = get_pair_count_upper_triangle_with_diagonal();
	uint k = 0;
	uint NTcheck = 0;
	uint NFcheck = 0;
	for (uint i = 0; i < ndom; ++i)
		{
		for (uint j = i; j < ndom; ++j)
			{
			assert(k < m_pair_count);
			assert(k == triangle_ij_to_k(i, j, ndom));
			if (is_tp_ij(i, j))
				{
				if (i != j)
					++NTcheck;
				}
			else
				++NFcheck;
			++k;
			}
		}
	asserta(k == m_pair_count);
	asserta(NTcheck*2 == m_NT);
	asserta(NFcheck*2 == m_NF);
	stats();
	}

void lookup::fill_fam()
	{
	const uint ndom = uint(m_doms.size());
	const uint nfam = uint(m_fams.size());
	m_famidx2ndom.clear();
	m_famidx2ndom.resize(nfam);
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		uint famidx = m_domidx2famidx[domidx];
		asserta(famidx < nfam);
		++m_famidx2ndom[famidx];
		}

	m_NT = 0;
	m_NF = 0;
	for (uint famidx = 0; famidx < nfam ; ++famidx)
		{
		uint famndom = m_famidx2ndom[famidx];
		asserta(famndom > 0);
		m_NT += (famndom*(famndom - 1))/2;
		}
	m_NT *= 2;
	m_NF = ndom*(ndom-1) - m_NT;
	m_NI = 0;

	m_pair_count = get_pair_count_upper_triangle_with_diagonal();
	uint k = 0;
	uint NTcheck = 0;
	uint NFcheck = 0;
	for (uint i = 0; i < ndom; ++i)
		{
		for (uint j = i; j < ndom; ++j)
			{
			assert(k < m_pair_count);
			assert(k == triangle_ij_to_k(i, j, ndom));
			if (is_tp_ij(i, j))
				{
				if (i != j)
					++NTcheck;
				}
			else
				++NFcheck;
			++k;
			}
		}
	asserta(k == m_pair_count);
	asserta(NTcheck*2 == m_NT);
	asserta(NFcheck*2 == m_NF);
	stats();
	}

void lookup::fill_sfx()
	{
	const uint ndom = uint(m_doms.size());
	const uint nsf = uint(m_sfs.size());
	m_sfidx2ndom.clear();
	m_sfidx2ndom.resize(nsf);
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		uint sfidx = m_domidx2sfidx[domidx];
		asserta(sfidx < nsf);
		++m_sfidx2ndom[sfidx];
		}

	const uint nfold = uint(m_folds.size());
	m_foldidx2ndom.clear();
	m_foldidx2ndom.resize(nfold);
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		uint foldidx = m_domidx2foldidx[domidx];
		asserta(foldidx < nfold);
		++m_foldidx2ndom[foldidx];
		}

	m_NT = 0;
	m_NF = 0;
	m_NI = 0;
	m_pair_count = 0;
	for (uint i = 0; i < ndom; ++i)
		{
		for (uint j = i; j < ndom; ++j)
			{
			if (is_ignored_ij(i, j))
				continue;
			if (is_tp_ij(i, j))
				{
				if (i != j)
					++m_NT;
				}
			else
				++m_NF;
			++m_pair_count;
			}
		}
	stats();
	}

uint lookup::get_singleton_count() const
	{
	uint n = 0;
	switch (m_LT)
		{
	case LT_SAME_FAM:
		{
		vector<uint> famidx2ndom(m_fams.size());
		for (size_t domidx = 0; domidx < m_domidx2famidx.size(); ++domidx)
			{
			uint famidx = m_domidx2famidx[domidx];
			famidx2ndom[famidx] += 1;
			}
		for (size_t famidx = 0; famidx < m_fams.size(); ++famidx)
			if (famidx2ndom[famidx] == 1) ++n;
		return n;
		}

	case LT_SAME_SF:
	case LT_TOP_SF:
		{
		vector<uint> sfidx2ndom(m_sfs.size());
		for (size_t domidx = 0; domidx < m_domidx2sfidx.size(); ++domidx)
			{
			uint sfidx = m_domidx2sfidx[domidx];
			sfidx2ndom[sfidx] += 1;
			}
		for (size_t sfidx = 0; sfidx < m_sfs.size(); ++sfidx)
			if (sfidx2ndom[sfidx] == 1) ++n;
		return n;
		}

	case LT_SAME_FOLD:
	case LT_TOP_FOLD:
		{
		vector<uint> foldidx2ndom(m_folds.size());
		for (size_t domidx = 0; domidx < m_domidx2foldidx.size(); ++domidx)
			{
			uint foldidx = m_domidx2foldidx[domidx];
			foldidx2ndom[foldidx] += 1;
			}
		for (size_t foldidx = 0; foldidx < m_folds.size(); ++foldidx)
			if (foldidx2ndom[foldidx] == 1) ++n;
		return n;
		}

	case LT_DIFF_FAM_SAME_SF:
		{
		vector<set<uint> > sfidx2famset(m_sfs.size());
		for (size_t domidx = 0; domidx < m_domidx2foldidx.size(); ++domidx)
			{
			uint sfidx = m_domidx2sfidx[domidx];
			uint famidx = m_domidx2famidx[domidx];
			sfidx2famset[sfidx].insert(famidx);
			}
		for (size_t domidx = 0; domidx < m_domidx2foldidx.size(); ++domidx)
			{
			uint sfidx = m_domidx2sfidx[domidx];
			if (sfidx2famset[sfidx].size() == 1) ++n;
			}
		return n;
		}

	case LT_DIFF_SF_SAME_FOLD:
		{
		vector<set<uint> > foldidx2sfset(m_folds.size());
		for (size_t domidx = 0; domidx < m_domidx2foldidx.size(); ++domidx)
			{
			uint foldidx = m_domidx2foldidx[domidx];
			uint sfidx = m_domidx2sfidx[domidx];
			foldidx2sfset[foldidx].insert(sfidx);
			}
		for (size_t domidx = 0; domidx < m_domidx2foldidx.size(); ++domidx)
			{
			uint foldidx = m_domidx2foldidx[domidx];
			if (foldidx2sfset[foldidx].size() == 1) ++n;
			}
		return n;
		}

	default:
		break;
		}

	asserta(false);
	return 0;
	}

void lookup::fill_foldx()
	{
	const uint ndom = uint(m_doms.size());
	const uint nsf = uint(m_sfs.size());
	m_sfidx2ndom.clear();
	m_sfidx2ndom.resize(nsf);
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		uint sfidx = m_domidx2sfidx[domidx];
		asserta(sfidx < nsf);
		++m_sfidx2ndom[sfidx];
		}

	const uint nfold = uint(m_folds.size());
	m_foldidx2ndom.clear();
	m_foldidx2ndom.resize(nfold);
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		uint foldidx = m_domidx2foldidx[domidx];
		asserta(foldidx < nfold);
		++m_foldidx2ndom[foldidx];
		}

	m_NT = 0;
	m_NF = 0;
	m_NI = 0;
	m_pair_count = 0;
	for (uint i = 0; i < ndom; ++i)
		{
		for (uint j = i; j < ndom; ++j)
			{
			if (is_ignored_ij(i, j))
				continue;
			if (is_tp_ij(i, j))
				{
				if (i != j)
					++m_NT;
				}
			else
				++m_NF;
			++m_pair_count;
			}
		}
	stats();
	}

void lookup::fill_fold()
	{
	const uint ndom = uint(m_doms.size());
	const uint nfold = uint(m_folds.size());
	m_foldidx2ndom.clear();
	m_foldidx2ndom.resize(nfold);
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		uint foldidx = m_domidx2foldidx[domidx];
		asserta(foldidx < nfold);
		++m_foldidx2ndom[foldidx];
		}

	m_NT = 0;
	m_NF = 0;
	m_NI = 0;
	for (uint foldidx = 0; foldidx < nfold ; ++foldidx)
		{
		uint foldndom = m_foldidx2ndom[foldidx];
		asserta(foldndom > 0);
		m_NT += (foldndom*(foldndom - 1))/2;
		}
	m_NT *= 2;
	m_NF = ndom*(ndom-1) - m_NT;

	m_pair_count = get_pair_count_upper_triangle_with_diagonal();
	uint k = 0;
	uint NTcheck = 0;
	uint NFcheck = 0;
	for (uint i = 0; i < ndom; ++i)
		{
		for (uint j = i; j < ndom; ++j)
			{
			assert(k < m_pair_count);
			assert(k == triangle_ij_to_k(i, j, ndom));
			if (is_tp_ij(i, j))
				{
				if (i != j)
					++NTcheck;
				}
			else
				{
				++NFcheck;
				}
			++k;
			}
		}
	asserta(k == m_pair_count);
	asserta(NTcheck*2 == m_NT);
	asserta(NFcheck*2 == m_NF);
	stats();
	}

void lookup::fill()
	{
	if (m_LT == LT_SAME_FAM)
		fill_fam();
	else if (m_LT == LT_SAME_SF || m_LT == LT_TOP_SF)
		fill_sf();
	else if (m_LT == LT_DIFF_FAM_SAME_SF)
		fill_sfx();
	else if (m_LT == LT_SAME_FOLD || m_LT == LT_TOP_FOLD)
		fill_fold();
	else if (m_LT == LT_DIFF_SF_SAME_FOLD)
		fill_foldx();
	else
		Die("fill");
	}

void lookup::stats()
	{
	ProgressLog("lookup: doms=%u", uint(m_doms.size()));
	ProgressLog(" SFs=%u", uint(m_sfs.size()));
	ProgressLog(" NT=%u", m_NT);
	ProgressLog(" NF=%u", m_NF);
	ProgressLog("\n");
	}

void lookup::sort_chains(
	const vector<flat_chain_t *> &chains,
	vector<flat_chain_t *> &sorted_chains) const
	{
	const uint ndom = get_ndom();
	const uint nchain = uint(chains.size());
	sorted_chains.clear();
	sorted_chains.resize(ndom, nullptr);
	uint nfound = 0;
	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		flat_chain_t *chain = chains[chainidx];
		const string &label = chain->m_label;
		uint domidx = get_domidx(label, true);
		if (domidx == UINT_MAX) continue;
		++nfound;
		asserta(sorted_chains[domidx] == 0);
		sorted_chains[domidx] = chain;
		}
	if (nfound != ndom)
		Die("nfound=%u ndom=%u", nfound, ndom);
	}
