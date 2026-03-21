#include "myutils.h"
#include "lookup.h"

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

		const string &scopid = flds[1];
		Split(scopid, flds2, '.');
		asserta(flds2.size() >= 3);
		const string sf =
			flds2[0] + string(".") + 
			flds2[1] + string(".") + flds2[2];
		uint sfidx = UINT_MAX;
		unordered_map<string, uint>::const_iterator iter =
			m_sf2idx.find(sf);
		if (iter == m_sf2idx.end())
			{
			sfidx = uint(m_sfs.size());
			m_sfs.push_back(sf);
			m_sf2idx[sf] = sfidx;
			}
		else
			sfidx = iter->second;
		asserta(sfidx < m_sfs.size());

		m_doms.push_back(dom);
		m_dom2idx[dom] = domidx;
		m_domidx2sfidx.push_back(sfidx);
		}
	size_t ndom = m_doms.size();
	asserta(m_dom2idx.size() == ndom);
	asserta(m_sf2idx.size() == m_sfs.size());

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

void lookup::fill()
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

	m_pair_count = get_pair_count_upper_triangle_with_diagonal();
	if (m_tpvec != 0) myfree(m_tpvec);
	m_tpvec = myalloc(bool, m_pair_count);
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
				m_tpvec[k] = true;
				if (i != j)
					++NTcheck;
				}
			else
				{
				m_tpvec[k] = false;
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

void lookup::stats()
	{
	ProgressLog("lookup: doms=%u", uint(m_doms.size()));
	ProgressLog(" SFs=%u", uint(m_sfs.size()));
	ProgressLog(" NT=%u", m_NT);
	ProgressLog(" NF=%u", m_NF);
	ProgressLog("\n");
	}
