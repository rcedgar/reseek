#include "myutils.h"
#include "collect.h"

void collect::from_file(const string &fn)
	{
	vector<string> lines;
	vector<string> flds;
	ReadLinesFromFile(fn, lines);
	const size_t N = lines.size();
	size_t i = 0;
	while (i < N)
		{
		const string &line = lines[i];
		Split(lines[i], flds, '\t');
		asserta(flds.size() == 3);
		asserta(flds[0] == "@");
		const string &name = flds[1];
		uint n = StrToUint(flds[2]);
		asserta(i + n <= N);
		vector<string> name_lines;
		for (uint j = 0; j < n; ++j)
			name_lines.push_back(lines[i+j+1]);
		m_name2lines[name] = name_lines;
		i += n + 1;
		}
	}

const vector<string> &collect::get_lines(const string &name) const
	{
	unordered_map<string, vector<string> >::const_iterator
		iter = m_name2lines.find(name);
	if (iter == m_name2lines.end())
		Die("collect::get_lines(%s)", name.c_str());
	return iter->second;
	}
