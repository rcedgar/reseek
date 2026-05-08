#pragma once


class collect
	{
public:
	string m_name;
	unordered_map<string, vector<string> > m_name2lines;

public:
	void from_file(const string &fn);
	const vector<string> &get_lines(const string &name) const;
	};