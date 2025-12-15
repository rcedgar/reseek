#pragma once

#include <set>

class PWAlnDB
	{
public:
	vector<string> m_Labels;
	map<string, uint> *m_LabelToIdx = 0;
	vector<uint> m_LabelIdx1s;
	vector<uint> m_LabelIdx2s;
	vector<string> m_Row1s;
	vector<string> m_Row2s;
	vector<vector<uint> > m_ColToPos1Vec;
	vector<vector<uint> > m_ColToPos2Vec;
	bool m_BuildMap = false;
	set<string> m_MissingLabels;

public:
	void Load(const string &FA2FN,
		map<string, uint> &LabelToIdx,
		bool BuildMap);
	uint GetLabelIdx(const string &Label);
	void AddAln(
		const string &Label1, const string &Row1,
		const string &Label2, const string &Row2);
	uint GetTotalColCount() const;
	uint GetAlnCount() const { return SIZE(m_ColToPos1Vec); }
	uint GetLabelIdx1(uint AlnIdx) const;
	uint GetLabelIdx2(uint AlnIdx) const;
	const vector<uint> &GetColToPosVec1(uint AlnIdx) const;
	const vector<uint> &GetColToPosVec2(uint AlnIdx) const;
	};
