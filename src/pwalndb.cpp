#include "myutils.h"
#include "pwalndb.h"
#include "sfasta.h"
#include "alpha.h"

void TruncLabel(string &Label)
	{
	size_t n = Label.find(' ');
	if (n != string::npos)
		Label.resize(n);
	n = Label.find('|');
	if (n != string::npos)
		Label.resize(n);
	n = Label.find('/');
	if (n != string::npos)
		Label.resize(n);
	}

uint PWAlnDB::GetLabelIdx(const string &Label)
	{
	asserta(!Label.empty());
	map<string, uint>::const_iterator iter =
		m_LabelToIdx->find(Label);
	if (iter == m_LabelToIdx->end())
		{
		if (m_BuildMap)
			{
			uint Idx = SIZE(m_Labels);
			m_Labels.push_back(Label);
			return Idx;
			}
		else
			{
			m_MissingLabels.insert(Label);
			return UINT_MAX;
			}
		}
	return iter->second;
	}

void PWAlnDB::Load(const string &FA2FN,
	map<string, uint> &LabelToIdx,
	bool BuildMap)
	{
	bool MapIsEmpty = LabelToIdx.empty();
	asserta(BuildMap == MapIsEmpty);
	m_BuildMap = BuildMap;
	if (BuildMap)
		m_LabelToIdx = new map<string, uint>;
	else
		m_LabelToIdx = &LabelToIdx;

	SFasta SF;
	SF.Open(FA2FN);
	SF.m_AllowGaps = true;

	bool FirstRow = true;
	string Label1;
	string Row1;
	uint LabelIdx1 = UINT_MAX;
	for (;;)
		{
		const char *Seq = SF.GetNextSeq();
		if (Seq == 0)
			break;
		const unsigned L = SF.GetSeqLength();
		asserta(L != 0);
		string Row;
		Row.reserve(L);
		for (uint i = 0; i < L; ++i)
			Row.push_back(Seq[i]);
		string Label = SF.GetLabel();
		TruncLabel(Label);
		uint LabelIdx = GetLabelIdx(Label);
		if (FirstRow)
			{
			Row1 = Row;
			Label1 = Label;
			LabelIdx1 = LabelIdx;
			}
		else
			AddAln(Label1, Row1, Label, Row);
		FirstRow = !FirstRow;
		}
	asserta(FirstRow);
	uint N = SIZE(m_Row1s);
	asserta(N%2 == 0);
	asserta(SIZE(m_Row2s) == N);
	asserta(SIZE(m_LabelIdx1s) == N);
	asserta(SIZE(m_LabelIdx2s) == N);
	asserta(SIZE(m_ColToPos1Vec) == N);
	asserta(SIZE(m_ColToPos2Vec) == N);
	}

void PWAlnDB::AddAln(
	const string &Label1, const string &Row1,
	const string &Label2, const string &Row2)
	{
	uint LabelIdx1 = GetLabelIdx(Label1);
	uint LabelIdx2 = GetLabelIdx(Label2);
	if (LabelIdx1 == UINT_MAX || LabelIdx2 == UINT_MAX)
		return;
	const uint ColCount = SIZE(Row1);
	asserta(SIZE(Row2) == ColCount);
	vector<uint> ColToPos1;
	vector<uint> ColToPos2;
	uint Pos1 = 0;
	uint Pos2 = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c1 = Row1[Col];
		char c2 = Row2[Col];

		if (isupper(c1) && isupper(c2))
			{
			ColToPos1.push_back(Pos1);
			ColToPos2.push_back(Pos2);
			}

		if (!isgap(c1))
			++Pos1;
		if (!isgap(c2))
			++Pos2;
		}
	m_LabelIdx1s.push_back(LabelIdx1);
	m_LabelIdx2s.push_back(LabelIdx2);
	m_Row1s.push_back(Row1);
	m_Row2s.push_back(Row2);
	m_ColToPos1Vec.push_back(ColToPos1);
	m_ColToPos2Vec.push_back(ColToPos2);
	}

uint PWAlnDB::GetTotalColCount() const
	{
	uint Total = 0;
	for (uint i = 0; i < SIZE(m_ColToPos1Vec); ++i)
		Total += SIZE(m_ColToPos1Vec[i]);
	return Total;
	}

uint PWAlnDB::GetLabelIdx1(uint AlnIdx) const
	{
	asserta(AlnIdx < SIZE(m_LabelIdx1s));
	return m_LabelIdx1s[AlnIdx];
	}

uint PWAlnDB::GetLabelIdx2(uint AlnIdx) const
	{
	asserta(AlnIdx < SIZE(m_LabelIdx2s));
	return m_LabelIdx2s[AlnIdx];
	}

const vector<uint> &PWAlnDB::GetColToPosVec1(uint AlnIdx) const
	{
	asserta(AlnIdx < SIZE(m_ColToPos1Vec));
	return m_ColToPos1Vec[AlnIdx];
	}

const vector<uint> &PWAlnDB::GetColToPosVec2(uint AlnIdx) const
	{
	asserta(AlnIdx < SIZE(m_ColToPos2Vec));
	return m_ColToPos2Vec[AlnIdx];
	}
