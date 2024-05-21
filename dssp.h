#pragma once

class SSEl
	{
public:
	char m_SS = '?';
	string m_Seq;
	uint m_Top = UINT_MAX;
	uint m_Left = UINT_MAX;
	bool m_Fwd = false;
	};

class DSSP
	{
public:
	string m_Label;
	vector<uint> m_SequentialResNrs;
	vector<string> m_OriginalResNrs;
	char m_ChainId = 0;
	string m_Seq;
	string m_SS;	// HBEGIPTSL
	vector<int> m_BetaPartnerResNr1s;
	vector<int> m_BetaPartnerResNr2s;
	vector<uint> m_SolventAccs;

public:
	void Clear()
		{
		m_Label.clear();
		m_SequentialResNrs.clear();
		m_OriginalResNrs.clear();
		m_ChainId = 0;
		m_Seq.clear();
		m_SS.clear();
		m_BetaPartnerResNr1s.clear();
		m_BetaPartnerResNr2s.clear();
		m_SolventAccs.clear();
		}

	void LogMe() const;
	void PrintSeq(FILE *f) const;
	void FromLines(const string &Label,
	  const vector<string> &Lines);
	uint GetSeqLength() const { return SIZE(m_Seq); }
	void GetSSEls(vector<SSEl *> &Els) const;
	void GetBetaPairs(vector<uint> &Start1s, vector<uint> &End1s,
	  vector<uint> &Start2s, vector<uint> &End2s) const;
	void LogBetaPairs(const vector<uint> &Start1s, const vector<uint> &End1s,
	  const vector<uint> &Start2s, const vector<uint> &End2s) const;

public:
	static char GetChainIdFromLine(const string &Line);
	};
