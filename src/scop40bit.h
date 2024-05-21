#pragma once

#include <map>

class SCOP40Bit
	{
public:
	vector<string> m_Doms;
	vector<string> m_Fams;
	map<string, uint> m_FamToIdx;
	map<string, uint> m_DomToIdx;
	map<string, string> m_DomToFam;
	vector<uint> m_DomIdxToFamIdx;

	vector<uint> m_DomIdx1s;
	vector<uint> m_DomIdx2s;
	vector<float> m_Scores;

	bool m_ScoresAreEvalues = false;

public:
	SCOP40Bit()
		{
		m_ScoresAreEvalues = opt_scores_are_evalues;
		}

	uint GetHitCount() const;
	uint GetDomCount() const;
	uint GetFamCount() const;
	void ReadDomInfo(const string &ChainsFN);
	void ReadHits_Tsv(const string &Algo);
	void ReadHits_Tsv_DSS();
	void WriteHits_Bin(const string &FileName) const;
	void ReadHits_Bin(const string &FileName);
	bool IsT(uint DomIdx1, uint DomIdx2) const;
	void GetTFs(vector<bool> &TFs) const;
	void GetFamSizes(vector<uint> &FamSizes) const;
	void GetFamSizes_Present(vector<uint> &FamSizes) const;
	void CalcNXs(uint &NT, uint &NF) const;
	void GetOrder(vector<uint> &Order) const;
	void GetROCSteps(vector<double> &ScoreSteps,
	   vector<uint> &TPCounts, vector<uint> &FPCounts) const;
	const char *GetDomName(uint DomIdx) const;
	uint GetDomIdx(const string &DomName) const;
	uint GetSens1stFP(const vector<float> &Scores) const;
	void CalcNXs_Present(uint &NT, uint &NF) const;
	};
