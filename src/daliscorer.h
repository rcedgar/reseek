#pragma once

class DALIScorer
	{
public:
	string m_Name;
	const SeqDB *m_MSA = 0;
	vector<PDBChain *> m_Chains;
	map<string, uint> m_SeqToChainIdx;

	vector<uint> m_SeqIdxToChainIdx;
	vector<string> m_MSALabels;
	vector<string> m_ChainLabels;
	set<string> m_NotFoundLabels;
	vector<bool> m_ColIsCore;

	vector<string> m_MSALabels1;
	vector<string> m_MSALabels2;
	vector<double> m_Scores;
	vector<double> m_Scores_core;
	vector<double> m_Zs;
	vector<double> m_Zs_core;
	double m_SumScore = DBL_MAX;
	double m_SumScore_core = DBL_MAX;
	double m_MeanScore = DBL_MAX;
	double m_MeanScore_core = DBL_MAX;
	double m_MeanZ = DBL_MAX;
	double m_MeanZ_core = DBL_MAX;

	vector<vector<uint> > m_ColToPosVec;
	vector<double> m_ColScores;

public:
	void ClearMSA()
		{
		m_MSA = 0;
		m_Name.clear();
		m_SeqIdxToChainIdx.clear();
		m_MSALabels.clear();
		m_ChainLabels.clear();
		m_NotFoundLabels.clear();
		m_ColIsCore.clear();

		m_MSALabels1.clear();
		m_MSALabels2.clear();
		m_Scores.clear();
		m_Scores_core.clear();
		m_Zs.clear();
		m_Zs_core.clear();
		m_SumScore = DBL_MAX;
		m_SumScore_core = DBL_MAX;
		m_MeanScore = DBL_MAX;
		m_MeanScore_core = DBL_MAX;
		m_MeanZ = DBL_MAX;
		m_MeanZ_core = DBL_MAX;
		}

	void LoadChains(const string &FN);
	void Run(const string &Name, const SeqDB &MSA);
	void RunCols(const string &Name, const SeqDB &MSA, bool DoCore);
	void RunRows(const string &Name, const SeqDB &MSA, bool DoCore);
	void SetSeqIdxToChainIdx();
	void SetCore();
	void ToFev(FILE *f) const;
	void SetColToPosVec(bool Core);
	void GetColToPos(uint SeqIdx, vector<uint> &ColToPos, bool Core);
	uint GetColCount() const { return m_MSA->GetColCount(); }
	uint GetSeqCount() const { return m_MSA->GetSeqCount(); }
	double GetDALIColScore(uint Col1, uint Col2) const;
	double GetDALIPosPairScore(
	  const PDBChain &ChainX, uint PosX1, uint PosX2,
	  const PDBChain &ChainY, uint PosY1, uint PosY2) const;
	double GetDiagScore() const;
	double GetDiagScoreSeqPair(uint SeqIdx1, uint SeqIdx2) const;
	};

double DALI_dpscorefun(double a, double b);
double GetDALIScore(const PDBChain &Q, const PDBChain &T,
	const vector<uint> &PosQs, const vector<uint> &PosTs);
double GetDALIZFromScoreAndLengths(double DALIScore, uint QL, uint TL);
extern float g_DALI_Theta;