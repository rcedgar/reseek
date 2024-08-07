#pragma once

enum SYMMETRY
	{
	SYMM_First,
	SYMM_Both,
	SYMM_Either
	};

class DALIScorer
	{
public:
	string m_Name;
	const SeqDB *m_MSA = 0;
	bool m_DoCore = false;
	vector<PDBChain *> m_Chains;
	map<string, uint> m_SeqToChainIdx;

	vector<uint> m_SeqIdxToChainIdx;
	set<string> m_NotFoundLabels;
	vector<bool> m_ColIsCore;
	uint m_CoreColCount = 0;

	vector<vector<uint> > m_ColToPosVec;

	vector<vector<vector<double> > > m_DistMxVec;
	vector<double> m_LDDT_thresholds;
	double m_LDDT_R0 = DBL_MAX;
	SYMMETRY m_LDDT_symm = SYMM_First;

	double m_DALI_D = 20.0;
	double m_DALI_d0 = 0.2;
	double m_DALI_Theta = 0.2;
	double m_DALI_R0 = DBL_MAX;

public:
	DALIScorer()
		{
		m_LDDT_R0 = 15;

		m_LDDT_thresholds.push_back(0.5);
		m_LDDT_thresholds.push_back(1);
		m_LDDT_thresholds.push_back(2);
		m_LDDT_thresholds.push_back(4);
		}

	void ClearMSA()
		{
		m_MSA = 0;
		m_Name.clear();
		m_SeqIdxToChainIdx.clear();
		m_NotFoundLabels.clear();
		m_ColIsCore.clear();
		}

	void LoadChains(const string &FN);
	void SetMSA(const string &Name, const SeqDB &MSA,
	  bool DoCore, bool MissingSeqOk);
	double GetZ() const;
	double GetZ_Rows() const;
	double GetSumScore_Cols() const;
	double GetSumScore_Rows() const;
	void SetSeqIdxToChainIdx(bool MissingSeqOk);
	void SetCore();
	void SetColToPosVec(bool Core);
	void GetColToPos(uint SeqIdx, vector<uint> &ColToPos, bool Core);
	uint GetColCount() const { return m_MSA->GetColCount(); }
	uint GetSeqCount() const { return m_MSA->GetSeqCount(); }
	bool GetDALIRowPair(uint SeqIdx1, uint SeqIdx2, double &Score, double &Z) const;
	double GetDALIScoreColPair(uint Col1, uint Col2) const;
	double GetDALIScorePosPair_ById(
	  uint ChainIdxX, uint PosX1, uint PosX2,
	  uint ChainIdxY, uint PosY1, uint PosY2) const;
	double GetDiagScore() const;
	double GetDiagScoreSeqPair(uint SeqIdx1, uint SeqIdx2) const;
	double GetDist(uint ChainId, uint Pos1, uint Pos2) const;
	const vector<vector<double> > &GetDistMx(uint ChainIdx) const;
	void SetDistMx(uint ChainId);
	void SetDistMxs();

	double GetDALIScore_ChainPair(uint ChainIdx1, uint ChainIdx2,
	  const vector<uint> &Pos1s, vector<uint> &Pos2s) const;
	double GetDALIScore_OffDiag(uint ChainIdx1, uint ChainIdx2,
	  const vector<uint> &PosQs, const vector<uint> &PosTs) const;

	double GetLDDT_foldmason() const;
	double GetLDDTColScore_foldmason(uint Col) const;

	double GetLDDT_muscle() const;
	double GetLDDTPair_muscle(uint ChainId1, uint ChainId2,
	  const vector<uint> &Pos1s, const vector<uint> &Pos2s) const;
	};

double DALI_dpscorefun(double a, double b);
//double GetDALIScore(const PDBChain &Q, const PDBChain &T,
//	const vector<uint> &PosQs, const vector<uint> &PosTs);
double GetDALIZFromScoreAndLengths(double DALIScore, uint QL, uint TL);
extern float g_DALI_Theta;