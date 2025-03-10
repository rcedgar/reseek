#pragma once

#include "pdbchain.h"
#include "seqdb.h"

class CMProf
	{
public:
	vector<vector<float> > m_MeanDistMx;
	vector<vector<float> > m_StdDevs;
	const SeqDB *m_MSA = 0;

// Training only
	vector<bool> m_ColIsCore;
	vector<uint> m_CoreCols;
	vector<vector<vector<float> > > m_DistMxVec;
	map<string, uint> m_UngappedSeqToIdx;


public:
	CMProf()
		{
		Clear();
		};

	void Clear()
		{
		m_ColIsCore.clear();
		m_MeanDistMx.clear();
		m_StdDevs.clear();
		m_DistMxVec.clear();
		}
	
	void SetMSA(const SeqDB &MSA);
	uint GetColCount() const { return m_MSA->GetColCount(); }
	uint GetCoreColCount() const { return SIZE(m_CoreCols); }

	void ToFile(const string &FileName) const;
	void MxToFile(FILE *f, const string &Name,
	  const vector<vector<float> > &Mx) const;
	void MxFromFile(FILE *f, string &Name, uint CoreColCount,
	  vector<vector<float> > &Mx);
	void FromFile(FILE *f);
	void FromFile(const string &FileName);

// For training
public:
	void InitTrain()
		{
		m_DistMxVec.clear();
		}

	void FinalizeTrain();
	bool TrainChain(const PDBChain &Chain);
	void GetDistMx(const PDBChain &Chain, const vector<uint> &PosVec,
	  vector<vector<float> > &DistMx);
	void GetMeanStdDev(uint i, uint j,
	  float &Mean, float &StdDev) const;
	};

float GetNormal(float Mu, float Sigma, float x);
