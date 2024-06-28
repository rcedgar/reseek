#ifndef pdbchain_h
#define pdbchain_h

#include "myutils.h"

class PDBChain
	{
public:
	string m_Label;
	string m_Seq;
	vector<double> m_Xs;
	vector<double> m_Ys;
	vector<double> m_Zs;

public:
	~PDBChain() { Clear(); }
	void Clear()
		{
		m_Label.clear();
		m_Seq.clear();
		m_Xs.clear();
		m_Ys.clear();
		m_Zs.clear();
		}

	uint GetSeqLength() const;
	char FromPDBLines(const string &Label,
	  const vector<string> &Lines);
	void FromCal(const string &FileName);
	void FromCalLines(const vector<string> &Lines);
	void ToCal(FILE *f) const;
	void ToFasta(const string &FileName) const;
	void ToFasta(FILE *f) const;
	void ToCalSeg(FILE *f, uint Pos, uint n) const;
	void ToCal(const string &FileName) const;
	void GetXFormChain_tR(
	  const vector<double> &t,
	  const vector<vector<double> > &R,
	  PDBChain &XChain) const;
	void LogMe(bool WithCoords = false) const;
	void GetXYZ(uint Pos, double &x, double &y, double &z) const;
	void GetPt(uint Pos, vector<double> &Pt) const;
	void SetPt(uint Pos, const vector<double> &Pt);
	double GetDist(uint Pos1, uint Pos2) const;
	double GetDist2(uint Pos1, uint Pos2) const;
	void GetSS(string &SS) const;

public:
	static bool IsATOMLine(const string &Line);
	static void GetXYZFromATOMLine(const string &InputLine,
	  double &x, double &y, double &z);
	static void SetXYZInATOMLine(const string &InputLine,
	  double x, double y, double z, string &OutputLine);
	static bool GetFieldsFromATOMLine(const string &Line,
	  double &X, double &Y, double &Z, char &aa);
	};

void ReadChains(const string &FileName, vector<PDBChain *> &Chains);

#endif // pdbchain_h
