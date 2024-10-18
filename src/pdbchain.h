#ifndef pdbchain_h
#define pdbchain_h

#include "myutils.h"
#include "features.h"

class DSS;

class PDBChain
	{
public:
	string m_Label;
	string m_Seq;
	vector<float> m_Xs;
	vector<float> m_Ys;
	vector<float> m_Zs;
	uint m_Idx = UINT_MAX;
	vector<string> m_Lines;

public:
	~PDBChain() { Clear(); }
	void Clear()
		{
		m_Label.clear();
		m_Seq.clear();
		m_Xs.clear();
		m_Ys.clear();
		m_Zs.clear();
		m_Lines.clear();
		}

	uint GetSeqLength() const;
	char FromPDBLines(const string &Label,
	  const vector<string> &Lines, bool SaveLines = false);
	void FromCal(const string &FileName);
	void FromCalLines(const vector<string> &Lines);
	void ToCal(FILE *f) const;
	void ToFasta(const string &FileName) const;
	void ToFasta(FILE *f) const;
	void ToFeatureFasta(FILE *f, DSS &D, FEATURE Feat) const;
	void ToCalSeg(FILE *f, uint Pos, uint n) const;
	void ToCal(const string &FileName) const;
	void GetXFormChain_tR(
	  const vector<double> &t,
	  const vector<vector<double> > &R,
	  PDBChain &XChain) const;
	void GetXFormChain_tR(
	  const vector<float> &t,
	  const vector<vector<float> > &R,
	  PDBChain &XChain) const;
	void LogMe(bool WithCoords = false) const;
	void GetXYZ(uint Pos, float &x, float &y, float &z) const;
	void GetPt(uint Pos, vector<float> &Pt) const;
	void GetPt(uint Pos, vector<double> &Pt) const;
	void SetPt(uint Pos, const vector<float> &Pt);
	float GetDist(uint Pos1, uint Pos2) const;
	float GetDist2(uint Pos1, uint Pos2) const;
	void GetSS(string &SS) const;
	void GetICs(vector<uint16_t> &ICs) const;
	void CoordsFromICs(const vector<uint16_t> &ICs);
	void CoordsFromICs(const uint16_t *ICs, uint L);
	void Reverse();

public:
	static bool IsATOMLine(const string &Line);
	static void GetXYZFromATOMLine(const string &InputLine,
	  float &x, float &y, float &z);
	static void SetXYZInATOMLine(const string &InputLine,
	  float x, float y, float z, string &OutputLine);
	static bool GetFieldsFromATOMLine(const string &Line,
	  float &X, float &Y, float &Z, char &aa);
	static uint16_t CoordToIC(float X) { return uint16_t((X + 1000)*10 + 0.5); }
	static float ICToCoord(uint16_t IC) { return float(IC/10.0f) - 1000; }
	};

void ReadChains(const string &FileName, vector<PDBChain *> &Chains);

#endif // pdbchain_h
