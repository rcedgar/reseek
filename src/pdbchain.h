#ifndef pdbchain_h
#define pdbchain_h

// More memory, no speed advantage
#define CACHE_DIST_MAX	0

#include "myutils.h"
#include "features.h"
#include "flatmx.h"
#include "coords.h"

class DSS;

class PDBChain
	{
public:
	string m_Label;

// CIF allows chain strings
// PDB must be single char
	bool m_HasChainStr = false;
	string m_ChainStr;
	string m_Seq;
	vector<float> m_Xs;
	vector<float> m_Ys;
	vector<float> m_Zs;
	uint m_Idx = UINT_MAX;
	vector<string> m_Lines;
#if CACHE_DIST_MAX
	FlatMx<float> *m_DistMx = 0;
#endif

//#if _MSC_VER && DEBUG
//	void *operator new(size_t n)
//		{
//		void *p = _malloc_dbg(n, _CLIENT_BLOCK|(1<<16), __FILE__, __LINE__);
//		return p;
//		}
//
//	void operator delete(void *p) { _free_dbg(p, _CLIENT_BLOCK|(1<<16)); }
//#endif // DEBUG

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
#if CACHE_DIST_MAX
		if (m_DistMx != 0) delete m_DistMx;
#endif
		}

	uint GetSeqLength() const;
	bool FromPDBLines(const string &Label, const vector<string> &Lines,
	  bool SaveLines = false);
	void FromCal(const string &FileName);
	void FromCalLines(const vector<string> &Lines);
	void ToCal(FILE *f) const;
	void ToFasta(const string &FileName) const;
	void ToFasta(FILE *f) const;
	void ToFeatureFasta(FILE *f, DSS &D, FEATURE Feat) const;
	void ToCalSeg(FILE *f, uint Pos, uint n) const;
	void ToCal(const string &FileName) const;
	void ToPDB(FILE *f, bool TruncateAtZ = false) const;
	void ToPDB(const string &FileName) const;
	void GetXFormChain_tR(
	  const vector<float> &t,
	  const vector<vector<float> > &R,
	  PDBChain &XChain) const;
	void LogMe(bool WithCoords = false) const;
	void GetXYZ(uint Pos, float &x, float &y, float &z) const;
	void GetPt(uint Pos, vector<float> &Pt) const;
	void SetPt(uint Pos, const vector<float> &Pt);
	float GetDist(uint Pos1, uint Pos2) const;
	void GetSS(string &SS) const;
	void GetICs(vector<uint16_t> &ICs) const;
	void CoordsFromICs(const vector<uint16_t> &ICs);
	void CoordsFromICs(const uint16_t *ICs, uint L);
	void Reverse();
	void Flip();
	void GetReverse(PDBChain &Rev) const;
	void GetCoords(uint i, coords &v) const
		{
		assert(i < SIZE(m_Xs));
		v.x = m_Xs[i];
		v.y = m_Ys[i];
		v.z = m_Zs[i];
		}

#if CACHE_DIST_MAX
	void SetDistMx();
	void ClearDistMx();
#endif

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
