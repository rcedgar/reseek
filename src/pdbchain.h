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
	vector<int> m_ResNrs;
	vector<vector<string> > m_ATOMs;
	vector<uint> m_MotifPosVec;
	string m_SS;

public:
	void Clear()
		{
		m_Label.clear();
		m_Seq.clear();
		m_Xs.clear();
		m_Ys.clear();
		m_Zs.clear();
		m_ResNrs.clear();
		m_MotifPosVec.clear();
		m_ATOMs.clear();
		m_SS.clear();
		}

	void GetReverse(PDBChain &Chain) const;
	void WriteSeqWithCoords(FILE *f) const;
	void SetMotifPosVec(uint PosA, uint PosB, uint PosC);
	uint GetSeqLength() const;
	void PrintSeqCoords(FILE *f) const;
	char FromPDBLines(const string &Label,
	  const vector<string> &Lines);
	void FromCal(const string &FileName);
	void FromCalLines(const vector<string> &Lines);
	void ParseCalLabelLine(const string &Line);
	void RenumberResidues(uint Start);
	void ToCal(FILE *f) const;
	void ToCalSeg(FILE *f, uint Pos, uint n) const;
	void ToCal(const string &FileName) const;
	void ToPDB(const string &FileName,
	  const vector<string> *ptrRemarks = 0) const;
	void ToPDB(FILE *f, const vector<string> *ptrRemarks = 0) const;
	void GetXFormChain_tR(
	  const vector<double> &t,
	  const vector<vector<double> > &R,
	  PDBChain &XChain) const;
	void LogMe(bool WithCoords = false) const;
	void GetSubSeq(uint Pos, uint n, string &s) const;
	void GetXYZ(uint Pos, double &x, double &y, double &z) const;
	double GetX(uint Pos) const;
	double GetY(uint Pos) const;
	double GetZ(uint Pos) const;
	bool Get_CB_XYZ(uint Pos, double &x, double &y, double &z) const;
	bool Get_CB_Pt(uint Pos, vector<double> &Pt) const;
	double GetCoord(uint Axis, uint Pos) const;
	void GetPt(uint Pos, vector<double> &Pt) const;
	void SetPt(uint Pos, const vector<double> &Pt);
	double GetDist(uint Pos1, uint Pos2) const;
	double GetDist2(uint Pos1, uint Pos2) const;
	void GetSubSeq(uint StartPos, uint n,
	  bool FailOnOverflow, string &MotifSeq) const;
	void GetSS(string &SS) const;
	void SetSS()
		{
		GetSS(m_SS);
		}
	double GetSmoothedCoord(uint Axis, uint i, uint N, uint w) const;
	const char *GetAcc(string &Acc) const;
	void GetDistMx(uint Pos, uint L, vector<vector<double> > &Mx) const;
	void GetRange(uint Lo, uint Hi, PDBChain &Chain) const;
	void GetCAAtomLine(uint Pos, string &Line) const;
	void GetATOMLines(uint Pos, vector<string> &Lines) const;
	int GetResidueNr(uint Pos, int ValueIfNotFound = INT_MAX) const;
	void GetResidueRange(uint PosLo, uint ResidueCount, int &ResLo,
	  int &ResHi) const;
	void GetSphere(uint Pos, double Radius,
	  uint MinPos, uint MaxPos,
	  vector<uint> &PosVec) const;
	void GetResidueAtomsInfo(uint Pos,
	  vector<double> &Xs,
	  vector<double> &Ys,
	  vector<double> &Zs,
	  vector<string> &ElementNames,
	  vector<string> &AtomNames,
	  vector<string> &Lines) const;

public:
	static void ChainsFromLines(const string &Label,
	  const vector<string> &Lines, vector<PDBChain *> &Chains);
	static void ReadChainsFromFile(const string &FileName,
	  vector<PDBChain *> &Chains);
	static void AppendChainToLabel(string &Label, char Chain);
	static char GetChainCharFromATOMLine(const string &Line);
	static bool IsATOMLine(const string &Line);
	static int GetResidueNrFromATOMLine(const string &Line);
	static void HackHETAMLine(string &Line);
	static void GetAtomNameFromATOMLine(const string &Line,
	  string &AtomName);
	static void GetElementNameFromATOMLine(const string &InputLine,
	  string &ElementName);
	static void SetResidueNrInATOMLine(const string &InputLine,
	  uint ResidueNr, string &OutputLine);
	static void SetAtomNrInATOMLine(const string &InputLine,
	  uint AtomNr, string &OutputLine);
	static void GetXYZFromATOMLine(const string &InputLine,
	  double &x, double &y, double &z);
	static void SetXYZInATOMLine(const string &InputLine,
	  double x, double y, double z, string &OutputLine);
	static bool GetFieldsFromResidueATOMLines(const vector<string> &Lines,
	  double &X, double &Y, double &Z, char &aa, int &ResNr);
	};

void ReadChains(const string &FileName,
  vector<PDBChain *> &Chains);
void ReadChainsFromFileNameVec(const vector<string> &FileNames,
  vector<PDBChain *> &Chains);
void ReadChainsFromDirectory(const string &DirName,
  vector<PDBChain *> &Chains, bool Recursive);
void GetLabelFromFileName(const string &FileName, string &Label);

#endif // pdbchain_h
