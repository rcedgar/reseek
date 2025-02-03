#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"
#include "abcxyz.h"

void GetThreeFromOne(char aa, string &AAA);

/***
PDBChain ATOM record format
http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.

MODEL     1                                                                     

          1         2         3         4         5         6         7         
01234567890123456789012345678901234567890123456789012345678901234567890123456789
                              xxxxxxxxyyyyyyyyzzzzzzzz
ATOM      1  N   PHE A   1      34.582  19.022  -8.646  1.00 24.35           N  
ATOM      2  CA  PHE A   1      33.319  19.558  -8.153  1.00 24.35           C  
ATOM      3  C   PHE A   1      32.243  19.483  -9.229  1.00 24.35           C  
ATOM      4  CB  PHE A   1      33.494  21.006  -7.685  1.00 24.35           C  
***/

/***
In PDB files x,y,z coordinates must be -999.999 .. 9999.999
due to fixed-width formatting.
Multiply by 10 for integer representation is
-9999 to 99999 = 199998 range, log2(199998)=17.6, 2^18 = 262144, 3*18 = 54 bits
5 bits for 0..20 amino alphabet
leaves 64-59 = 5 bits unused in 64-bit integer.

Range in PDB: -963.585 .. 3017.753

16-bit unsigned int = 65535 
CoordToIC(X) = int((X + 1000)*10 + 0.5f)
ICToCoord(IC) = float(IC/10.0f) - 1000.0f
C:\src\py\reseek_integer_coords.py
# reseek_integer_coords.py
Min coord -1000.0
Max coord 5553.5
 -999.9           1   -999.9
X=9999.9 IC 109999 out of range
 9999.9      109999   9999.9
-1000.0           0  -1000.0
 -100.0        9000   -100.0
  -50.0        9500    -50.0
    0.0       10000      0.0
   50.0       10500     50.0
  100.0       11000    100.0
 1000.0       20000   1000.0
***/

char GetFeatureChar(byte Letter, uint AlphaSize)
	{
	asserta(AlphaSize <= 36);
	if (Letter == UINT_MAX)
		return '*';
	if (Letter < 26)
		return 'A' + Letter;
	else if (Letter < 36)
		return 'a' + (Letter - 26);
	asserta(false);
	return '!';
	}

void PDBChain::LogMe(bool WithCoords) const
	{
	const size_t L = m_Xs.size();
	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	Log("\n");
	Log(">%s\n", m_Label.c_str());
	Log("%s\n", m_Seq.c_str());

	if (WithCoords)
		{
		Log("\n");
		Log(" Pos  S         X         Y         Z\n");
		for (size_t i = 0; i < L; ++i)
			{
			Log("%4u", i);
			Log("  %c", m_Seq[i]);
			Log("  %8.3f", m_Xs[i]);
			Log("  %8.3f", m_Ys[i]);
			Log("  %8.3f", m_Zs[i]);
			Log("\n");
			}
		}
	}

void PDBChain::ToFeatureFasta(FILE *f, DSS &D, FEATURE Feat) const
	{
	if (f == 0)
		return;
	D.Init(*this);
	const uint L = GetSeqLength();
	const uint AlphaSize = D.GetAlphaSize(Feat);
	string Seq;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		uint Letter = D.GetFeature(Feat, Pos);
		char c = GetFeatureChar(Letter, AlphaSize);
		Seq += c;
		}
	SeqToFasta(f, m_Label, Seq);
	}

void PDBChain::ToFasta(FILE *f) const
	{
	if (f == 0)
		return;
	SeqToFasta(f, m_Label, m_Seq);
	}

void PDBChain::ToFasta(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	ToFasta(f);
	CloseStdioFile(f);
	}

void PDBChain::ToCal(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	ToCal(f);
	CloseStdioFile(f);
	}

void PDBChain::ToCal(FILE* f) const
	{
	if (f == 0)
		return;
	uint QL = GetSeqLength();
	ToCalSeg(f, 0, QL);
	}

void PDBChain::ToCalSeg(FILE *f, uint Pos, uint n) const
	{
	if (f == 0)
		return;
	if (n == 0)
		return;

	const size_t L = m_Xs.size();
	if (L == 0)
		return;

	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	fprintf(f, ">%s\n", m_Label.c_str());

	for (size_t i = Pos; i < Pos + n; ++i)
		{
		if (i < SIZE(m_Seq))
			fputc(m_Seq[i], f);
		else
			fputc('*', f);
		fprintf(f, "\t%.1f", m_Xs[i]);
		fprintf(f, "\t%.1f", m_Ys[i]);
		fprintf(f, "\t%.1f", m_Zs[i]);
		fputc('\n', f);
		}
	}

static float GetFloatFromString(const string &s, uint Pos, uint n)
	{
	string t = s.substr(Pos, n);
	StripWhiteSpace(t);
	float Value = (float) StrToFloat(t);
	return Value;
	}

// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
// 47 - 57        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
void PDBChain::GetXYZFromATOMLine(const string &InputLine,
  float &x, float &y, float &z)
	{
	x = GetFloatFromString(InputLine, 30, 8);
	y = GetFloatFromString(InputLine, 38, 8);
	z = GetFloatFromString(InputLine, 46, 8);
	}

void PDBChain::SetXYZInATOMLine(const string &InputLine,
  float x, float y, float z, string &OutputLine)
	{
	string sx;
	string sy;
	string sz;
	Ps(sx, "%8.3f", x);
	Ps(sy, "%8.3f", y);
	Ps(sz, "%8.3f", z);
	asserta(SIZE(sx) == 8);
	asserta(SIZE(sy) == 8);
	asserta(SIZE(sz) == 8);

	OutputLine = InputLine;
	for (uint i = 0; i < 8; ++i)
		{
		OutputLine[30+i] = sx[i];
		OutputLine[38+i] = sy[i];
		OutputLine[46+i] = sz[i];
		}
	}

bool PDBChain::GetFieldsFromATOMLine(const string &Line,
  float &X, float &Y, float &Z, char &aa)
	{
	aa = 'X';
	X = -999;
	Y = -999;
	Z = -999;
	string AtomName = Line.substr(12, 4);
	StripWhiteSpace(AtomName);
	if (AtomName != "CA")
		return false;
	char AltLoc = Line[16];
	if (AltLoc != ' ' && AltLoc != 'A' && AltLoc != '1')
		return false;

	string AAA = Line.substr(17, 3);
	aa = GetOneFromThree(AAA);

	string sX, sY, sZ;
	sX = Line.substr(30, 8);
	sY = Line.substr(38, 8);
	sZ = Line.substr(46, 8);

	StripWhiteSpace(sX);
	StripWhiteSpace(sY);
	StripWhiteSpace(sZ);

	X = StrToFloatf(sX);
	Y = StrToFloatf(sY);
	Z = StrToFloatf(sZ);

	return true;
	}

void PDBChain::GetPt(uint Pos, vector<float> &Pt) const
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));

	Resize3(Pt);
	Pt[X] = m_Xs[Pos];
	Pt[Y] = m_Ys[Pos];
	Pt[Z] = m_Zs[Pos];
	}

void PDBChain::GetPt(uint Pos, vector<double> &Pt) const
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));

	Resize3(Pt);
	Pt[X] = m_Xs[Pos];
	Pt[Y] = m_Ys[Pos];
	Pt[Z] = m_Zs[Pos];
	}

void PDBChain::SetPt(uint Pos, const vector<float> &Pt)
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));

	m_Xs[Pos] = Pt[X];
	m_Ys[Pos] = Pt[Y];
	m_Zs[Pos] = Pt[Z];
	}

void PDBChain::GetXYZ(uint Pos, float &x, float &y, float &z) const
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));
	x = m_Xs[Pos];
	y = m_Ys[Pos];
	z = m_Zs[Pos];
	}

float PDBChain::GetDist(uint Pos1, uint Pos2) const
	{
	float x1, y1, z1;
	float x2, y2, z2;
	GetXYZ(Pos1, x1, y1, z1);
	GetXYZ(Pos2, x2, y2, z2);
	float d = GetDist3D(x1, y1, z1, x2, y2, z2);
	return d;
	}

float PDBChain::GetDist2(uint Pos1, uint Pos2) const
	{
	float x1 = m_Xs[Pos1];
	float y1 = m_Ys[Pos1];
	float z1 = m_Zs[Pos1];

	float x2 = m_Xs[Pos2];
	float y2 = m_Ys[Pos2];
	float z2 = m_Zs[Pos2];

	float dx = x1 - x2;
	float dy = y1 - y2;
	float dz = z1 - z2;

	float d = GetDist(Pos1, Pos2);

	float d2 = dx*dx + dy*dy + dz*dz;
	asserta(feq(d*d, d2));
	return d2;
	}

uint PDBChain::GetSeqLength() const
	{
	return SIZE(m_Seq);
	}

void PDBChain::GetXFormChain_tR(
  const vector<float> &t,
  const vector<vector<float> > &R,
  PDBChain &XChain) const
	{
	XChain.Clear();
	XChain.m_Label = m_Label;
	XChain.m_Seq = m_Seq;

	const uint N = SIZE(m_Seq);
	asserta(SIZE(m_Xs) == N);
	asserta(SIZE(m_Ys) == N);
	asserta(SIZE(m_Zs) == N);

	vector<float> Pt(3);
	vector<float> XPt(3);
	for (uint Pos = 0; Pos < N; ++Pos)
		{
		GetPt(Pos, Pt);
		XFormPt(Pt, t, R, XPt);

		float x = XPt[X];
		float y = XPt[Y];
		float z = XPt[Z];

		XChain.m_Xs.push_back(x);
		XChain.m_Ys.push_back(y);
		XChain.m_Zs.push_back(z);
		}
	}

void PDBChain::GetXFormChain_tR(
  const vector<double> &t,
  const vector<vector<double> > &R,
  PDBChain &XChain) const
	{
	XChain.Clear();
	XChain.m_Label = m_Label;
	XChain.m_Seq = m_Seq;

	const uint N = SIZE(m_Seq);
	asserta(SIZE(m_Xs) == N);
	asserta(SIZE(m_Ys) == N);
	asserta(SIZE(m_Zs) == N);

	vector<double> Pt(3);
	vector<double> XPt(3);
	for (uint Pos = 0; Pos < N; ++Pos)
		{
		GetPt(Pos, Pt);
		XFormPt(Pt, t, R, XPt);

		float x = (float) XPt[X];
		float y = (float) XPt[Y];
		float z = (float) XPt[Z];

		XChain.m_Xs.push_back(x);
		XChain.m_Ys.push_back(y);
		XChain.m_Zs.push_back(z);
		}
	}

bool PDBChain::IsATOMLine(const string &Line)
	{
	if (SIZE(Line) < 27)
		return false;
// Insertion code is old PDB hack to preserve numbering from
// a reference sequence like NAST or ICTV numbers -- should
// just ignore residue numbering.
	//if (Line[26] != ' ') // insertion code
	//	return false;
	if (strncmp(Line.c_str(), "ATOM  ", 6) == 0)
		return true;
	return false;
	}

// flattened as x0, y0, z0, x1, y1, z1 ...
void PDBChain::GetICs(vector<uint16_t> &ICs) const
	{
	const uint L = SIZE(m_Xs);
	asserta(SIZE(m_Ys) == L);
	asserta(SIZE(m_Zs) == L);
	ICs.clear();
	ICs.reserve(3*L);
	for (uint i = 0; i < L; ++i)
		{
		ICs.push_back(CoordToIC(m_Xs[i]));
		ICs.push_back(CoordToIC(m_Ys[i]));
		ICs.push_back(CoordToIC(m_Zs[i]));
		}
	}

void PDBChain::CoordsFromICs(const uint16_t *ICs, uint L)
	{
	m_Xs.clear();
	m_Ys.clear();
	m_Zs.clear();
	m_Xs.reserve(L);
	m_Ys.reserve(L);
	m_Zs.reserve(L);
	for (uint i = 0; i < L; ++i)
		{
		m_Xs.push_back(ICToCoord(ICs[3*i]));
		m_Ys.push_back(ICToCoord(ICs[3*i+1]));
		m_Zs.push_back(ICToCoord(ICs[3*i+2]));
		}
	}

void PDBChain::CoordsFromICs(const vector<uint16_t> &ICs)
	{
	const uint L3 = SIZE(ICs);
	asserta(L3%3 == 0);
	const uint L = L3/3;
	m_Xs.clear();
	m_Ys.clear();
	m_Zs.clear();
	m_Xs.reserve(L);
	m_Ys.reserve(L);
	m_Zs.reserve(L);
	for (uint i = 0; i < L; ++i)
		{
		m_Xs.push_back(ICToCoord(ICs[3*i]));
		m_Ys.push_back(ICToCoord(ICs[3*i+1]));
		m_Zs.push_back(ICToCoord(ICs[3*i+2]));
		}
	}

void PDBChain::Reverse()
	{
	std::reverse(m_Seq.begin(), m_Seq.end());
	std::reverse(m_Xs.begin(), m_Xs.end());
	std::reverse(m_Ys.begin(), m_Ys.end());
	std::reverse(m_Zs.begin(), m_Zs.end());
	}

void PDBChain::GetReverse(PDBChain &Rev) const
	{
	Rev = *this;
	Rev.Reverse();
	Rev.m_Label += ".rev";
	}

void PDBChain::Flip()
	{
	const uint L = GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		m_Xs[i] = -m_Xs[i];
		m_Ys[i] = -m_Ys[i];
		m_Zs[i] = -m_Zs[i];
		}
	}
