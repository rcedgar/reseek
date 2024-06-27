#include "myutils.h"
#include "pdbchain.h"
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
		fprintf(f, "\t%.3f", m_Xs[i]);
		fprintf(f, "\t%.3f", m_Ys[i]);
		fprintf(f, "\t%.3f", m_Zs[i]);
		fputc('\n', f);
		}
	}

static double GetFloatFromString(const string &s, uint Pos, uint n)
	{
	string t = s.substr(Pos, n);
	StripWhiteSpace(t);
	double Value = StrToFloat(t);
	return Value;
	}

// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
// 47 - 57        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
void PDBChain::GetXYZFromATOMLine(const string &InputLine,
  double &x, double &y, double &z)
	{
	x = GetFloatFromString(InputLine, 30, 8);
	y = GetFloatFromString(InputLine, 38, 8);
	z = GetFloatFromString(InputLine, 46, 8);
	}

void PDBChain::SetXYZInATOMLine(const string &InputLine,
  double x, double y, double z, string &OutputLine)
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

// 77 - 78        LString(2)    element      Element symbol, right-justified.
void PDBChain::GetElementNameFromATOMLine(const string &Line,
  string &ElementName)
	{
	asserta(SIZE(Line) >= 78);
	ElementName = Line.substr(76, 2);
	StripWhiteSpace(ElementName);
	}

void PDBChain::GetAtomNameFromATOMLine(const string &Line,
  string &AtomName)
	{
	asserta(SIZE(Line) > 40);
	AtomName = Line.substr(12, 4);
	StripWhiteSpace(AtomName);
	}

// Residue nr in Cols 23-26 (1-based)
int PDBChain::GetResidueNrFromATOMLine(const string &Line)
	{
	asserta(SIZE(Line) > 60);
	string s;
	s += Line[22];
	s += Line[23];
	s += Line[24];
	s += Line[25];
	StripWhiteSpace(s);
	int Nr = atoi(s.c_str());
	return Nr;
	}

// Atom nr in Cols 7-11 (1-based)
//   7 - 11        Integer       serial       Atom  serial number.
void PDBChain::SetAtomNrInATOMLine(const string &InputLine,
  uint AtomNr, string &OutputLine)
	{
	string s;
	Ps(s, "%5u", AtomNr);
	uint n = SIZE(s);
	if (n != 5)
		Die("Residue number %d overflow (max 9999)", AtomNr);
	OutputLine.clear();
	uint N = SIZE(InputLine);
	for (uint i = 0; i < 6; ++i)
		OutputLine += InputLine[i];
	for (uint i = 0; i < 5; ++i)
		OutputLine += s[i];
	for (uint i = 6+5; i < N; ++i)
		OutputLine += InputLine[i];
	asserta(SIZE(OutputLine) == SIZE(InputLine));
	}

// Residue nr in Cols 23-26 (1-based)
void PDBChain::SetResidueNrInATOMLine(const string &InputLine,
  uint ResidueNr, string &OutputLine)
	{
	string s;
	Ps(s, "%4u", ResidueNr);
	uint n = SIZE(s);
	if (n != 4)
		Die("Residue number %d overflow (max 9999)", ResidueNr);
	OutputLine.clear();
	uint N = SIZE(InputLine);
	for (uint i = 0; i < 22; ++i)
		OutputLine += InputLine[i];
	for (uint i = 0; i < 4; ++i)
		OutputLine += s[i];
	for (uint i = 22+4; i < N; ++i)
		OutputLine += InputLine[i];
	asserta(SIZE(OutputLine) == SIZE(InputLine));
	}

void PDBChain::HackHETAMLine(string &Line)
	{
	if (strncmp(Line.c_str(), "HETATM", 6) == 0 &&
		Line.substr(17, 3) == "MSE")
		{
		Line[0] = 'A';
		Line[1] = 'T';
		Line[2] = 'O';
		Line[3] = 'M';
		Line[4] = ' ';
		Line[5] = ' ';

		Line[17] = 'M';
		Line[18] = 'E';
		Line[19] = 'T';
		}
	}

char PDBChain::FromPDBLines(const string &Label,
  const vector<string> &Lines)
	{
	Clear();
	m_Label = Label;
	const uint N = SIZE(Lines);
	char ChainChar = 0;
	uint ResidueCount = 0;
	int CurrentResidueNumber = INT_MAX;
	vector<string> ATOMLines;
	for (uint LineNr = 0; LineNr < N; ++LineNr)
		{
		const string &Line = Lines[LineNr];
		const size_t L = Line.size();

		char LineChain = Line[21];
		if (ChainChar == 0)
			ChainChar = LineChain;
		else if (ChainChar != LineChain)
			Die("PDBChain::FromPDBLines() two chains %c, %c",
			  ChainChar, LineChain);

		char aa;
		double X, Y, Z;
		bool IsCA = GetFieldsFromATOMLine(Line, X, Y, Z, aa);
		if (!IsCA)
			continue;

		m_Seq.push_back(aa);
		m_Xs.push_back(X);
		m_Ys.push_back(Y);
		m_Zs.push_back(Z);
		}

	if (ChainChar != 0 && !isspace(ChainChar))
		{
		m_Label += optset_chainsep ? string(opt_chainsep) : ":";
		m_Label += ChainChar;
		}

	return ChainChar;
	}

bool PDBChain::GetFieldsFromATOMLine(const string &Line,
  double &X, double &Y, double &Z, char &aa)
	{
	aa = 'X';
	X = -999;
	Y = -999;
	Z = -999;
	string AtomName = Line.substr(12, 4);
	StripWhiteSpace(AtomName);
	if (AtomName != "CA")
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

	X = StrToFloat(sX);
	Y = StrToFloat(sY);
	Z = StrToFloat(sZ);

	return true;
	}

bool PDBChain::GetFieldsFromResidueATOMLines(const vector<string> &Lines,
  double &X, double &Y, double &Z, char &aa, int &ResNr)
	{
	aa = 'X';
	ResNr = -999;
	X = -999;
	Y = -999;
	Z = -999;
	const uint N = SIZE(Lines);
	for (uint i = 0; i < SIZE(Lines); ++i)
		{
		const string &Line = Lines[i];
		string AtomName = Line.substr(12, 4);
		StripWhiteSpace(AtomName);
		if (AtomName == "CA")
			{
			string AAA = Line.substr(17, 3);
			aa = GetOneFromThree(AAA);

			string sX, sY, sZ;
			sX = Line.substr(30, 8);
			sY = Line.substr(38, 8);
			sZ = Line.substr(46, 8);

			StripWhiteSpace(sX);
			StripWhiteSpace(sY);
			StripWhiteSpace(sZ);

			X = StrToFloat(sX);
			Y = StrToFloat(sY);
			Z = StrToFloat(sZ);

			ResNr = GetResidueNrFromATOMLine(Line);
			return true;
			}
		}
	return false;
	}

void PDBChain::GetSubSeq(uint Pos, uint n, string &s) const
	{
	if (Pos == UINT_MAX)
		{
		s = ".";
		return;
		}

	s.clear();
	size_t L = m_Seq.size();
	asserta(Pos + n <= L);

	for (uint i = 0; i < n; ++i)
		{
		char c = m_Seq[Pos+i];
		s += c;
		}
	}

double PDBChain::GetCoord(uint Axis, uint Pos) const
	{
	switch (Axis)
		{
	case X: return m_Xs[Pos];
	case Y: return m_Ys[Pos];
	case Z: return m_Zs[Pos];
		}
	asserta(false);
	return 0;
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

void PDBChain::SetPt(uint Pos, const vector<double> &Pt)
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));

	m_Xs[Pos] = Pt[X];
	m_Ys[Pos] = Pt[Y];
	m_Zs[Pos] = Pt[Z];
	}

double PDBChain::GetX(uint Pos) const
	{
	assert(Pos < SIZE(m_Xs));
	return m_Xs[Pos];
	}

double PDBChain::GetY(uint Pos) const
	{
	assert(Pos < SIZE(m_Ys));
	return m_Ys[Pos];
	}

double PDBChain::GetZ(uint Pos) const
	{
	assert(Pos < SIZE(m_Zs));
	return m_Zs[Pos];
	}

void PDBChain::GetXYZ(uint Pos, double &x, double &y, double &z) const
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));
	x = m_Xs[Pos];
	y = m_Ys[Pos];
	z = m_Zs[Pos];
	}

void PDBChain::GetDistMx(uint Pos, uint L,
  vector<vector<double> > &Mx) const
	{
	Mx.resize(L);
	for (uint i = 0; i < L; ++i)
		{
		Mx[i].resize(L);
		for (uint j = 0; j < L; ++j)
			{
			double d = GetDist(Pos+i, Pos+j);
			Mx[i][j] = d;
			}
		}
	}

double PDBChain::GetDist(uint Pos1, uint Pos2) const
	{
	double x1, y1, z1;
	double x2, y2, z2;
	GetXYZ(Pos1, x1, y1, z1);
	GetXYZ(Pos2, x2, y2, z2);
	double d = GetDist3D(x1, y1, z1, x2, y2, z2);
	return d;
	}

double PDBChain::GetDist2(uint Pos1, uint Pos2) const
	{
	double x1 = m_Xs[Pos1];
	double y1 = m_Ys[Pos1];
	double z1 = m_Zs[Pos1];

	double x2 = m_Xs[Pos2];
	double y2 = m_Ys[Pos2];
	double z2 = m_Zs[Pos2];

	double dx = x1 - x2;
	double dy = y1 - y2;
	double dz = z1 - z2;

	double d = GetDist(Pos1, Pos2);

	double d2 = dx*dx + dy*dy + dz*dz;
	asserta(feq(d*d, d2));
	return d2;
	}

void PDBChain::WriteSeqWithCoords(FILE *f) const
	{
	if (f == 0)
		return;
	const uint L = GetSeqLength();
	uint StartPos = 0;
	for (;;)
		{
		uint EndPos = StartPos + 99;
		if (EndPos >= L)
			EndPos = L - 1;

		fprintf(f, "\n");
		for (uint Pos = StartPos; Pos <= EndPos; Pos += 10)
			fprintf(f, "%-10u", Pos);
		fprintf(f, "\n");
		for (uint Pos = StartPos; Pos <= EndPos; ++Pos)
			{
			if (Pos%10 == 0)
				fprintf(f, " ");
			else
				fprintf(f, "%u", Pos%10);
			}
		fprintf(f, "\n");
		for (uint Pos = StartPos; Pos <= EndPos; ++Pos)
			fprintf(f, "%c", m_Seq[Pos]);
		fprintf(f, "\n");

		StartPos = EndPos + 1;
		if (StartPos >= L)
		break;
		}
	fprintf(f, "\n");
	}

void PDBChain::PrintSeqCoords(FILE *f) const
	{
	if (f == 0)
		return;
	const uint L = GetSeqLength();

	fprintf(f, "\n");
	fprintf(f, ">%s\n", m_Label.c_str());
	for (uint i = 0; i < L; ++i)
		{
		if ((i+1)%10 == 0)
			fprintf(f, "%10u", (i+1)/10);
		}
	fprintf(f, "\n");

	for (uint i = 0; i < L; ++i)
		fprintf(f, "%u", (i+1)%10);
	fprintf(f, "\n");

	fprintf(f, "%s\n", m_Seq.c_str());
	}

uint PDBChain::GetSeqLength() const
	{
	return SIZE(m_Seq);
	}

void PDBChain::GetSubSeq(uint MotifStartPos, uint n,
  bool FailOnOverflow, string &MotifSeq) const
	{
	MotifSeq.clear();
	const uint L = GetSeqLength();
	for (uint i = 0; i < n; ++i)
		{
		uint Pos = MotifStartPos + i;
		if (Pos < 0 || Pos >= L)
			{
			if (FailOnOverflow)
				Die("'%s' GetMotifSeqFromMidPos overflow",
				  m_Label.c_str());
			MotifSeq += '!';
			}
		else
			MotifSeq += m_Seq[Pos];
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

		double x = XPt[X];
		double y = XPt[Y];
		double z = XPt[Z];

		XChain.m_Xs.push_back(x);
		XChain.m_Ys.push_back(y);
		XChain.m_Zs.push_back(z);
		}
	}

bool PDBChain::IsATOMLine(const string &Line)
	{
	if (SIZE(Line) < 27)
		return false;
	if (Line[26] != ' ') // insertion code
		return false;
	if (strncmp(Line.c_str(), "ATOM  ", 6) == 0)
		return true;
	return false;
	}

char PDBChain::GetChainCharFromATOMLine(const string &Line)
	{
	if (!IsATOMLine(Line))
		return 0;
	if (Line.size() < 22)
		return 0;
	return Line[21];
	}

double PDBChain::GetSmoothedCoord(uint Axis, uint i, uint N, uint w) const
	{
	double SumCoord = 0;
	double SumWeight = 0;
	double dN = (double) N;
	double dw = (double) w;
	int iQL = (int) GetSeqLength();
	double dQL = (double) iQL;
	for (int k = -1; k <= int(w); ++k)
		{
		double dPos = (int(i) + k)*dQL/dN;
		int iPos = int(dPos + 0.5);
		if (iPos < 0 || iPos >= iQL)
			continue;
		double Weight = dw - fabs(dPos - iPos);
		double Coord = GetCoord(Axis, (uint) iPos);
		SumCoord += Weight*Coord;
		SumWeight += Weight;
		}
	asserta(SumWeight > 0);
	double SmoothedCoord = SumCoord/SumWeight;
	return SmoothedCoord;
	}

static uint GetDiffs3(const char *s, const char *t)
	{
	uint n = 0;
	for (uint i = 0; i < 3; ++i)
		if (s[i] != t[i])
			++n;
	return n;
	}

const char *PDBChain::GetAcc(string &Acc) const
	{
	size_t n = m_Label.find(' ');
	if (n > 0)
		Acc = m_Label.substr(0, n);
	else
		Acc = m_Label;
	return Acc.c_str();
	}

void PDBChain::GetRange(uint Lo, uint Hi, PDBChain &Chain) const
	{
	Chain.Clear();
	asserta(Lo <= Hi);
	uint L = GetSeqLength();
	asserta(Hi < L);

	Chain.m_Label = m_Label;

	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
#define c(x)	Chain.m_##x.push_back(m_##x[Pos])
		c(Seq);
		c(Xs);
		c(Ys);
		c(Zs);
#undef c
		}
	}

void PDBChain::GetSphere(uint Pos, double Radius,
  uint MinPos, uint MaxPos, vector<uint> &PosVec) const
	{
	PosVec.clear();
	const uint L = GetSeqLength();
	asserta(MinPos <= MaxPos && MaxPos < L);
	for (uint Pos2 = MinPos; Pos2 <= MaxPos; ++Pos2)
		{
		if (Pos2 == Pos)
			continue;
		double d = GetDist(Pos, Pos2);
		if (d <= Radius)
			PosVec.push_back(Pos2);
		}
	}
