#include "myutils.h"
#include "paralign.h"
#include "nu.h"
#include "seqdb.h"
#include "alpha.h"
#include "cigar.h"

void ExpandParaCigar_reverseDI(const string &s, string &Path);
void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);
float SWFast_SubstMx(XDPMem &Mem,
	const byte *A, uint LA, const byte *B, uint LB,
	const vector<vector<float> > &SubstMx,
	float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
	string &Path);

extern parasail_matrix_t parasail_mu_matrix;
extern int Blosum62_int[20][20];
extern int Mu_S_k_i8[36*36];
extern int Mu_hjmumx[36*36];
extern int parasail_mu_8[36*36];

parasail_matrix_t Paralign::m_matrix;
int Paralign::m_Open = INT_MAX;	// penalty > 0
int Paralign::m_Ext = INT_MAX;	// penalty > 0
int Paralign::m_SaturatedScore = INT_MAX;
int Paralign::m_Bits = 16;
string Paralign::m_SubstMxName = "_NOT_SET_";
vector<vector<float> > Paralign::m_SWFastSubstMx;
atomic<uint> Paralign::m_Count8;
atomic<uint> Paralign::m_Count16;
atomic<uint> Paralign::m_CountSWFast;
atomic<uint> Paralign::m_SaturatedCount;

bool Paralign::m_GapLengthDist = false;
uint Paralign::m_MaxGapLength = UINT_MAX;
vector<uint> Paralign::m_GapLengthToCount;
omp_lock_t Paralign::m_GapLengthLock;
     
// Ye olde BLOSUM62 as used by NCBI BLAST (1/2-bit units)
// alphabetical order, no wildcards or stop codon
static const int Blosum62_Open = 11;
static const int Blosum62_Ext = 1;
static const int Blosum62_SaturatedScore = 999;

void Paralign::LogGapLengthDist()
	{
	ProgressLog("Paralign::LogGapLengthDist() max=%u\n", m_MaxGapLength);
	for (uint L = 1; L <= m_MaxGapLength; ++L)
		ProgressLog("%3u  %u\n", L, m_GapLengthToCount[L]);
	}

void Paralign::InitGapLengthDist(uint MaxLen)
	{
	asserta(m_MaxGapLength == UINT_MAX);
	m_MaxGapLength = MaxLen;
	m_GapLengthToCount.clear();
	m_GapLengthToCount.resize(m_MaxGapLength+1);
	omp_init_lock(&m_GapLengthLock);
	m_GapLengthDist = true;
	}

void Paralign::SetSWFastSubstMx_FromParasailMx()
	{
	const uint AS = GetAlphaSize();
	m_SWFastSubstMx.resize(AS);
	for (uint i = 0; i < AS; ++i)
		{
		m_SWFastSubstMx[i].resize(AS);
		for (uint j = 0; j < AS; ++j)
			m_SWFastSubstMx[i][j] = (float) GetSubstScore(i, j);
		}
	}

void Paralign::SetSWFastSubstMx(const vector<vector<float> > &Mx,
	int Open, int Ext, bool DisableParasail)
	{
	m_Open = Open;
	m_Ext = Ext;
	const uint AS = SIZE(Mx[0]);
	m_matrix.size = AS;
	m_SWFastSubstMx.resize(AS);
	for (uint i = 0; i < AS; ++i)
		{
		m_SWFastSubstMx[i].resize(AS);
		for (uint j = 0; j < AS; ++j)
			m_SWFastSubstMx[i][j] = Mx[i][j];
		}
	if (DisableParasail)
		memset(&m_matrix, 0, sizeof(m_matrix));
	}

void Paralign::UpdateGapLengthDist(const string &Path)
	{
	asserta(SIZE(m_GapLengthToCount) == m_MaxGapLength + 1);
	uint L = 0;
	const uint ColCount = SIZE(Path);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M')
			{
			if (L != 0)
				{
				omp_set_lock(&m_GapLengthLock);
				m_GapLengthToCount[min(L, m_MaxGapLength)] += 1;
				omp_unset_lock(&m_GapLengthLock);
				}
			L = 0;
			continue;
			}
		else
			++L;
		}
	if (L != 0)
		{
		omp_set_lock(&m_GapLengthLock);
		m_GapLengthToCount[min(L, m_MaxGapLength)] += 1;
		omp_unset_lock(&m_GapLengthLock);
		}
	}

void Paralign::Align_SWFast(const string &LabelT, const byte *T, uint LT)
	{
	ClearResult();

	const uint AS = GetAlphaSize();
	m_LabelT = LabelT;
	m_T = T;
	m_LT = LT;

	float Open = -float(m_Open);
	float Ext = -float(m_Ext);

	uint LoQ, LenQ, LoT, LenT;
	m_SWFastScore = SWFast_SubstMx(m_Mem, m_Q, m_LQ, m_T, m_LT,
		m_SWFastSubstMx, Open, Ext, LoQ, LenQ, LoT, LenT, m_SWFastPath);
	if (m_GapLengthDist)
		UpdateGapLengthDist(m_SWFastPath);
	m_SWFastScoreInt = int(round(m_SWFastScore));
	++m_CountSWFast;
	}

int Paralign::GetSubstScore(uint LetterQ, uint LetterT)
	{
	const uint AS = GetAlphaSize();
	asserta(LetterQ < AS);
	asserta(LetterT < AS);
	return m_matrix.matrix[LetterQ*AS + LetterT];
	}

const char *Paralign::GetLetterToChar() const
	{
	return m_matrix.alphabet;
	}

void Paralign::WriteAln(FILE *f) const
	{
	if (f == 0)
		return;
	if (m_SemiGlobalPath.empty())
		return;

	const char *LetterToChar = GetLetterToChar();
	const uint AS = GetAlphaSize();

	string RowQ;
	string RowT;
	string Annot;
	uint PosQ = m_LoQ;
	uint PosT = m_LoT;
	int Score = 0;
	const uint ColCount = SIZE(m_SemiGlobalPath);
	bool InGap = false;
	const int *mx = Paralign::m_matrix.matrix;
	uint FirstM = UINT_MAX;
	uint LastM = UINT_MAX;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = m_SemiGlobalPath[Col];
		if (c == 'M')
			{
			if (FirstM == UINT_MAX)
				FirstM = Col;
			LastM = Col;
			}
		}
	for (uint Col = 0; Col < FirstM; ++Col)
		{
		if (m_SemiGlobalPath[Col] == 'D')
			++PosQ;
		else if (m_SemiGlobalPath[Col] == 'I')
			++PosT;
		else
			asserta(false);
		}

	for (uint Col = FirstM; Col <= LastM; ++Col)
		{
		char c = m_SemiGlobalPath[Col];
		switch (c)
			{
		case 'M':
			{
			asserta(PosQ < m_LQ);
			asserta(PosT < m_LT);
			byte iq = m_Q[PosQ];
			byte it = m_T[PosT];
			asserta(iq < AS);
			asserta(it < AS);
			char cq = LetterToChar[iq];
			char ct = LetterToChar[it];
			RowQ += cq;
			RowT += ct;
			if (cq == ct)
				Annot += '|';
			else
				{
				int s = GetSubstScore(iq, it);
				if (s > 0)
					Annot += '+';
				else
					Annot += ' ';
				}
			++PosQ;
			++PosT;
			break;
			}

		case 'D':
			{
			asserta(PosQ < m_LQ);
			byte iq = m_Q[PosQ];
			asserta(iq < AS);
			char cq = LetterToChar[iq];
			RowQ += cq;
			RowT += '-';
			Annot += ' ';
			++PosQ;
			break;
			}

		case 'I':
			{
			asserta(PosT < m_LT);
			byte it = m_T[PosT];
			asserta(it < AS);
			char ct = LetterToChar[it];
			RowT += ct;
			RowQ += '-';
			Annot += ' ';
			++PosT;
			break;
			}
		default:
			asserta(false);
			}
		}
	fprintf(f, "\n");
	fprintf(f, "%s  %s\n", RowQ.c_str(), m_LabelQ.c_str());
	fprintf(f, "%s\n", Annot.c_str());
	fprintf(f, "%s  %s\n", RowT.c_str(), m_LabelT.c_str());
	fprintf(f, "  Score %d\n", m_Score);
	}

int Paralign::ScoreAln(bool Trace) const
	{
	if (m_SemiGlobalPath.empty())
		return 0;
	uint PosQ = m_LoQ;
	uint PosT = m_LoT;
	int Score = 0;
	const uint ColCount = SIZE(m_SemiGlobalPath);
	bool InGap = false;
	const int *mx = Paralign::m_matrix.matrix;
	uint FirstM = UINT_MAX;
	uint LastM = UINT_MAX;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = m_SemiGlobalPath[Col];
		if (c == 'M')
			{
			if (FirstM == UINT_MAX)
				FirstM = Col;
			LastM = Col;
			}
		}
	for (uint Col = 0; Col < FirstM; ++Col)
		{
		if (m_SemiGlobalPath[Col] == 'D')
			++PosQ;
		else if (m_SemiGlobalPath[Col] == 'I')
			++PosT;
		else
			asserta(false);
		}

	for (uint Col = FirstM; Col <= LastM; ++Col)
		{
		char c = m_SemiGlobalPath[Col];
		if (c == 'M')
			{
			uint lettera = m_Q[PosQ];
			uint letterb = m_T[PosT];
			InGap = false;
			int SubstScore = mx[lettera*36 + letterb];
			Score += SubstScore;
			if (Trace)
				{
				Log("M PosQ=%u PosT=%u", PosQ, PosT);
				Log(" %u/%c", lettera, g_LetterToCharMu[lettera]);
				Log(" %u/%c", letterb, g_LetterToCharMu[letterb]);
				Log(" %+d", SubstScore);
				Log(" %d\n", Score);
				}
			++PosQ;
			++PosT;
			}
		else if (c == 'D')
			{
			if (InGap)
				{
				if (Trace)
					{
					Log("D -%d\n", m_Ext);
					Log(" %d\n", Score);
					}
				Score -= m_Ext;
				}
			else
				{
				if (Trace)
					{
					Log("D -%d\n", m_Open);
					Log(" %d\n", Score);
					}
				InGap = true;
				Score -= m_Open;
				}
			++PosQ;
			}
		else if (c == 'I')
			{
			if (InGap)
				{
				if (Trace)
					{
					Log("I -%d\n", m_Ext);
					Log(" %d\n", Score);
					}
				Score -= m_Ext;
				}
			else
				{
				if (Trace)
					{
					Log("I -%d\n", m_Open);
					Log(" %d\n", Score);
					}
				InGap = true;
				Score -= m_Open;
				}
			++PosT;
			}
		else
			asserta(false);
		}
	return Score;
	}

/***
2025-10_reseek_tune/2025-11-27_Mu_S_k_i8_gap_sweep
open30.ext4  1.259
open20.ext8  1.264
open24.ext6  1.265
open28.ext8  1.265
open29.ext5  1.266
open24.ext9  1.267
open30.ext7  1.267
open23.ext7  1.268
open30.ext5  1.268
open25.ext8  1.269
open29.ext7  1.269
open23.ext8  1.270
open24.ext7  1.270
open24.ext8  1.270
open28.ext6  1.270
open30.ext6  1.270
open29.ext6  1.271 <<<
***/
void Paralign::Set_Mu_S_k_i8()
	{
	m_Open = 29;
	m_Ext = 6;
	if (optset_intopen)
		m_Open = opt(intopen);
	if (optset_intext)
		m_Ext = opt(intext);

	int MinScore = 0;
	int MaxScore = 0;
	for (uint i = 0; i < 36*36; ++i)
		{
		int Score = Mu_S_k_i8[i];
		if (i == 0 || Score < MinScore) MinScore = Score;
		if (i == 0 || Score > MaxScore) MaxScore = Score;
		}
	m_matrix.size = 36;
	m_matrix.length = 36;
	m_matrix.type = PARASAIL_MATRIX_TYPE_SQUARE;
	m_matrix.matrix = Mu_S_k_i8;
	m_matrix.min = MinScore;
	m_matrix.max = MaxScore;
	int *Mapper = myalloc(int, 256);
	memset(Mapper, 0, 256*sizeof(int));
	for (int i = 0; i < 36; ++i)
		Mapper[i] = i;
	m_matrix.mapper = Mapper;

	m_Bits = 16;
	SetSWFastSubstMx_FromParasailMx();
	}

void Paralign::SetMu_parasail_mu_8()
	{
	m_Open = 2;
	m_Ext = 1;
	if (optset_intopen)
		m_Open = opt(intopen);
	if (optset_intext)
		m_Ext = opt(intext);

	int MinScore = 0;
	int MaxScore = 0;
	for (uint i = 0; i < 36*36; ++i)
		{
		int Score = Mu_S_k_i8[i];
		if (i == 0 || Score < MinScore) MinScore = Score;
		if (i == 0 || Score > MaxScore) MaxScore = Score;
		}
	m_matrix.size = 36;
	m_matrix.length = 36;
	m_matrix.type = PARASAIL_MATRIX_TYPE_SQUARE;
	m_matrix.matrix = parasail_mu_8;
	m_matrix.min = MinScore;
	m_matrix.max = MaxScore;
	int *Mapper = myalloc(int, 256);
	memset(Mapper, 0, 256*sizeof(int));
	for (int i = 0; i < 36; ++i)
		Mapper[i] = i;
	m_matrix.mapper = Mapper;

	m_Bits = 8;
	SetSWFastSubstMx_FromParasailMx();
	}

void Paralign::Set_Mu_hjmux()
	{
	m_Open = 4;
	m_Ext = 2;
	if (optset_intopen)
		m_Open = opt(intopen);
	if (optset_intext)
		m_Ext = opt(intext);

	int MinScore = 0;
	int MaxScore = 0;
	for (uint i = 0; i < 36*36; ++i)
		{
		int Score = Mu_hjmumx[i];
		if (i == 0 || Score < MinScore) MinScore = Score;
		if (i == 0 || Score > MaxScore) MaxScore = Score;
		}
	m_matrix.size = 36;
	m_matrix.length = 36;
	m_matrix.type = PARASAIL_MATRIX_TYPE_SQUARE;
	m_matrix.matrix = Mu_hjmumx;
	m_matrix.min = MinScore;
	m_matrix.max = MaxScore;
	int *Mapper = myalloc(int, 256);
	memset(Mapper, 0, 256*sizeof(int));
	for (int i = 0; i < 36; ++i)
		Mapper[i] = i;
	m_matrix.mapper = Mapper;

	m_Bits = 16;
	SetSWFastSubstMx_FromParasailMx();
	}

/***
2025-10_reseek_tune/2025-11-27_musubstmx_gap_sweep

open3.ext2.log:SEPQ0.1=0.129 SEPQ1=0.216 SEPQ10=0.395 Area0=0.549 Sum3=0.977
open4.ext2.log:SEPQ0.1=0.124 SEPQ1=0.216 SEPQ10=0.391 Area0=0.558 Sum3=0.963
open3.ext3.log:SEPQ0.1=0.123 SEPQ1=0.217 SEPQ10=0.388 Area0=0.560 Sum3=0.959
open4.ext1.log:SEPQ0.1=0.127 SEPQ1=0.205 SEPQ10=0.390 Area0=0.536 Sum3=0.952
open5.ext1.log:SEPQ0.1=0.126 SEPQ1=0.206 SEPQ10=0.389 Area0=0.530 Sum3=0.949
open3.ext1.log:SEPQ0.1=0.128 SEPQ1=0.202 SEPQ10=0.387 Area0=0.531 Sum3=0.947
open5.ext2.log:SEPQ0.1=0.119 SEPQ1=0.214 SEPQ10=0.384 Area0=0.548 Sum3=0.943
open2.ext1.log:SEPQ0.1=0.128 SEPQ1=0.196 SEPQ10=0.377 Area0=0.518 Sum3=0.927
open6.ext2.log:SEPQ0.1=0.116 SEPQ1=0.210 SEPQ10=0.377 Area0=0.536 Sum3=0.925
open8.ext1.log:SEPQ0.1=0.115 SEPQ1=0.204 SEPQ10=0.379 Area0=0.527 Sum3=0.916
open4.ext4.log:SEPQ0.1=0.113 SEPQ1=0.208 SEPQ10=0.370 Area0=0.527 Sum3=0.907
open8.ext2.log:SEPQ0.1=0.108 SEPQ1=0.203 SEPQ10=0.367 Area0=0.514 Sum3=0.888
open8.ext4.log:SEPQ0.1=0.095 SEPQ1=0.185 SEPQ10=0.341 Area0=0.463 Sum3=0.806
***/
void Paralign::SetMu_musubstmx()
	{
	int Open = 3;
	int Ext = 2;
	if (optset_intopen)
		Open = opt(intopen);
	if (optset_intext)
		Ext = opt(intext);
	extern float musubstmx[36][36];
	vector<vector<float> > ScoreMx(36);
	for (uint i = 0; i < 36; ++i)
		{
		ScoreMx[i].resize(36);
		for (uint j = 0; j < 36; ++j)
			ScoreMx[i][j] = musubstmx[i][j];
		}

	/////////////////////////////////////////////////
	// Supports SWFast only, disables parasail matrix
	/////////////////////////////////////////////////
	SetSWFastSubstMx(ScoreMx, Open, Ext, false);
	}

// Low accuracy 
// C:\src\notebooks\2025-11-24_mu_variants.txt
// scop40.tm0.6_0.8:TrainedMx R=0.9478 << Sum3=0.54 (subset a)
//   with void Paralign::SetMu_scop40_tm0_6_0_8_fa2()
static int scop40_tm0_6_0_8_fa2[36*36] = {
//    0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26    27    28    29    30    31    32    33    34    35 
     10,   -5,    2,    5,  -15,   -7,    6,   -9,   -4,    4,   -7,   -2,   -1,  -16,   -9,    1,  -11,   -6,   -2,   -8,   -5,   -5,  -21,  -11,   -4,  -12,   -8,   -5,  -10,   -8,   -8,  -19,  -15,   -7,  -14,  -12 ,
     -5,   16,    7,   -7,    4,    1,   -4,    6,    2,   -8,    8,    0,  -11,   -1,   -3,   -7,    2,   -3,   -6,    4,   -2,  -12,   -2,   -4,  -11,    0,   -4,  -10,   -7,   -5,  -16,   -6,   -9,  -12,   -4,   -6 ,
      2,    7,   11,   -1,   -4,    3,    1,    0,    5,   -2,    1,    3,   -5,   -8,   -2,   -3,   -3,    0,   -5,   -3,   -1,   -6,   -7,   -4,   -6,   -5,   -3,   -8,   -6,   -4,  -11,  -12,   -7,   -9,   -9,   -7 ,
      5,   -7,   -1,   13,  -11,    0,    8,   -8,   -2,   -1,   -9,   -4,    7,  -11,   -3,    2,  -10,   -5,   -5,  -10,   -8,    2,  -12,   -6,   -3,  -10,   -7,   -8,  -17,  -11,   -3,  -17,  -10,   -7,  -18,  -12 ,
    -15,    4,   -4,  -11,   10,    3,  -12,    7,   -1,  -18,   -1,  -11,  -13,    3,   -3,  -15,    0,   -7,  -17,   -3,  -12,  -12,    0,   -6,  -16,   -2,   -9,  -23,   -6,  -15,  -19,   -2,   -9,  -18,   -6,  -13 ,
     -7,    1,    3,    0,    3,   10,   -3,    3,    6,  -11,   -5,   -4,   -5,   -1,    3,   -7,   -2,   -1,  -12,   -5,   -7,   -9,   -4,   -2,  -10,   -4,   -4,  -15,   -8,  -10,  -13,   -7,   -5,  -14,   -9,   -8 ,
      6,   -4,    1,    8,  -12,   -3,   10,   -6,    0,    1,   -7,   -3,    3,  -13,   -5,    4,   -7,   -3,   -3,   -8,   -5,   -1,  -13,   -7,   -1,   -9,   -5,   -7,  -12,   -9,   -5,  -16,  -11,   -4,  -12,   -9 ,
     -9,    6,    0,   -8,    7,    3,   -6,   11,    4,  -14,    0,   -6,  -11,    0,   -3,  -10,    3,   -2,  -13,   -2,   -8,  -12,   -2,   -5,  -11,    0,   -5,  -18,   -7,  -13,  -17,   -4,   -9,  -14,   -3,   -8 ,
     -4,    2,    5,   -2,   -1,    6,    0,    4,    8,   -8,   -4,   -2,   -6,   -6,    0,   -4,   -1,    1,  -10,   -6,   -5,  -10,   -7,   -4,   -7,   -4,   -2,  -13,   -9,   -8,  -12,  -11,   -7,  -10,   -7,   -6 ,
      4,   -8,   -2,   -1,  -18,  -11,    1,  -14,   -8,    7,   -4,    0,    0,  -16,   -8,    3,  -10,   -5,    3,   -6,   -2,   -2,  -14,   -9,   -1,  -10,   -6,    1,   -7,   -4,   -4,  -17,  -12,   -3,  -12,   -9 ,
     -7,    8,    1,   -9,   -1,   -5,   -7,    0,   -4,   -4,   11,    3,  -10,    2,   -3,   -6,    3,   -2,   -5,    7,    1,  -10,    0,   -4,   -8,    2,   -4,   -7,    4,   -3,  -10,   -2,   -7,   -9,   -2,   -6 ,
     -2,    0,    3,   -4,  -11,   -4,   -3,   -6,   -2,    0,    3,    6,   -5,   -8,   -2,   -2,   -3,    1,   -2,    0,    3,   -6,   -8,   -4,   -4,   -5,   -2,   -4,   -2,   -1,   -8,   -9,   -6,   -6,   -7,   -4 ,
     -1,  -11,   -5,    7,  -13,   -5,    3,  -11,   -6,    0,  -10,   -5,    9,   -8,   -1,    4,   -9,   -4,   -2,   -8,   -7,    6,   -9,   -3,    1,   -8,   -5,   -5,  -10,   -9,    5,  -12,   -5,   -1,  -12,   -7 ,
    -16,   -1,   -8,  -11,    3,   -1,  -13,    0,   -6,  -16,    2,   -8,   -8,    9,    3,  -12,    4,   -3,  -14,   -1,   -9,   -9,    5,   -1,  -13,    0,   -6,  -17,   -4,  -12,  -11,    5,   -3,  -14,    0,   -8 ,
     -9,   -3,   -2,   -3,   -3,    3,   -5,   -3,    0,   -8,   -3,   -2,   -1,    3,    7,   -4,    1,    2,   -9,   -5,   -4,   -2,    0,    3,   -7,   -3,   -1,  -11,   -7,   -7,   -4,    0,    2,   -7,   -3,   -3 ,
      1,   -7,   -3,    2,  -15,   -7,    4,  -10,   -4,    3,   -6,   -2,    4,  -12,   -4,    6,   -6,   -1,   -1,   -7,   -4,    2,  -11,   -5,    3,   -8,   -4,   -3,   -9,   -7,    0,  -12,   -7,    1,   -8,   -6 ,
    -11,    2,   -3,  -10,    0,   -2,   -7,    3,   -1,  -10,    3,   -3,   -9,    4,    1,   -6,    8,    2,  -11,    1,   -5,  -10,    1,   -2,   -8,    4,   -2,  -13,   -2,   -8,  -10,    1,   -4,   -9,    4,   -3 ,
     -6,   -3,    0,   -5,   -7,   -1,   -3,   -2,    1,   -5,   -2,    1,   -4,   -3,    2,   -1,    2,    4,   -7,   -4,   -2,   -5,   -5,   -1,   -4,   -2,    1,  -10,   -5,   -5,   -6,   -5,   -2,   -5,   -2,   -1 ,
     -2,   -6,   -5,   -5,  -17,  -12,   -3,  -13,  -10,    3,   -5,   -2,   -2,  -14,   -9,   -1,  -11,   -7,    9,   -2,    3,    2,  -10,   -5,    4,   -6,   -3,    4,   -3,    0,   -1,  -13,   -8,    0,   -8,   -5 ,
     -8,    4,   -3,  -10,   -3,   -5,   -8,   -2,   -6,   -6,    7,    0,   -8,   -1,   -5,   -7,    1,   -4,   -2,   13,    5,   -7,    5,    0,   -5,    7,    1,   -5,    8,    1,   -6,   -1,   -4,   -6,    1,   -2 ,
     -5,   -2,   -1,   -8,  -12,   -7,   -5,   -8,   -5,   -2,    1,    3,   -7,   -9,   -4,   -4,   -5,   -2,    3,    5,    8,   -3,   -4,    1,   -1,   -1,    3,   -2,    1,    3,   -5,   -9,   -4,   -4,   -4,   -2 ,
     -5,  -12,   -6,    2,  -12,   -9,   -1,  -12,  -10,   -2,  -10,   -6,    6,   -9,   -2,    2,  -10,   -5,    2,   -7,   -3,   12,   -2,    3,    6,   -4,   -1,   -2,   -8,   -5,    7,   -8,   -1,    2,   -8,   -4 ,
    -21,   -2,   -7,  -12,    0,   -4,  -13,   -2,   -7,  -14,    0,   -8,   -9,    5,    0,  -11,    1,   -5,  -10,    5,   -4,   -2,   12,    5,   -7,    7,    0,  -13,   -1,   -7,   -8,    6,    0,   -9,    1,   -4 ,
    -11,   -4,   -4,   -6,   -6,   -2,   -7,   -5,   -4,   -9,   -4,   -4,   -3,   -1,    3,   -5,   -2,   -1,   -5,    0,    1,    3,    5,    9,   -1,    3,    4,   -9,   -4,   -4,   -1,    0,    4,   -5,   -3,   -1 ,
     -4,  -11,   -6,   -3,  -16,  -10,   -1,  -11,   -7,   -1,   -8,   -4,    1,  -13,   -7,    3,   -8,   -4,    4,   -5,   -1,    6,   -7,   -1,    8,   -4,    1,   -1,   -6,   -4,    2,  -11,   -5,    4,   -6,   -2 ,
    -12,    0,   -5,  -10,   -2,   -4,   -9,    0,   -4,  -10,    2,   -5,   -8,    0,   -3,   -8,    4,   -2,   -6,    7,   -1,   -4,    7,    3,   -4,   11,    3,   -9,    1,   -4,   -7,    2,   -1,   -6,    5,    0 ,
     -8,   -4,   -3,   -7,   -9,   -4,   -5,   -5,   -2,   -6,   -4,   -2,   -5,   -6,   -1,   -4,   -2,    1,   -3,    1,    3,   -1,    0,    4,    1,    3,    6,   -7,   -3,   -2,   -4,   -5,   -1,   -3,   -1,    1 ,
     -5,  -10,   -8,   -8,  -23,  -15,   -7,  -18,  -13,    1,   -7,   -4,   -5,  -17,  -11,   -3,  -13,  -10,    4,   -5,   -2,   -2,  -13,   -9,   -1,   -9,   -7,    9,   -1,    3,    2,  -11,   -5,    4,   -8,   -3 ,
    -10,   -7,   -6,  -17,   -6,   -8,  -12,   -7,   -9,   -7,    4,   -2,  -10,   -4,   -7,   -9,   -2,   -5,   -3,    8,    1,   -8,   -1,   -4,   -6,    1,   -3,   -1,   12,    5,   -4,    3,   -1,   -3,    5,    1 ,
     -8,   -5,   -4,  -11,  -15,  -10,   -9,  -13,   -8,   -4,   -3,   -1,   -9,  -12,   -7,   -7,   -8,   -5,    0,    1,    3,   -5,   -7,   -4,   -4,   -4,   -2,    3,    5,    7,   -2,   -5,    0,    0,   -1,    2 ,
     -8,  -16,  -11,   -3,  -19,  -13,   -5,  -17,  -12,   -4,  -10,   -8,    5,  -11,   -4,    0,  -10,   -6,   -1,   -6,   -5,    7,   -8,   -1,    2,   -7,   -4,    2,   -4,   -2,   12,   -5,    3,    6,   -5,    0 ,
    -19,   -6,  -12,  -17,   -2,   -7,  -16,   -4,  -11,  -17,   -2,   -9,  -12,    5,    0,  -12,    1,   -5,  -13,   -1,   -9,   -8,    6,    0,  -11,    2,   -5,  -11,    3,   -5,   -5,   11,    4,   -7,    7,    0 ,
    -15,   -9,   -7,  -10,   -9,   -5,  -11,   -9,   -7,  -12,   -7,   -6,   -5,   -3,    2,   -7,   -4,   -2,   -8,   -4,   -4,   -1,    0,    4,   -5,   -1,   -1,   -5,   -1,    0,    3,    4,    9,   -1,    2,    4 ,
     -7,  -12,   -9,   -7,  -18,  -14,   -4,  -14,  -10,   -3,   -9,   -6,   -1,  -14,   -7,    1,   -9,   -5,    0,   -6,   -4,    2,   -9,   -5,    4,   -6,   -3,    4,   -3,    0,    6,   -7,   -1,    8,   -3,    1 ,
    -14,   -4,   -9,  -18,   -6,   -9,  -12,   -3,   -7,  -12,   -2,   -7,  -12,    0,   -3,   -8,    4,   -2,   -8,    1,   -4,   -8,    1,   -3,   -6,    5,   -1,   -8,    5,   -1,   -5,    7,    2,   -3,   10,    3 ,
    -12,   -6,   -7,  -12,  -13,   -8,   -9,   -8,   -6,   -9,   -6,   -4,   -7,   -8,   -3,   -6,   -3,   -1,   -5,   -2,   -2,   -4,   -4,   -1,   -2,    0,    1,   -3,    1,    2,    0,    0,    4,    1,    3,    6 ,
};

void Paralign::SetMu_scop40_tm0_6_0_8_fa2()
	{
	m_Open = 16;
	m_Ext = 5;
	int MinScore = 0;
	int MaxScore = 0;
	for (uint i = 0; i < 36*36; ++i)
		{
		int Score = Mu_S_k_i8[i];
		if (i == 0 || Score < MinScore) MinScore = Score;
		if (i == 0 || Score > MaxScore) MaxScore = Score;
		}
	m_matrix.size = 36;
	m_matrix.length = 36;
	m_matrix.type = PARASAIL_MATRIX_TYPE_SQUARE;
	m_matrix.matrix = scop40_tm0_6_0_8_fa2;
	m_matrix.min = MinScore;
	m_matrix.max = MaxScore;
	int *Mapper = myalloc(int, 256);
	memset(Mapper, 0, 256*sizeof(int));
	for (int i = 0; i < 36; ++i)
		Mapper[i] = i;
	m_matrix.mapper = Mapper;
	m_Bits = 16;
	SetSWFastSubstMx_FromParasailMx();
	}

void Paralign::SetBlosum62()
	{
	vector<vector<int> > ScoreMx(20);
	for (uint i = 0; i < 20; ++i)
		{
		ScoreMx[i].resize(20);
		for (uint j = 0; j < 20; ++j)
			ScoreMx[i][j] = Blosum62_int[i][j];
		}
	SetMatrix(ScoreMx, Blosum62_Open, Blosum62_Ext, Blosum62_SaturatedScore);
	}

void Paralign::LogMatrix()
	{
	int AS = m_matrix.size;
	Log("Paralign::LogMatrix() open %d, ext %d, min %d, max %d, size %d, length %d\n",
		m_Open,
		m_Ext,
		m_matrix.min,
		m_matrix.max,
		m_matrix.size,
		m_matrix.length);

	Log("int Paralign_IntMx[36*36] = {\n");
	for (int i = 0; i < AS; ++i)
		{
		for (int j = 0; j < AS; ++j)
			Log("%3d,", m_matrix.matrix[i*AS + j]);
		Log("  // %d\n", i);
		}
	Log("};\n");
	}

void Paralign::LogSWFastMatrix()
	{
	int AS = m_matrix.size;

	Log("float Paralign_SWFastMx[36*36] = {\n");
	for (int i = 0; i < AS; ++i)
		{
		for (int j = 0; j < AS; ++j)
			Log("%7.3g,", m_SWFastSubstMx[i][j]);
		Log("  // %d\n", i);
		}
	Log("};\n");
	}

void Paralign::SetSubstMxByName(const string &Name)
	{
	m_SubstMxName = Name;
	if (Name == "Mu_S_k_i8")
		Set_Mu_S_k_i8();
	else if (Name == "Mu_hjmux")
		Set_Mu_hjmux();
	else if (Name == "Mu_scop40_tm0_6_0_8_fa2")
		SetMu_scop40_tm0_6_0_8_fa2();
	else if (Name == "musubstmx")
		SetMu_musubstmx();
	else if (Name == "parasail_mu_8")
		SetMu_parasail_mu_8();
	else
		Die("SetSubstMx(%s)", Name.c_str());
	}

void Paralign::SetMatrix(
	const vector<vector<int> > &ScoreMx,
	int Open, int Ext, int SaturatedScore,
	bool SetSWFastMatrix)
	{
	asserta(!ScoreMx.empty());
	memset(&m_matrix, 0, sizeof(m_matrix));
	m_Open = Open;
	m_Ext = Ext;
	m_SaturatedScore = SaturatedScore;
	int MinScore = 0;
	int MaxScore = 0;
	int AS = SIZE(ScoreMx[0]);
	int *ScoreVec = myalloc(int, AS*AS);
	int *Mapper = myalloc(int, 256);
	char *Alphabet = myalloc(char, AS);
	memset(Mapper, 0, 256*sizeof(int));
	for (int i = 0; i < AS; ++i)
		{
		asserta(SIZE(ScoreMx[i]) == AS);
		Mapper[i] = i;

	// Up to 64 ASCII
		if (i < 26)
			Alphabet[i] = 'A' + i;
		else if (i < 26*2)
			Alphabet[i] = 'a' + i - 26;
		else if (i < 26*2 + 10)
			Alphabet[i] = '0' + i - 26*2;
		else if (i == 62)
			Alphabet[i] = '@';
		else if (i == 63)
			Alphabet[i] = '$';
		else
			Alphabet[i] = 0;

		for (int j = 0; j < AS; ++j)
			{
			int Score = ScoreMx[i][j];
			if (i == 0 && j == 0)
				{
				MinScore = Score;
				MaxScore = Score;
				}
			else
				{
				MinScore = min(Score, MinScore);
				MaxScore = max(Score, MaxScore);
				}
			ScoreVec[AS*i + j] = Score;
			}
		}

	m_matrix.name = "Paralign::m_matrix";
	m_matrix.matrix = ScoreVec;
	m_matrix.mapper = Mapper;
	m_matrix.size = AS;
	m_matrix.min = MinScore;
	m_matrix.max = MaxScore;
	m_matrix.user_matrix = 0;
	m_matrix.type = PARASAIL_MATRIX_TYPE_SQUARE;
	m_matrix.length = AS;
	m_matrix.alphabet = Alphabet;
	m_matrix.query = 0;
	if (SetSWFastMatrix)
		SetSWFastSubstMx_FromParasailMx();
	}

void Paralign::SetQueryNoProfile(const string &LabelQ, const byte *Q, uint LQ)
	{
	m_LabelQ = LabelQ;
	m_Q = Q;
	m_LQ = LQ;
	}

void Paralign::SetQueryProfile(const string &LabelQ, const byte *Q, uint LQ)
	{
	m_LabelQ = LabelQ;
	m_Q = Q;
	m_LQ = LQ;
	if (m_ProfQ != 0)
		parasail_profile_free(m_ProfQ);
	switch (m_Bits)
		{
	case 8:
		m_ProfQ = parasail_profile_create_avx_256_8((const char *) Q, LQ, &m_matrix);
		break;

	case 16:
		m_ProfQ = parasail_profile_create_avx_256_16((const char *) Q, LQ, &m_matrix);
		break;

	default:
		asserta(false);
		}
#if 0
	{
	Log("SetQuery\n");
	Log_parasail_mu_matrix(m_matrix);
	Log("QL %u\n", m_LQ);
	Log("Q: ");
	for (uint i = 0; i < m_LQ; ++i)
		Log(" %u", m_Q[i]);
	Log("\n");
	Log("ProfQ:");
	const byte *ptrProf = (const byte *) m_ProfQ;
	for (uint i = 0; i < sizeof(*m_ProfQ); ++i)
		Log(" %02x", ptrProf[i]);
	Log("\n");
	}
#endif
	}

void Paralign::Align_ScoreOnly(const string &LabelT, const byte *T, uint LT)
	{
	ClearResult();
	m_T = T;
	m_LT = LT;
	m_LabelT = LabelT;
	if (m_result != 0)
		parasail_result_free(m_result);

	switch (m_Bits)
		{
	case 8:
		m_result = parasail_sw_striped_profile_avx2_256_8(
			m_ProfQ, (const char *) T, LT, m_Open, m_Ext);
		++m_Count8;
		break;

	case 16:
		m_result = parasail_sw_striped_profile_avx2_256_16(
			m_ProfQ, (const char *) T, LT, m_Open, m_Ext);
		++m_Count16;
		break;

	default:
		asserta(false);
		}

#if 0
	{
	Log("QL %u, TL %u\n", m_LQ, m_LT);
	Log("Q: ");
	for (uint i = 0; i < m_LQ; ++i)
		Log(" %u", m_Q[i]);
	Log("\n");
	Log("T: ");
	for (uint i = 0; i < m_LT; ++i)
		Log(" %u", m_T[i]);
	Log("\n");
	Log("ProfQ:");
	const parasail_profile_t *prof =
		(parasail_profile_t *) m_ProfQ;
	Log_parasail_profile(*prof);
	Log("\n");
	Log("score %d\n", m_result->score);
	}
#endif
	if (m_result->flag & PARASAIL_FLAG_SATURATED)
		{
		m_Score = m_SaturatedScore;
		++m_SaturatedCount;
		}
	else
		m_Score = m_result->score;
	}

bool Paralign::Align_Path(const string &LabelT, const byte *T, uint LT)
	{
	ClearResult();
	m_LabelT = LabelT;
	m_T = T;
	m_LT = LT;
	if (m_result != 0)
		parasail_result_free(m_result);
	m_result = parasail_sw_trace_striped_profile_avx2_256_8(
		m_ProfQ, (const char *) T, LT, m_Open, m_Ext);
	if (m_result->flag & PARASAIL_FLAG_SATURATED)
		{
		m_Score = m_SaturatedScore;
		return false;
		}

	m_Score = m_result->score;
	parasail_cigar_t* cig = parasail_result_get_cigar_extra(
		m_result,
		(const char *) m_Q, m_LQ,
		(const char *) m_T, m_LT,
		&m_matrix, 1, 0);

	char *cig_str = parasail_cigar_decode(cig);
	m_LoQ = (uint) cig->beg_query;
	m_LoT = (uint) cig->beg_ref;
	ExpandParaCigar_reverseDI(cig_str, m_SemiGlobalPath);
	free(cig_str);
	parasail_cigar_free(cig);

#if DEBUG
	{
	uint M, D, I;
	GetPathCounts(m_SemiGlobalPath, M, D, I);
	asserta(m_SemiGlobalPath.back() == 'M');
	asserta(M + D <= m_LQ);
	asserta(M + I <= m_LT);
	}
#endif
	return true;
	}

void Paralign::SetCompoundMx(
	const vector<FEATURE> &Fs, const vector<float> &Weights,
	int ScaleFactor, int Open, int Ext, int SaturatedScore)
	{
	const uint NF = SIZE(Fs);
	asserta(SIZE(Weights) == NF);

	vector<const float * const *> ScoreMxs;
	for (uint i = 0; i < NF; ++i)
		ScoreMxs.push_back(DSSParams::GetScoreMx(Fs[i]));

	Nu TheNu;
	TheNu.SetComponents(Fs, Weights);
	uint AS = TheNu.GetAlphaSize();

	m_SWFastSubstMx.clear();
	m_SWFastSubstMx.resize(AS);
	vector<vector<int> > IntScoreMx(AS);
	for (uint i = 0; i < AS; ++i)
		{
		m_SWFastSubstMx[i].resize(AS);
		IntScoreMx[i].resize(AS);

		vector<byte> Lettersi;
		TheNu.NuLetterToComponentLetters(i, Lettersi);
		asserta(SIZE(Lettersi) == NF);


		for (uint j = 0; j < AS; ++j)
			{
			vector<byte> Lettersj;
			TheNu.NuLetterToComponentLetters(j, Lettersj);
			asserta(SIZE(Lettersj) == NF);

			float SumScore = 0;
			for (uint k = 0; k < NF; ++k)
				{
				float wk = NF*Weights[k];
				uint Letteri = Lettersi[k];
				uint Letterj = Lettersj[k];
				SumScore += wk*ScoreMxs[k][Letteri][Letterj];
				}
			float Score = ScaleFactor*SumScore;
			m_SWFastSubstMx[i][j] = Score;
			IntScoreMx[i][j] = int(round(Score));
			}
		}
	bool SetSWFastMx = false;
	if (opt(roundmx))
		SetSWFastMx = true;
	SetMatrix(IntScoreMx, Open, Ext, SaturatedScore, SetSWFastMx);
	}

#if 0
void cmd_paralign_test()
	{
	const string &FastaFN = g_Arg1;
	SeqDB Seqs;
	Seqs.FromFasta(FastaFN);
	FILE *fOut = CreateStdioFile(opt(output));
	FILE *fAln = CreateStdioFile(opt(aln));

	const uint SeqCount = Seqs.GetSeqCount();
	vector<vector<byte> > ByteSeqs(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		const string &Seq = Seqs.GetSeq(i);
		const uint L = Seqs.GetSeqLength(i);
		vector<byte> &ByteSeq = ByteSeqs[i];
		ByteSeq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			char c = Seq[Pos];
			byte Letter = g_CharToLetterMu[c];
			ByteSeq.push_back(Letter);
			}
		}

	Paralign PA;
	Paralign::Set_Mu_S_k_i8();
	Paralign::SetSWFastSubstMx_FromParasailMx();
	uint Pairs = 0;
	uint ScoreDiffs = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		ProgressStep(i, SeqCount, "Aligning");
		const string &Label_i = Seqs.GetLabel(i);
		const string &Seq_i = Seqs.GetSeq(i);
		const byte *ByteSeq_i = ByteSeqs[i].data();
		uint L_i = Seqs.GetSeqLength(i);
		PA.SetQueryProfile(Label_i, ByteSeq_i, L_i);
		for (uint j = i+1; j < SeqCount; ++j)
			{
			++Pairs;
			const string &Label_j = Seqs.GetLabel(j);
			const string &Seq_j = Seqs.GetSeq(j);
			const byte *ByteSeq_j = ByteSeqs[j].data();
			uint L_j = Seqs.GetSeqLength(j);
			PA.Align_ScoreOnly(Label_j, ByteSeq_j, L_j);
			int Score = PA.m_Score;
			bool Ok = PA.Align_Path(Label_j, ByteSeq_j, L_j);
			int LoQ = PA.m_LoQ;
			int LoT = PA.m_LoT;
			int Score2 = PA.ScoreAln();
			string CIGAR;
			LocalPathToCIGAR(PA.m_SemiGlobalPath.c_str(), 0, 0, CIGAR, false);
			if (fAln != 0)
				PA.WriteAln(fAln);
			PA.Align_SWFast(Label_j, ByteSeq_j, L_j);
			float Score3 = PA.m_SWFastScore;
			int Score4 = PA.m_SWFastScoreInt;
			char isdiff = 'n';
			if (Score + Score2 + Score3 + Score4 != 4*Score)
				{
				isdiff = 'Y';
				++ScoreDiffs;
				}
			fprintf(fOut, "%s", Label_i.c_str());
			fprintf(fOut, "\t%s", Label_j.c_str());
			fprintf(fOut, "\t%c", isdiff);
			fprintf(fOut, "\t%d", Score);
			fprintf(fOut, "\t%d", Score2);
			fprintf(fOut, "\t%.1f", Score3);
			fprintf(fOut, "\t%d", Score4);
			fprintf(fOut, "\t%d", LoQ);
			fprintf(fOut, "\t%d", LoT);
			fprintf(fOut, "\t%s", CIGAR.c_str());
			fprintf(fOut, "\n");
			}
		}
	CloseStdioFile(fOut);
	CloseStdioFile(fAln);
	ProgressLog("%u / %u score diffs\n", ScoreDiffs, Pairs);
	}
#endif