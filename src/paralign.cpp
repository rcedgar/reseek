#include "myutils.h"
#include "paralign.h"
#include "seqdb.h"
#include "alpha.h"
#include "cigar.h"
#include "flat_params.h"

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

parasail_matrix_t Paralign::m_matrix;
int Paralign::m_Open = INT_MIN;	// penalty > 0
int Paralign::m_Ext = INT_MIN;	// penalty > 0
int Paralign::m_SaturatedScore = INT_MAX;
int Paralign::m_Bits = 16;
string Paralign::m_SubstMxName = "_NOT_SET_";
vector<vector<float> > Paralign::m_SWFastSubstMx;
atomic<uint> Paralign::m_Count8;
atomic<uint> Paralign::m_Count16;
atomic<uint> Paralign::m_Count8_rev;
atomic<uint> Paralign::m_Count16_rev;
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

// Note paralign generates paths with DI and ID transitions,
// these are scored with gapopen
int Paralign::score_nu_path(
	const string &label_i, const uint8_t *nu_codeseq_i, uint lo_i, uint L_i,
	const string &label_j, const uint8_t *nu_codeseq_j, uint lo_j, uint L_j,
	const string &path)
	{
	uint pos_i = lo_i;
	uint pos_j = lo_j;
	uint score = 0;
	const uint ncol = uint(path.size());
	char last_c = 'M';
	for (uint col = 0; col < ncol; ++col)
		{
		char c = path[col];
		if (c == 'M')
			{
			asserta(pos_i < L_i);
			asserta(pos_j < L_j);
			uint8_t code_i = nu_codeseq_i[pos_i];
			uint8_t code_j = nu_codeseq_j[pos_j];
			score += m_matrix.matrix[code_i*256 + code_j];
			++pos_i;
			++pos_j;
			last_c = 'M';
			}
		else if (c == 'D')
			++pos_i;
		else if (c == 'I')
			++pos_j;
		if (c == 'D' || c == 'I')
			{
			if (last_c == c)
				score -= m_Ext;
			else
				score -= m_Open;
			}
		last_c = c;
		}
	return score;
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

#include "final_nu_matrix.h"

/***
$src/reseek_tune2/bash/final_nu_matrix.bash
                                 vvvvvvvvvvvvvv--- scale pre-built into matrix
intopen=2.90E+01;intext=3.00E+00;scale=8.81E+00;aa4=5.15E-01;pm2=2.84E-01;sec32=2.00E-01;
***/
void Paralign::set_final_nu()
	{
	m_Open = 29;
	m_Ext = 3;
	int MinScore = 0;
	int MaxScore = 0;
	for (uint i = 0; i < 256*256; ++i)
		{
		int Score = s_final_nu_matrix[i];
		if (i == 0 || Score < MinScore) MinScore = Score;
		if (i == 0 || Score > MaxScore) MaxScore = Score;
		}
	m_matrix.size = 256;
	m_matrix.length = 256;
	m_matrix.type = PARASAIL_MATRIX_TYPE_SQUARE;
	m_matrix.matrix = s_final_nu_matrix;
	m_matrix.min = MinScore;
	m_matrix.max = MaxScore;
	int *Mapper = myalloc(int, 256);
	memset(Mapper, 0, 256*sizeof(int));
	for (int i = 0; i < 256; ++i)
		Mapper[i] = i;
	m_matrix.mapper = Mapper;
	m_Bits = 16;
	}

static bool do_init()
	{
	Paralign::set_final_nu();
	return true;
	}
static bool s_init_done = do_init();

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

	Log("int Paralign_IntMx[%u*%u] = {\n", AS, AS);
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

void Paralign::SetQueryProfile_rev(const byte *Q, uint LQ)
	{
	if (!m_DoReverse)
		return;
	if (m_Q_rev != 0)
		myfree(m_Q_rev);
	if (m_ProfQ_rev != 0)
		parasail_profile_free(m_ProfQ_rev);

	m_Q_rev = myalloc(byte, LQ);
	for (uint i = 0; i < LQ; ++i)
		m_Q_rev[i] = Q[LQ-i-1];

	switch (m_Bits)
		{
	case 8:
		m_ProfQ_rev = parasail_profile_create_avx_256_8((const char *) m_Q_rev, LQ, &m_matrix);
		break;

	case 16:
		m_ProfQ_rev = parasail_profile_create_avx_256_16((const char *) m_Q_rev, LQ, &m_matrix);
		break;

	default:
		asserta(false);
		}
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
	if (m_DoReverse)
		SetQueryProfile_rev(Q, LQ);

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

void Paralign::Align_ScoreOnly_rev(const string &LabelT, const byte *T, uint LT)
	{
	ClearResult();
	if (m_result != 0)
		parasail_result_free(m_result);

	switch (m_Bits)
		{
	case 8:
		m_result = parasail_sw_striped_profile_avx2_256_8(
			m_ProfQ_rev, (const char *) T, LT, m_Open, m_Ext);
		++m_Count8_rev;
		break;

	case 16:
		m_result = parasail_sw_striped_profile_avx2_256_16(
			m_ProfQ_rev, (const char *) T, LT, m_Open, m_Ext);
		++m_Count16_rev;
		break;

	default:
		asserta(false);
		}

	if (m_result->flag & PARASAIL_FLAG_SATURATED)
		{
		m_Score_rev = m_SaturatedScore;
		++m_SaturatedCount;
		}
	else
		m_Score_rev = m_result->score;
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
	asserta(!m_DoReverse);
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

void Paralign::set_flat_compound(
	flat_params &params,
	const unordered_map<string, float> &name2weight,
	float ScaleFactor,
	int Open,
	int Ext,
	int SaturatedScore)
	{
	asserta(ScaleFactor > 0.1);
	params.apply_weights(name2weight);
	const uint compound_alpha_size = params.get_compound_alpha_size();

	m_SWFastSubstMx.clear();
	m_SWFastSubstMx.resize(compound_alpha_size);
	asserta(compound_alpha_size <= 256);
	vector<vector<int> > IntScoreMx(compound_alpha_size);
	for (uint i = 0; i < compound_alpha_size; ++i)
		{
		m_SWFastSubstMx[i].resize(compound_alpha_size);
		IntScoreMx[i].resize(compound_alpha_size);
		const uint8_t code_i = uint8_t(i);
		for (uint j = 0; j < compound_alpha_size; ++j)
			{
			const uint8_t code_j = uint8_t(j);
			float Score = ScaleFactor*
				params.get_compound_subst_score_slow(code_i, code_j);
			m_SWFastSubstMx[i][j] = Score;
			int IntScore = int(round(Score));
			IntScoreMx[i][j] = IntScore;
			}
		}
	bool SetSWFastMx = false;
	if (opt(roundmx))
		SetSWFastMx = true;
	SetMatrix(IntScoreMx, Open, Ext, SaturatedScore, SetSWFastMx);
	if (opt(logmx))
		{
		LogMatrix();
		Die("-logmx");
		}
	}
