#include "myutils.h"
#include "paralign.h"
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

parasail_matrix_t Paralign::m_matrix;
int Paralign::m_Open = INT_MAX;	// penalty > 0
int Paralign::m_Ext = INT_MAX;	// penalty > 0
int Paralign::m_SaturatedScore = INT_MAX;
uint Paralign::m_MaxLength = 9999;
int Paralign::m_Bits = 16;
vector<vector<float> > Paralign::m_SWFastSubstMx;
atomic<uint> Paralign::m_Count8;
atomic<uint> Paralign::m_Count16;
atomic<uint> Paralign::m_TooLongCount;
atomic<uint> Paralign::m_SaturatedCount;
     
// Ye olde BLOSUM62 as used by NCBI BLAST (1/2-bit units)
// alphabetical order, no wildcards or stop codon
static const int Blosum62_Open = 11;
static const int Blosum62_Ext = 1;
static const int Blosum62_SaturatedScore = 999;

void Paralign::SetSWFastSubstMx()
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

void Paralign::Align_SWFast()
	{
	const uint AS = GetAlphaSize();
	asserta(SIZE(m_SWFastSubstMx) == AS);
	float Open = -float(m_Open);
	float Ext = -float(m_Ext);
	uint LoQ, LenQ, LoT, LenT;
	m_SWFastScore = SWFast_SubstMx(m_Mem, m_Q, m_LQ, m_T, m_LT,
		m_SWFastSubstMx, Open, Ext, LoQ, LenQ, LoT, LenT, m_SWFastPath);
	m_SWFastScoreInt = int(round(m_SWFastScore));
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

void Paralign::SetMu()
	{
	//m_Open = 2;
	//m_Ext = 1;
	//m_SaturatedScore = 777;
	//memcpy(&m_matrix, &parasail_mu_matrix, sizeof(parasail_matrix_t));

// $src/2025-10_reseek_tune/2025-11-22_paralign_scop40_mu_16bit_gap_parameter_sweep
// open20.ext8=0.724
// open24.ext10=0.724
// open30.ext5=0.724
// open35.ext5=0.724
// open18.ext10=0.726
// open22.ext10=0.726
// open22.ext8=0.726
// open24.ext8=0.726
// open20.ext10=0.727
	m_Open = 24;
	m_Ext = 8;
	extern int Mu_S_k_i8[36*36];
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

void Paralign::SetMatrix(
	const vector<vector<int> > &ScoreMx,
	int Open, int Ext, int SaturatedScore)
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
	}

void Paralign::SetQuery(const string &LabelQ, const byte *Q, uint LQ)
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
	if (m_LQ > m_MaxLength || m_LT > m_MaxLength)
		{
		m_Score = -999;//@@TODO
		++m_TooLongCount;
		return;
		}
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
	if (m_LQ > m_MaxLength || m_LT > m_MaxLength)
		{
		m_Score = -999;//@@TODO
		return false;
		}
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
	Paralign::SetMu();
	Paralign::SetSWFastSubstMx();
	uint Pairs = 0;
	uint ScoreDiffs = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		ProgressStep(i, SeqCount, "Aligning");
		const string &Label_i = Seqs.GetLabel(i);
		const string &Seq_i = Seqs.GetSeq(i);
		const byte *ByteSeq_i = ByteSeqs[i].data();
		uint L_i = Seqs.GetSeqLength(i);
		PA.SetQuery(Label_i, ByteSeq_i, L_i);
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
			PA.Align_SWFast();
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
