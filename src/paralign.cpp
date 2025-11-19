#include "myutils.h"
#include "paralign.h"
#include "seqdb.h"
#include "alpha.h"

void ExpandParaCigar_reverseDI(const string &s, string &Path);

extern parasail_matrix_t parasail_mu_matrix;
extern int Blosum62_int[20][20];

parasail_matrix_t Paralign::m_matrix;
int Paralign::m_Open = INT_MAX;
int Paralign::m_Ext = INT_MAX;
int Paralign::m_SaturatedScore = INT_MAX;
uint Paralign::m_MaxLength = 415;

// Ye olde BLOSUM62 as used by NCBI BLAST (1/2-bit units)
// alphabetical order, no wildcards or stop codon
static const int Blosum62_Open = 11;
static const int Blosum62_Ext = 1;
static const int Blosum62_SaturatedScore = 999;

void Paralign::SetMu()
	{
	m_Open = 2;
	m_Ext = 1;
	m_SaturatedScore = 777;
	memcpy(&m_matrix, &parasail_mu_matrix, sizeof(parasail_matrix_t));
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
	m_matrix.max = MinScore;
	m_matrix.min = MaxScore;
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
	m_ProfQ = parasail_profile_create_avx_256_8((const char *) Q, LQ, &m_matrix);
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
		return;
		}
	if (m_result != 0)
		parasail_result_free(m_result);
	m_result = parasail_sw_striped_profile_avx2_256_8(
		  m_ProfQ, (const char *) T, LT, m_Open, m_Ext);
	if (m_result->flag & PARASAIL_FLAG_SATURATED)
		m_Score = m_SaturatedScore;
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
	ExpandParaCigar_reverseDI(cig_str, m_Path);
	free(cig_str);
	parasail_cigar_free(cig);

	return true;
	}

#if 0
void cmd_paralign_test()
	{
	const string Label2("d1ujva_/b.36.1.1");
	const char *Seq2 = "RRQQQQQQQQNNNffifijjhbdjgfZiigjXaXXRRIQaaaMPOOOOQRRIHEFFOMOLORRaaaXVYPMSYSbdgfgggffNQNNNOOOONOOO";

	const string Label3("d1vkya_/e.53.1.1");
	const char *Seq3 = "OROONNNRRPROOggijgjghRLRQQQQQRPRRO";

	const uint L2 = ustrlen(Seq2);
	const uint L3 = ustrlen(Seq3);
	byte *ByteSeq2 = myalloc(byte, L2);
	byte *ByteSeq3 = myalloc(byte, L3);
	for (uint i = 0; i < L2; ++i)
		ByteSeq2[i] = g_CharToLetterMu[Seq2[i]];
	for (uint i = 0; i < L2; ++i)
		ByteSeq3[i] = g_CharToLetterMu[Seq3[i]];

	Paralign::SetMu();
	Paralign PA;
	PA.SetQuery(Label3, ByteSeq3, L3);
	PA.Align_ScoreOnly(Label2, ByteSeq2, L2);
	ProgressLog("Score %d\n", PA.m_Score);
	}
#else
void cmd_paralign_test()
	{
	const string &FastaFN = g_Arg1;
	SeqDB Seqs;
	Seqs.FromFasta(FastaFN);
	FILE *fOut = CreateStdioFile(opt(output));

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
			const string &Label_j = Seqs.GetLabel(j);
			const string &Seq_j = Seqs.GetSeq(j);
			const byte *ByteSeq_j = ByteSeqs[j].data();
			uint L_j = Seqs.GetSeqLength(j);
			//bool Ok = PA.Align_Path(Label_j, ByteSeq_j, L_j);
			PA.Align_ScoreOnly(Label_j, ByteSeq_j, L_j);
			fprintf(fOut, "%s", Label_i.c_str());
			fprintf(fOut, "\t%s", Label_j.c_str());
			fprintf(fOut, "\t%d", PA.m_Score);
			fprintf(fOut, "\n");
			}
		}
	CloseStdioFile(fOut);
	}
#endif
