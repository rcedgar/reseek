#include "myutils.h"
#include "parasail.h"
#include "cigar.h"

extern int8_t IntScoreMx_Combo[36][36];
extern parasail_matrix_t parasail_combo_matrix;

char GetComboChar(byte Letter)
	{
	if (Letter < 26)
		return 'A' + Letter;
	else if (Letter < 36)
		return 'a' + (Letter - 26);
	asserta(false);
	return '!';
	}

void LogAlnComboLetters(
  const vector<byte> &A, uint LoA,
  const vector<byte> &B, uint LoB,
  const string &Path)
	{
	uint LA = SIZE(A);
	uint LB = SIZE(B);
	const uint ColCount = SIZE(Path);
	uint PosA = LoA;
	uint PosB = LoB;
	string RowA;
	string RowB;
	string RowM;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			{
			asserta(PosA < LA);
			asserta(PosB < LB);
			byte a = A[PosA++];
			byte b = B[PosB++];
			char m = (a == b ? '|' : ' ');
			if (a == b)
				m = '|';
			else if (IntScoreMx_Combo[a][b] > 0)
				m = '+';
			else
				m = ' ';
			RowM += m;

			RowA += GetComboChar(a);
			RowB += GetComboChar(b);
			break;
			}

		case 'I':
			{
			asserta(PosA < LA);
			RowA += GetComboChar(A[PosA++]);
			RowB += '-';
			RowM += ' ';
			break;
			}

		case 'D':
			{
			asserta(PosB < LB);
			RowA += '-';
			RowB += GetComboChar(B[PosB++]);
			RowM += ' ';
			break;
			}

		default:
			asserta(false);
			}
		}
	asserta(PosA <= LA);
	asserta(PosB <= LB);
	Log("%s\n", RowA.c_str());
	Log("%s\n", RowM.c_str());
	Log("%s\n", RowB.c_str());
	}

static void Test1()
	{
	vector<byte> S1;
	vector<byte> S2;
	for (uint i = 0; i < 10; ++i)
		{
		if (i > 0)
			S1.push_back(i);
		if (i != 7)
			S2.push_back(i);
		}
	S1[3] = 0;

	const int s1Len = (int) S1.size();
	const int s2Len = (int) S2.size();

	const int open = 3;
	const int ext = 1;

	const char *s1 = (const char *) S1.data();
	const char *s2 = (const char *) S2.data();
	
	const parasail_profile_t* profile1 = 
	  parasail_profile_create_avx_256_8(s1, s1Len, &parasail_combo_matrix);
	
	parasail_result_t *result1 =
	  parasail_sw_striped_profile_avx2_256_8(profile1, s2, s2Len, open, ext);

	parasail_result_t *result2 =
	  parasail_sw_trace_striped_profile_avx2_256_8(profile1, s2, s2Len, open, ext);

	const parasail_result_extra_trace_t *trace = result2->trace;
	parasail_cigar_t* cig = parasail_result_get_cigar_extra(result2, s1, s1Len, s2, s2Len,
	  &parasail_combo_matrix, 1, 0);
	char *cig_str = parasail_cigar_decode(cig);
	string Path;
	ExpandParaCigar(cig_str, Path);

	uint Lo1 = (uint) cig->beg_query;
	uint Lo2 = (uint) cig->beg_ref;
	Log("score = %d, %d lo %d, %d path %s\n",
	  result1->score, result2->score, Lo1, Lo2, Path.c_str());
	LogAlnComboLetters(S1, Lo1, S2, Lo2, Path);
	}

static byte GetRandLetter()
	{
	return (byte) randu32()%36;
	}

static void GetRandSeq(uint L, vector<byte> &Letters)
	{
	Letters.clear();
	for (uint i = 0; i < L; ++i)
		Letters.push_back(GetRandLetter());
	}

static bool TrueProb(double Prob)
	{
	const uint M = 1024*1024*1024;
	uint k = randu32()%M;
	bool IsTrue = (k <= uint(M*Prob));
	return IsTrue;
	}

static void GetMutatedSeq(vector<byte> &Letters,
  uint LeftSkip, uint RightSkip, 
  double SubFreq, double InsFreq, double DelFreq,
  vector<byte> &MutatedLetters)
	{
	MutatedLetters.clear();
	const uint L = SIZE(Letters);
	for (uint Pos = LeftSkip; Pos + RightSkip < L; ++Pos)
		{
		if (TrueProb(InsFreq))
			MutatedLetters.push_back(GetRandLetter());
		else if (TrueProb(DelFreq))
			continue;
		byte Letter = Letters[Pos];
		if (TrueProb(SubFreq))
			{
			byte NewLetter = GetRandLetter();
			if (NewLetter == Letter)
				Letter = (NewLetter + 1)%36;
			else
				Letter = NewLetter;
			}
		MutatedLetters.push_back(Letter);
		}
	}

int ParasailAlign(const vector<byte> &Q, const vector<byte> &T,
  int GapOpen, int GapExt, uint &LoQ, uint &LoT, string &Path)
	{
	const uint QL = SIZE(Q);
	const uint TL = SIZE(T);
	const char *ptrQ = (const char *) Q.data();
	const char *ptrT = (const char *) T.data();

	parasail_profile_t *QP = 
	  parasail_profile_create_avx_256_8(ptrQ, QL, &parasail_combo_matrix);

	parasail_result_t *result =
	  parasail_sw_trace_striped_profile_avx2_256_8(QP, ptrT, TL, GapOpen, GapExt);
	int score = result->score;
	if (parasail_result_is_saturated(result))
		{
		Warning("parasail_result_is_saturated()");
		return -777;
		}

	parasail_cigar_t* cig = parasail_result_get_cigar_extra(
	  result, ptrQ, QL, ptrT, TL, &parasail_combo_matrix, 1, 0);

	char *cig_str = parasail_cigar_decode(cig);
	ExpandParaCigar(cig_str, Path);

	LoQ = (uint) cig->beg_query;
	LoT = (uint) cig->beg_ref;

	free(cig_str);
	parasail_cigar_free(cig);
	parasail_result_free(result);
	parasail_profile_free(QP);

	return score;
	}

static void Test2(bool Trace)
	{
	vector<byte> S1;
	vector<byte> S2;

	const int s1Len = 100 + randu32()%20;
	GetRandSeq(s1Len, S1);

	uint LeftSkip = randu32()%3;
	uint RightSkip = randu32()%3;
	GetMutatedSeq(S1, LeftSkip, RightSkip, 0.3, 0.1, 0.1, S2);
	const int s2Len = (int) S2.size();

	if (Trace)
		{
		Log("__________________________________________\n");
		Log("S1  ");
		for (int i = 0; i < s1Len; ++i)
			Log("%c", GetComboChar(S1[i]));
		Log("\n");

		Log("S2  ");
		for (int i = 0; i < s2Len; ++i)
			Log("%c", GetComboChar(S2[i]));
		Log("\n");
		}

	const int open = 3;
	const int ext = 1;

	const char *s1 = (const char *) S1.data();
	const char *s2 = (const char *) S2.data();
	
	parasail_profile_t* profile1 = 
	  parasail_profile_create_avx_256_8(s1, s1Len, &parasail_combo_matrix);
	
	parasail_result_t *result1 =
	  parasail_sw_striped_profile_avx2_256_8(profile1, s2, s2Len, open, ext);

	parasail_result_t *result2 =
	  parasail_sw_trace_striped_profile_avx2_256_8(profile1, s2, s2Len, open, ext);

	const parasail_result_extra_trace_t *trace = result2->trace;
	parasail_cigar_t* cig = parasail_result_get_cigar_extra(result2, s1, s1Len, s2, s2Len,
	  &parasail_combo_matrix, 1, 0);
	char *cig_str = parasail_cigar_decode(cig);
	string Path;
	ExpandParaCigar(cig_str, Path);

	uint Lo1 = (uint) cig->beg_query;
	uint Lo2 = (uint) cig->beg_ref;
	if (Trace)
		{
		Log("score = %d, %d lo %d, %d path %s\n",
		  result1->score, result2->score, Lo1, Lo2, Path.c_str());
		LogAlnComboLetters(S1, Lo1, S2, Lo2, Path);
		}

	string Path2;
	uint LoQ, LoT;
	int Score2 = ParasailAlign(S1, S2, open, ext, LoQ, LoT, Path2);
	asserta(Score2 == result1->score);
	asserta(Path2 == Path);

	free(cig_str);
	parasail_cigar_free(cig);
	parasail_profile_free(profile1);
	parasail_result_free(result1);
	parasail_result_free(result2);
	}

void cmd_test_para()
	{
	const uint N = 10000000;
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Working");
		Test2(false);
		}
	}
