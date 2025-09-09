#include "myutils.h"
#include "mermx.h"
#include "alpha.h"
#include "sort.h"
#include "prefiltermuparams.h"
#include "quarts.h"

static mutex g_Lock;
static int64 g_HoodCallCount;

const MerMx &GetMuMerMx(uint k);

void MerMx::KmerToLetters(uint Kmer, uint k, vector<byte> &Letters) const
	{
	Letters.clear();
	for (uint i = 0; i < k; ++i)
		{
		byte Letter = Kmer%m_AS;
		Letters.push_back(Letter);
		Kmer /= m_AS;
		}
	}

uint MerMx::StrToKmer(const string &s) const
	{
	uint k = SIZE(s);
	uint Kmer = 0;
	for (uint i = 0; i < k; ++i)
		{
		Kmer *= m_AS;
		byte c = s[i];
		uint Letter = m_CharToLetter[c];
		asserta(Letter < m_AS);
		Kmer += Letter;
		}
	string Tmp;
	return Kmer;
	}

const char *MerMx::KmerToStr(uint Kmer, uint k, string &s) const
	{
	s.clear();
	for (uint i = 0; i < k; ++i)
		{
		byte Letter = Kmer%m_AS;
		s.push_back(m_LetterToChar[Letter]);
		Kmer /= m_AS;
		}
	reverse(s.begin(), s.end());
	return s.c_str();
	}

int MerMx::GetScoreKmerPair(uint a_Kmer_i, uint a_Kmer_j) const
	{
	uint Kmer_i = a_Kmer_i;
	uint Kmer_j = a_Kmer_j;
	int sum = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		uint code1 = Kmer_i%m_AS;
		uint code2 = Kmer_j%m_AS;
		assert(code1 < 36);
		assert(code2 < 36);
		sum += m_Mx[code1][code2];
		Kmer_i /= m_AS;
		Kmer_j /= m_AS;
		}
	return sum;
	}

uint MerMx::GetMaxLetterCount(uint Kmer) const
	{
	uint Counts[36];
	zero_array(Counts, 36);
	uint maxn = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		uint code = Kmer%m_AS;
		uint n = ++(Counts[code]);
		maxn = max(n, maxn);
		Kmer /= m_AS;
		}
	return maxn;
	}
	
short MerMx::GetScore2merPair(uint a_Kmer_i, uint a_Kmer_j) const
	{
	uint Kmer_i = a_Kmer_i;
	uint Kmer_j = a_Kmer_j;
	assert(Kmer_i < m_AS2);
	assert(Kmer_j < m_AS2);

	uint Letter0_i = Kmer_i%m_AS;
	uint Letter0_j = Kmer_j%m_AS;
	short Score = m_Mx[Letter0_i][Letter0_j];

	uint Letter1_i = Kmer_i/m_AS;
	uint Letter1_j = Kmer_j/m_AS;
	Score += m_Mx[Letter1_i][Letter1_j];

#if DEBUG
	{
	vector<byte> Letters_i;
	vector<byte> Letters_j;
	KmerToLetters(a_Kmer_i, 2, Letters_i);
	KmerToLetters(a_Kmer_j, 2, Letters_j);
	assert(SIZE(Letters_i) == 2);
	assert(SIZE(Letters_j) == 2);
	assert(Letters_i[0] == Letter0_i);
	assert(Letters_j[0] == Letter0_j);
	assert(Letters_i[1] == Letter1_i);
	assert(Letters_j[1] == Letter1_j);
	short ScoreCheck = 0;
	for (uint Pos = 0; Pos < 2; ++Pos)
		ScoreCheck += m_Mx[Letters_i[Pos]][Letters_j[Pos]];
	assert(ScoreCheck == Score);
	}
#endif
	return Score;
	}

int16_t MerMx::GetSelfScore6mer(uint Kmer) const
	{
	assert(Kmer < m_AS_pow[6]);

	uint First3mer = Kmer/m_AS_pow[3];
	uint Second3mer = Kmer%m_AS_pow[3];

	assert(First3mer < m_AS_pow[3]);
	assert(Second3mer < m_AS_pow[3]);

	int FirstScore = m_Mx3[First3mer][First3mer];
	int SecondScore = m_Mx3[Second3mer][Second3mer];
	int Score = FirstScore + SecondScore;
	int16_t Score16 = int16_t(Score);
	assert(int(Score16) == Score);
	return Score16;
	}

int16_t MerMx::GetSelfScore5mer(uint Kmer) const
	{
	int16_t Score16 = 0;
	assert(m_k == 5);
	for (uint i = 0; i < m_k; ++i)
		{
		uint Letter = Kmer%m_AS;
		Score16 += m_Mx[Letter][Letter];
		Kmer /= m_AS;
		}
	return Score16;
	}

int16_t MerMx::GetSelfScoreKmer(uint Kmer) const
	{
	int16_t Score16 = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		uint Letter = Kmer%m_AS;
		Score16 += m_Mx[Letter][Letter];
		Kmer /= m_AS;
		}
#if KMER_COMPLEXITY
// SCOP40
//Max letters [1] = 237643 (12.6%)
//Max letters [2] = 1023722 (54.4%)
//Max letters [3] = 463046 (24.6%)
//Max letters [4] = 112057 (6.0%)
//Max letters [5] = 44543 (2.4%)

// Dictionary
//Max letters [1] = 45239040 (74.8%)
//Max letters [2] = 14779800 (24.4%)
//Max letters [3] = 441000 (0.7%)
//Max letters [4] = 6300 (0.0%)
//Max letters [5] = 36 (0.0%)
	uint8_t n = GetKmerMaxLetterCount(Kmer);
	if (n == 1)
		Score += 16;
	else if (n == 2)
		Score += 8;
	else if (n == 4)
		Score16 -= 8;
	else if (n == 5)
		Score16 -= 16;
#endif
	return Score16;
	}

short MerMx::GetScore3merPair(uint a_Kmer_i, uint a_Kmer_j) const
	{
	uint Kmer_i = a_Kmer_i;
	uint Kmer_j = a_Kmer_j;
	assert(Kmer_i < m_AS3);
	assert(Kmer_j < m_AS3);

	uint Letter0_i = Kmer_i%m_AS;
	uint Letter0_j = Kmer_j%m_AS;
	short Score = m_Mx[Letter0_i][Letter0_j];

	Kmer_i = (Kmer_i - Letter0_i)/m_AS;
	Kmer_j = (Kmer_j - Letter0_j)/m_AS;

	uint Letter1_i = Kmer_i%m_AS;
	uint Letter1_j = Kmer_j%m_AS;
	Score += m_Mx[Letter1_i][Letter1_j];

	Kmer_i = (Kmer_i - Letter1_i)/m_AS;
	Kmer_j = (Kmer_j - Letter1_j)/m_AS;

	assert(Kmer_i < m_AS);
	assert(Kmer_j < m_AS);
	Score += m_Mx[Kmer_i][Kmer_j];

#if DEBUG
	{
	vector<byte> Letters_i;
	vector<byte> Letters_j;
	KmerToLetters(a_Kmer_i, 3, Letters_i);
	KmerToLetters(a_Kmer_j, 3, Letters_j);
	assert(SIZE(Letters_i) == 3);
	assert(SIZE(Letters_j) == 3);

	assert(Letters_i[0] == Letter0_i);
	assert(Letters_j[0] == Letter0_j);
	assert(Letters_i[1] == Letter1_i);
	assert(Letters_j[1] == Letter1_j);
	assert(Letters_i[2] == Kmer_i);
	assert(Letters_j[2] == Kmer_j);

	short ScoreCheck = 0;
	for (uint Pos = 0; Pos < 3; ++Pos)
		ScoreCheck += m_Mx[Letters_i[Pos]][Letters_j[Pos]];
	assert(ScoreCheck == Score);
	}
#endif

	return Score;
	}

void MerMx::Init(const short * const *Mx, uint k, uint AS, uint n)
	{
	asserta(n == 2 || n == 3);
	asserta(k > 1 && k < 8);
	m_k = k;
	if (AS == 20)
		{
		m_CharToLetter = g_CharToLetterAmino;
		m_LetterToChar = g_LetterToCharAmino;
		}
	else if (AS == 36)
		{
		m_CharToLetter = g_CharToLetterMu;
		m_LetterToChar = g_LetterToCharMu;
		}
	else
		Die("MerMx::Init() AS=%d", AS);
	m_Mx = Mx;
	m_AS = AS;
	m_AS2 = AS*AS;
	m_AS3 = AS*AS*AS;

	uint AS_pow = 1;
	m_AS_pow = myalloc(uint, k+1);
	for (uint i = 0; i <= k; ++i)
		{
		m_AS_pow[i] = AS_pow;
		AS_pow *= AS;
		}
	asserta(m_AS2 == m_AS_pow[2]);
	asserta(m_AS3 == m_AS_pow[3]);

	m_Order = myalloc(uint, n == 2 ? m_AS2 : m_AS3);

	m_Mx2 = myalloc(short *, m_AS2);
	m_Scores1 = myalloc(short *, m_AS);
	m_Scores2 = myalloc(short *, m_AS2);

	for (uint Letter = 0; Letter < m_AS; ++Letter)
		m_Scores1[Letter] = myalloc(short, 2*m_AS);

	for (uint Letter = 0; Letter < m_AS; ++Letter)
		BuildRow1(Letter);

	for (uint Kmer = 0; Kmer < m_AS2; ++Kmer)
		{
		m_Mx2[Kmer] = myalloc(short, m_AS2);
		m_Scores2[Kmer] = myalloc(short, 2*m_AS2);
		}

	for (uint Twomer = 0; Twomer < m_AS2; ++Twomer)
		BuildRow2(Twomer);

	if (n == 2)
		return;

	m_Scores3 = myalloc(short *, m_AS3);
	m_Mx3 = myalloc(short *, m_AS3);

	for (uint Threemer = 0; Threemer < m_AS3; ++Threemer)
		{
		ProgressStep(Threemer, m_AS3, "3-mers Mx");
		m_Mx3[Threemer] = myalloc(short, m_AS3);
		m_Scores3[Threemer] = myalloc(short, 2*m_AS3);
		}

	for (uint Threemer = 0; Threemer < m_AS3; ++Threemer)
		{
		ProgressStep(Threemer, m_AS3, "3-mers sorted");
		BuildRow3(Threemer);
		}
	}

void MerMx::BuildRow1(uint Letter_i)
	{
	QuickSortOrderDesc<short>(m_Mx[Letter_i], m_AS, m_Order);
	short LastScore = SHRT_MAX;
	for (uint j = 0; j < m_AS; ++j)
		{
		uint Letter_j = m_Order[j];
		short Score = m_Mx[Letter_i][Letter_j];
		assert(Score <= LastScore);
		m_Scores1[Letter_i][2*j] = Score;
		m_Scores1[Letter_i][2*j+1] = Letter_j;
		LastScore = Score;
		}
	}

void MerMx::BuildRow2(uint Kmer_i)
	{
	vector<short> Scores;
	Scores.reserve(m_AS2);

	for (uint Kmer_j = 0; Kmer_j < m_AS2; ++Kmer_j)
		{
		short Score = GetScore2merPair(Kmer_i, Kmer_j);
		m_Mx2[Kmer_i][Kmer_j] = Score;
		Scores.push_back(Score);
		}

	QuickSortOrderDesc<short>(Scores.data(), m_AS2, m_Order);
	short LastScore = SHRT_MAX;
	for (uint i = 0; i < m_AS2; ++i)
		{
		uint Kmer_j = m_Order[i];
		short Score = Scores[Kmer_j];
		asserta(Score <= LastScore);
		assert(Score == m_Mx2[Kmer_i][Kmer_j]);

		m_Scores2[Kmer_i][2*i] = Score;
		m_Scores2[Kmer_i][2*i+1] = Kmer_j;

		LastScore = Score;
		}
	}

void MerMx::BuildRow3(uint Kmer_i)
	{
	vector<short> Scores;
	Scores.reserve(m_AS3);

	for (uint Kmer_j = 0; Kmer_j < m_AS3; ++Kmer_j)
		{
		short Score = GetScore3merPair(Kmer_i, Kmer_j);
		m_Mx3[Kmer_i][Kmer_j] = Score;
		Scores.push_back(Score);
		}

	QuickSortOrderDesc<short>(Scores.data(), m_AS3, m_Order);
	short LastScore = SHRT_MAX;
	for (uint i = 0; i < m_AS3; ++i)
		{
		uint Kmer_j = m_Order[i];
		short Score = Scores[Kmer_j];
		assert(Score <= LastScore);
		assert(Score == m_Mx3[Kmer_i][Kmer_j]);

		m_Scores3[Kmer_i][2*i] = Score;
		m_Scores3[Kmer_i][2*i+1] = Kmer_j;

		LastScore = Score;
		}
	}

void MerMx::LogMe() const
	{
	Log("AS %u, AS2 %u, AS3 %u\n", m_AS, m_AS2, m_AS3);

	uint m = min(m_AS, 4u);
	for (uint i = 0; i < m; ++i)
		{
		Log("Mx[%c] ", m_LetterToChar[i]);
		for (uint j = 0; j < m; ++j)
			Log(" %4d", m_Mx[i][j]);
		Log(" ... \n");
		}

	string s;
	Log("\n");
	for (uint i = 0; i < 2*m; ++i)
		{
		Log("Mx2[%s] ", KmerToStr(i, 2, s));
		for (uint j = 0; j < 2*m; ++j)
			Log(" %4d", m_Mx2[i][j]);
		Log(" ... \n");
		}

	Log("\n");
	for (uint i = 0; i < 2*m; ++i)
		{
		Log("Mx3[%s] ", KmerToStr(i, 3, s));
		for (uint j = 0; j < 2*m; ++j)
			Log(" %4d", m_Mx3[i][j]);
		Log(" ... \n");
		}

	Log("\n");
	Log("            ");
	for (uint i = 0; i < 2*m; ++i)
		Log("     %s", KmerToStr(i, 3, s));
	Log("\n");
	for (uint i = 0; i < 2*m; ++i)
		{
		Log("Scores2[%s] ", KmerToStr(i, 2, s));
		for (uint j = 0; j < 2*m; ++j)
			{
			uint Score = m_Scores2[i][2*j];
			short Kmer_j = m_Scores2[i][2*j+1];
			Log(" %s=%4d", KmerToStr(Kmer_j, 2, s), Score); 
			}
		Log(" ... \n");
		}
	Log("\n");

	Log("             ");
	for (uint i = 0; i < 2*m; ++i)
		Log("      %s", KmerToStr(i, 3, s));
	Log("\n");
	for (uint i = 0; i < 2*m; ++i)
		{
		Log("Scores3[%s] ", KmerToStr(i, 3, s));
		for (uint j = 0; j < 2*m; ++j)
			{
			uint Score = m_Scores3[i][2*j];
			short Kmer_j = m_Scores3[i][2*j+1];
			Log(" %s=%4d", KmerToStr(Kmer_j, 3, s), Score); 
			}
		Log(" ... \n");
		}
	}

uint MerMx::GetSubmer(uint a_Kmer, uint pos, uint s) const
	{
	string Tmp, Tmp1, Tmp2;
	assert(pos + s <= m_k);
	uint k_minus_pos = m_k - pos;
	uint K_minus_pos_mer = a_Kmer%m_AS_pow[m_k - pos];
	uint Smer = K_minus_pos_mer/m_AS_pow[k_minus_pos - s];
	return Smer;
	}

short MerMx::GetMaxPairScoreSubmer(uint Kmer, uint pos, uint s) const
	{
	uint Smer = GetSubmer(Kmer, pos, s);
	if (s == 1)
		{
		assert(Smer < m_AS);
		return m_Scores1[Smer][0];
		}
	else if (s == 2)
		{
		assert(Smer < m_AS2);
		return m_Scores2[Smer][0];
		}
	else if (s == 3)
		{
		assert(Smer < m_AS3);
		return m_Scores3[Smer][0];
		}
	else
		Die("GetSelfScoreSubmer(s=%u)", s);
	return 0;
	}

uint MerMx::GetHighScoring5mers(uint ABCDE, short MinScore, uint *Fivemers) const
	{
//  Query 5-mer is ABCDE
//  Neigbor 5-mer is abcde
//	Alignment is split into two 2-mers and one letter
//		AB CD E
//		ab cd e
	const uint AB = GetSubmer(ABCDE, 0, 2);
	const uint CD = GetSubmer(ABCDE, 2, 2);
	const uint E = GetSubmer(ABCDE, 4, 1);
#if 0
	{
	string Tmp;
	Log("AB=%s", KmerToStr(AB, 2, Tmp));
	Log(" CD=%s", KmerToStr(CD, 2, Tmp));
	Log(" E=%s", KmerToStr(E, 1, Tmp));
	Log("\n");
	}
#endif

// Max possible score of CD+cd
	const short MaxScore_CD_cd = GetMaxPairScoreSubmer(ABCDE, 2, 2);

// Max possible score of E+e
	const short MaxScore_E_e = m_Scores1[E][0];

// Max possible score for CDE+cde
	const short MaxScore_CDE_cde = MaxScore_CD_cd + MaxScore_E_e;

// Minimum score for AB+ab
	const short MinScore_AB_ab = MinScore - MaxScore_CDE_cde;
#if 0
	Log("MaxScore_CD_cd=%d, MaxScore_E_e=%d, MaxScore_CDE_cde=%d, MinScore_AB_ab=%d\n",
		MaxScore_CD_cd, MaxScore_E_e, MaxScore_CDE_cde, MinScore_AB_ab);
#endif

	const short *Row_AB = m_Scores2[AB];
	const uint ABmul = m_AS_pow[3];
	const uint CDmul = m_AS_pow[1];
	const short *Row_CD = m_Scores2[CD];
	uint n = 0;
	for (uint Idx_ab = 0; Idx_ab < m_AS2; ++Idx_ab)
		{
		short Score_AB_ab = Row_AB[2*Idx_ab];
		if (Score_AB_ab < MinScore_AB_ab)
			break;
		assert(Score_AB_ab + MaxScore_CDE_cde >= MinScore);

		uint ab = Row_AB[2*Idx_ab+1];
		const uint ABmul_ab = ABmul*ab;

		const short MinScore_CD_cd = MinScore - Score_AB_ab - MaxScore_E_e;
		for (uint Idx_cd = 0; Idx_cd < m_AS2; ++Idx_cd)
			{
			short Score_CD_cd = Row_CD[2*Idx_cd];
			if (Score_CD_cd < MinScore_CD_cd)
				break;
			assert(Score_AB_ab + Score_CD_cd + MaxScore_E_e >= MinScore);
			uint cd = Row_CD[2*Idx_cd+1];
			const uint CDmul_cd = CDmul*cd;
			const uint ABmul_ab_plus_CDmul_cd = ABmul_ab + CDmul_cd;
			const short MinScore_E_e = MinScore - Score_AB_ab - Score_CD_cd;
			const short *Scores1E = m_Scores1[E];
			for (uint Idx_e = 0; Idx_e < m_AS; ++Idx_e)
				{
				//short Score_E_e = m_Scores1[E][2*Idx_e];
				short Score_E_e = Scores1E[2*Idx_e];
				if (Score_E_e < MinScore_E_e)
					break;
				assert(Score_AB_ab + Score_CD_cd + Score_E_e >= MinScore);

				//uint e = m_Scores1[E][2*Idx_e+1];
				uint e = Scores1E[2*Idx_e+1];
				// uint abcde = ABmul*ab + CDmul*cd + e;
				uint abcde = ABmul_ab_plus_CDmul_cd + e;
#if 0
				{
				string Tmp;
				Log("ab=%s", KmerToStr(ab, 2, Tmp));
				Log(" cd=%s", KmerToStr(cd, 2, Tmp));
				Log(" e=%s", KmerToStr(e, 1, Tmp));
				Log(" abcde=%s", KmerToStr(abcde, 5, Tmp));
				int Score = GetScoreKmerPair(ABCDE, abcde);
				Log(" Score=%d (Min %d)", Score, MinScore);
				Log("\n");
				}
#endif
#if DEBUG
				{
				int Score = GetScoreKmerPair(ABCDE, abcde);
				assert(Score == Score_AB_ab + Score_CD_cd + Score_E_e);
				assert(Score >= MinScore);
				}
#endif
				Fivemers[n++] = abcde;
				}
			}
		}
	asserta(n <= MAX_HOOD_SIZE);
	return n;
	}

uint MerMx::GetHighScoring6mers(uint Sixmera, short MinScore, uint *Sixmers) const
	{
	const uint First_3mera = GetSubmer(Sixmera, 0, 3);
	const uint Second_3mera = GetSubmer(Sixmera, 3, 3);
	const short MaxScoreSecond_3mer = GetMaxPairScoreSubmer(Sixmera, 3, 3);
	const short MinScoreFirstThreemer = MinScore - MaxScoreSecond_3mer;
	uint n = 0;
	const uint First_mul = m_AS_pow[3];
	const short *Row1 = m_Scores3[First_3mera];
	for (uint i = 0; i < m_AS3; ++i)
		{
		short Score_First3merab = Row1[2*i];
		if (Score_First3merab < MinScoreFirstThreemer)
			break;
		const uint First_3merb = Row1[2*i+1];
		const short MinScore_Second3mer = MinScore - Score_First3merab;
		const short *Row2 = m_Scores3[Second_3mera];
		for (uint j = 0; j < m_AS3; ++j)
			{
			short Score_Second3merab = Row2[2*j];
			if (Score_Second3merab < MinScore_Second3mer)
				break;
			uint Second_3merb = Row2[2*j+1];
			uint ThisSixmer = First_3merb*First_mul + Second_3merb;
			Sixmers[n++] = ThisSixmer;
			}
		}
	return n;
	}

uint MerMx::GetHighScoringKmers(uint Kmer, short MinScore, uint *Kmers) const
	{
	g_Lock.lock();
	++g_HoodCallCount;
	g_Lock.unlock();
	switch (m_k)
		{
	case 4: return GetHighScoring4mers(Kmer, MinScore, Kmers);
	case 5: return GetHighScoring5mers(Kmer, MinScore, Kmers);
	case 6: return GetHighScoring6mers(Kmer, MinScore, Kmers);
		}
	asserta(false);
	return UINT_MAX;
	}

uint MerMx::GetHighScoring4mers(uint Fourmera, short MinScore, uint *Fourmers) const
	{
	const uint First_2mera = GetSubmer(Fourmera, 0, 2);
	const uint Second_2mera = GetSubmer(Fourmera, 2, 2);
	const short MaxScoreSecond_2mer = GetMaxPairScoreSubmer(Fourmera, 2, 2);
	const short MinScoreFirstThreemer = MinScore - MaxScoreSecond_2mer;
	uint n = 0;
	const uint First_mul = m_AS_pow[2];
	const short *Row1 = m_Scores2[First_2mera];
	for (uint i = 0; i < m_AS2; ++i)
		{
		short Score_First2merab = Row1[2*i];
		if (Score_First2merab < MinScoreFirstThreemer)
			break;
		const uint First_2merb = Row1[2*i+1];
		const short MinScore_Second2mer = MinScore - Score_First2merab;
		const short *Row2 = m_Scores2[Second_2mera];
		for (uint j = 0; j < m_AS2; ++j)
			{
			short Score_Second2merab = Row2[2*j];
			if (Score_Second2merab < MinScore_Second2mer)
				break;
			uint Second_2merb = Row2[2*j+1];
			uint ThisFourmer = First_2merb*First_mul + Second_2merb;
			Fourmers[n++] = ThisFourmer;
			}
		}
	return n;
	}

uint MerMx::GetHighScoring5mers_Brute(uint Fivemer, short MinScore, uint *Fivemers,
									  bool Trace) const
	{
	uint n = 0;
	uint DictSize = m_AS_pow[5];
	for (uint Fivemer2 = 0; Fivemer2 < DictSize; ++Fivemer2)
		{
		short Score = GetScoreKmerPair(Fivemer, Fivemer2);
		if (Score >= MinScore)
			{
			if (Trace)
				{
				string Tmp;
				Log(" %s", KmerToStr(Fivemer2, 5, Tmp));
				}
			Fivemers[n++] = Fivemer2;
			}
		}
	if (Trace)
		Log("\n");
	return n;
	}

uint MerMx::GetHighScoring6mers_Brute(uint Sixmer, short MinScore, uint *Sixmers,
									  bool Trace) const
	{
	uint n = 0;
	uint DictSize = m_AS_pow[6];
	for (uint Sixmer2 = 0; Sixmer2 < DictSize; ++Sixmer2)
		{
		short Score = GetScoreKmerPair(Sixmer, Sixmer2);
		if (Score >= MinScore)
			{
			if (Trace)
				{
				string Tmp;
				Log(" %s", KmerToStr(Sixmer2, 6, Tmp));
				}
			Sixmers[n++] = Sixmer2;
			}
		}
	if (Trace)
		Log("\n");
	return n;
	}

int16_t *MerMx::BuildSelfScores6mers() const
	{
	const uint AS6 = m_AS_pow[6];
	int16_t *SelfScores = myalloc(int16_t, AS6);
	for (uint Kmer = 0; Kmer < AS6; ++Kmer)
		SelfScores[Kmer] = GetSelfScore6mer(Kmer);
	return SelfScores;
	}

int16_t *MerMx::BuildSelfScores5mers() const
	{
	const uint AS5 = m_AS_pow[5];
	int16_t *SelfScores = myalloc(int16_t, AS5);
	for (uint Kmer = 0; Kmer < AS5; ++Kmer)
		SelfScores[Kmer] = GetSelfScore5mer(Kmer);
	return SelfScores;
	}

int16_t *MerMx::BuildSelfScores_Kmers() const
	{
	const uint ASk = m_AS_pow[m_k];
	int16_t *SelfScores = myalloc(int16_t, ASk);
	for (uint Kmer = 0; Kmer < ASk; ++Kmer)
		SelfScores[Kmer] = GetSelfScoreKmer(Kmer);
	return SelfScores;
	}

void LogHoodCount()
	{
	ProgressLog("%s hood calls\n", Int64ToStr(g_HoodCallCount));
	}

#define LOW_COMPLEXITY 0

/***
     60.5M  DICT_SIZE
       92G  Total size of all neighborhoods
      1.2M  Kmers with low self score (1.9%)
     41.3k  Max size 'BBBBB' (41293)
      1521  Mean
       690  Median
***/
void cmd_kmrnbh()
	{
	const MerMx &ScoreMx = GetMuMerMx(5);
	uint *Kmers = myalloc(uint, DICT_SIZE);
	uint MinScore = MIN_KMER_PAIR_SCORE;
	uint64 Sumn = 0;
	uint Maxn = 0;
	uint MaxKmer = UINT_MAX;
	const uint N = DICT_SIZE;
	vector<float> Sizes;
	Sizes.reserve(DICT_SIZE);
	uint M = 0;
	uint LowSelfScore = 0;
#if LOW_COMPLEXITY
	uint LowComplexity = 0;
#endif
	for (uint Kmer = 0; Kmer < N; ++Kmer)
		{
		ProgressStep(Kmer, DICT_SIZE, "Neighborhood");
#if LOW_COMPLEXITY
		uint mlc = ScoreMx.GetMaxLetterCount(Kmer);
		if (mlc >= 3)
			{
			++LowComplexity;
			continue;
			}
#endif
		uint n = ScoreMx.GetHighScoring5mers(Kmer, MinScore, Kmers);
		if (n == 0)
			{
			short SelfScore = ScoreMx.GetScoreKmerPair(Kmer, Kmer);
			asserta(SelfScore < MIN_KMER_PAIR_SCORE);
			++LowSelfScore;
			continue;
			}
		++M;
		if (n > Maxn)
			{
			Maxn = n;
			MaxKmer = Kmer;
			}
		Sumn += n;
		Sizes.push_back(float(n));
		}

	QuartsFloat Q;
	GetQuartsFloat(Sizes, Q);

	string TmpStr;
	ProgressLog("%10.10s  DICT_SIZE\n", IntToStr(DICT_SIZE));

	ProgressLog("%10.10s  Total size of all neighborhoods\n",
				Int64ToStr(Sumn));

#if LOW_COMPLEXITY
	ProgressLog("%10.10s  Kmers with low complexity (%.1f%%)\n",
				IntToStr(LowComplexity), GetPct(LowComplexity, N));
#endif

	ProgressLog("%10.10s  Kmers with low self score (%.1f%%)\n",
				IntToStr(LowSelfScore), GetPct(LowSelfScore, N));

	ProgressLog("%10.10s  Max size '%s' (%u)\n", 
				IntToStr(Maxn),
				ScoreMx.KmerToStr(MaxKmer, 5, TmpStr),
				Maxn);

	ProgressLog("%10.10s  Mean\n", IntToStr(uint(Q.Avg)));
	ProgressLog("%10.10s  Median\n", IntToStr(uint(Q.Med)));
	}
