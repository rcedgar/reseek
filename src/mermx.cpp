#include "myutils.h"
#include "mermx.h"
#include "alpha.h"
#include "sort.h"

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
		uint Letter = g_CharToLetterAmino[c];
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
		s.push_back(g_LetterToCharAmino[Letter]);
		Kmer /= m_AS;
		}
	reverse(s.begin(), s.end());
	return s.c_str();
	}

short MerMx::GetScoreKmerPair(uint a_Kmer_i, uint a_Kmer_j) const
	{
	uint Kmer_i = a_Kmer_i;
	uint Kmer_j = a_Kmer_j;
	int sum = 0;
	for (uint i = 0; i < 6; ++i)
		{
		uint code1 = Kmer_i%20;
		uint code2 = Kmer_j%20;
		sum += m_Mx[code1][code2];
		Kmer_i /= 20;
		Kmer_j /= 20;
		}
	return sum;
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

void MerMx::Init(short **Mx, uint k, uint AS)
	{
	asserta(k > 1 && k < 8);
	m_k = k;
	asserta(AS > 2 && AS < 4096);
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

	m_Order = myalloc(uint, m_AS3);

	m_Mx2 = myalloc(short *, m_AS2);
	m_Mx3 = myalloc(short *, m_AS3);
	m_Scores2 = myalloc(short *, m_AS2);
	m_Scores3 = myalloc(short *, m_AS3);

	for (uint Kmer = 0; Kmer < m_AS2; ++Kmer)
		{
		m_Mx2[Kmer] = myalloc(short, m_AS2);
		m_Scores2[Kmer] = myalloc(short, 2*m_AS2);
		}

	for (uint Kmer = 0; Kmer < m_AS3; ++Kmer)
		{
		ProgressStep(Kmer, m_AS3, "3-mers Mx");
		m_Mx3[Kmer] = myalloc(short, m_AS3);
		m_Scores3[Kmer] = myalloc(short, 2*m_AS3);
		}

	for (uint Kmer = 0; Kmer < m_AS2; ++Kmer)
		BuildRow2(Kmer);

	for (uint Kmer = 0; Kmer < m_AS3; ++Kmer)
		{
		ProgressStep(Kmer, m_AS3, "3-mers sorted");
		BuildRow3(Kmer);
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
		Log("Mx[%c] ", g_LetterToCharAmino[i]);
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
	if (s == 2)
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

// Work buffer stores possible first 3mers of high-scoring Sixmer
//	Work[2*i] is score of first 3mer against first 3mer in Sixmer
//	Work[2*i+1] is i'th first 3mer
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
#if DEBUG
			{
			short CheckScore = GetScoreKmerPair(Sixmer, ThisSixmer);
			assert(CheckScore == Score);
			}
#endif
			Sixmers[n++] = ThisSixmer;
			}
		}
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
