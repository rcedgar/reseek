#pragma once

// K-mer scoring matrix
class MerMx
	{
public:
// AS = alphabet size
	uint m_AS = UINT_MAX;
	uint m_AS2 = UINT_MAX;	// AS^2
	uint m_AS3 = UINT_MAX;	// AS^3

// AS x AS scoring matrix m_Mx[Letter_i][Letter_j]
	const short * const *m_Mx = 0;

// AS2 x AS2 scoring matrix m_Mx2[Kmer_i][Kmer_j]
	short **m_Mx2 = 0;

// AS3 x AS3 scoring matrix m_Mx3[Kmer_i][Kmer_j]
	short **m_Mx3 = 0;

// m_Scores2 sorted by decreasing score
// m_Scores2[2mer_i][2*i] = score of 2mer_i,2mer_j
// m_Scores2[2mer_i][2*i+1] = 2mer_j
	short **m_Scores2 = 0;

// m_Scores3 sorted by decreasing score
// m_Scores3[3mer_i][2*i] = score of 3mer_i,3mer_j
// m_Scores3[3mer_i][2*i+1] = 3mer_j
	short **m_Scores3 = 0;

	uint *m_Order = 0;

public:
	void Init(short **Mx, uint AS);
	void BuildRow2(uint Letter);
	void BuildRow3(uint Letter);
	short GetScore2(uint Kmer_i, uint Kmer_j) const;
	short GetScore3(uint Kmer_i, uint Kmer_j) const;
	void KmerToLetters(uint Kmer, uint k, vector<byte> &Letters) const;
	const char *KmerToStr(uint Kmer, uint k, string &s) const;
	void LogMe() const;
	};
