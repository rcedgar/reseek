#pragma once

// K-mer scoring matrix
// First letter in K-mer is in most significant bits = Kmer/(AS^(k-1))
// Last letter in K-mer is in least significant bits = Kmer%AS
class MerMx
	{
public:
// K-mer size
	uint m_k = 6;

// AS = alphabet size
	uint m_AS = UINT_MAX;
	uint m_AS2 = UINT_MAX;	// AS^2
	uint m_AS3 = UINT_MAX;	// AS^3
	uint *m_AS_pow = 0;		// m_AS_pow[i] == AS^i
	const byte *m_CharToLetter = 0;
	const byte *m_LetterToChar = 0;

// AS x AS scoring matrix m_Mx[Letter_i][Letter_j]
	const short * const *m_Mx = 0;

// AS2 x AS2 scoring matrix m_Mx2[Kmer_i][Kmer_j]
	short **m_Mx2 = 0;

// AS3 x AS3 scoring matrix m_Mx3[Kmer_i][Kmer_j]
	short **m_Mx3 = 0;

// m_Scores1 sorted by decreasing score
// m_Scores1[Letter_j][2*i] = score of Letter_i,Letter_j
// m_Scores1[Letter_j][2*i+1] = Letter_j
	short **m_Scores1 = 0;

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
	void Init(const short * const *Mx, uint k, uint AS, uint n);
	void BuildRow1(uint Letter);
	void BuildRow2(uint Twomer);
	void BuildRow3(uint Threemer);
	int GetScoreKmerPair(uint Kmer_i, uint Kmer_j) const;
	short GetScore2merPair(uint Kmer_i, uint Kmer_j) const;
	short GetScore3merPair(uint Kmer_i, uint Kmer_j) const;
	void KmerToLetters(uint Kmer, uint k, vector<byte> &Letters) const;
	const char *KmerToStr(uint Kmer, uint k, string &s) const;
	uint StrToKmer(const string &s) const;
	void LogMe() const;

// s-mer starting at position i in Kmer
	uint GetSubmer(uint Kmer, uint s, uint i) const;

// Max score of s-mer starting at position pos in Kmer
//   against another s-mer
	short GetMaxPairScoreSubmer(uint Kmer, uint pos, uint s) const;

	uint GetHighScoring5mers(uint Kmer, short MinScore, uint *Kmers) const;
	uint GetHighScoring6mers(uint Kmer, short MinScore, uint *Kmers) const;
	uint GetHighScoring5mers_Brute(uint Kmer, short MinScore, uint *Kmers,
								   bool Trace = false) const;
	uint GetHighScoring6mers_Brute(uint Kmer, short MinScore, uint *Kmers,
								   bool Trace = false) const;
	};
