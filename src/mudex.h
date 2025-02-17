#pragma once

class SeqDB;

// Mu 5-mer index
class MuDex
	{
public:
	static const uint32_t m_DictSize;
	static const uint32_t m_ItemSize;
	static const uint32_t m_k;
	static const uint32_t m_K;
	static const char m_Pattern[6];
	static const byte m_Offsets[5];

/***
Finger index.

Pass 1 counts K-mers.
After Pass 1, before adjustment:
	Size(Kmer)
		= m_Finger[Kmer+1]

After adjustment, before Pass 2:
	Data offset of first entry
		= m_Finger[Kmer+1]

	Size(Kmer)
		= m_Finger[Kmer+2] - m_Finger[Kmer+1]

Pass 2 builds the index.
During Pass 2:
	Initially m_Finger[Kmer+1] is data offset of first entry
	With each K-mer added, m_Finger[Kmer+1] is incremented
	Finally m_Finger[Kmer+1] is 1 + (data offset of last entry of Kmer)
	   = (data offset of first entry for Kmer+1, if any).
	   For this to work for the special case Kmer=0, an initial 0 in
	   m_Finger[0] is preserved and not touched throughout Passes 1 and 2.

After Pass 2:
	Data offset of first entry:
		m_Finger[Kmer]

	Size(Kmer)
		= m_Finger[Kmer+1] - m_Finger[Kmer]

***/
#if DEBUG
	vector<uint> m_KmerToCount1;
	vector<uint> m_KmerToCount2;
	vector<uint> m_KmerToDataStart;
#endif

	const SeqDB *m_SeqDB = 0;
	uint32_t m_Size = 0;
	uint32_t *m_Finger = 0;
	uint8_t *m_Data = 0;
	int16_t *m_KmerSelfScores = 0;
	int m_MinKmerSelfScore = 0;

public:
	void FromSeqDB(const SeqDB &Input);
	const char *KmerToStr(uint Kmer, string &s) const;
	uint StrToKmer(const string &s) const;
	uint BytesToKmer(const byte *s) const;
	uint StrToKmer(const char *s) const;
	void Alloc_Pass1();
	void AdjustFinger();
	void Alloc_Pass2();
	void AddSeq_Pass1(uint SeqIdx, const char *Label, const char *Seq, uint L);
	void AddSeq_Pass2(uint SeqIdx, const char *Label, const char *Seq, uint L);
	void LogStats() const;
	void Validate() const;
	void ValidateKmer(uint Kmer) const;
	uint GetSeqKmer(const byte *Seq, uint SeqPos, bool SelfScoreMask) const;
	void LogIndexKmer(uint Kmer) const;
	const uint32_t GetRowStart(uint Kmer) const
		{
		assert(Kmer < m_DictSize);
		return m_Finger[Kmer];
		}

	uint GetRowSize(uint Kmer) const;
	void Put(uint DataOffset, uint32_t SeqIdx, uint16_t SeqPos);
	void Get(uint DataOffset, uint32_t &SeqIdx, uint16_t &SeqPos) const;
	void GetKmers(const byte *Seq, uint L, vector<uint> &Kmers) const;

#if DEBUG
	void CheckAfterPass1() const;
	void CheckAfterAdjust() const;
	void CheckAfterPass2() const;
#endif
	};
