#pragma once

class SeqDB;

// Mu k-mer index
class MuDex
	{
public:
	static const uint32_t m_ItemSize = 6;	// 4 byte SeqIdx + 2 byte Pos
	static uint32_t m_DictSize; // 36*36*36*36*36;
	static uint32_t m_k; // 5;

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

public:
	static void Set_k(uint k);

public:
	void FromSeqDB(const SeqDB &Input);
	const char *KmerToStr(uint Kmer, string &s) const;
	uint StrToKmer(const string &s) const;
	uint StrToKmer(const char *s) const;
	void Alloc_Pass1();
	void AdjustFinger();
	void Alloc_Pass2();
	void AddSeq_Pass1(uint SeqIdx, const char *Label, const char *Seq, uint L);
	void AddSeq_Pass2(uint SeqIdx, const char *Label, const char *Seq, uint L);
	void LogStats() const;
	void Validate() const;
	void ValidateKmer(uint Kmer) const;
	uint GetSeqKmer(uint SeqIdx, uint SeqPos) const;
	void LogIndexKmer(uint Kmer) const;
	const uint32_t GetDataOffset(uint Kmer) const { return m_Finger[Kmer]; }
	const uint32_t GetRowStart(uint Kmer) const
		{
		assert(Kmer < m_DictSize);
		return m_Finger[Kmer];
		}

#if DEBUG
	void CheckAfterPass1() const;
	void CheckAfterAdjust() const;
	void CheckAfterPass2() const;
#endif

	uint GetRowSize(uint Kmer) const
		{
		assert(Kmer < m_DictSize);
		uint n = m_Finger[Kmer+1] - m_Finger[Kmer];
		assert(m_Finger[Kmer] + n <= m_Size);
		return n;
		}

	void Put(uint DataOffset, uint32_t SeqIdx, uint16_t SeqPos)
		{
		uint8_t *ptr = m_Data + m_ItemSize*DataOffset;
		*(uint32_t *) ptr = SeqIdx;
		*(uint16_t *) (ptr + 4) = SeqPos;
#if DEBUG
		{
		uint32_t Check_SeqIdx;
		uint16_t Check_SeqPos;
		Get(DataOffset, Check_SeqIdx, Check_SeqPos);
		assert(Check_SeqIdx == SeqIdx);
		assert(Check_SeqPos == SeqPos);
		}
#endif
		}

	void Get(uint DataOffset, uint32_t &SeqIdx, uint16_t &SeqPos) const
		{
		const uint8_t *ptr = m_Data + m_ItemSize*DataOffset;
		SeqIdx = *(uint32_t *) ptr;
		SeqPos = *(uint16_t *) (ptr + 4);
		}
	};
