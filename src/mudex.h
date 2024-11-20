#pragma once

class SeqDB;

// Mu 5-mer index
class MuDex
	{
public:
	static const uint32_t m_DictSize = 36*36*36*36*36;
	static const uint32_t m_ItemSize = 6;	// 4 byte SeqIdx + 2 byte Pos
	static const uint32_t m_k = 5;

/***
After Pass 1, before adjustment:
	Size(Kmer)
		= m_Finger[Kmer+1]

After adjustment, before Pass 2:
	Data offset of first entry
		= m_Finger[Kmer+1]

	Size(Kmer)
		= m_Finger[Kmer+2] - m_Finger[Kmer+1]

After Pass 2:
	Data offset of first entry:
		m_Finger[Kmer]

	Size(Kmer)
		= m_Finger[Kmer+1] - m_Finger[Kmer]

***/
#if DEBUG
	vector<uint> m_KmerToCount;
	vector<uint> m_KmerToDataOffset;
#endif

	SeqDB *m_SeqDB = 0;
	uint32_t m_Size = 0;
	uint32_t *m_Finger = 0;
	uint8_t *m_Data = 0;

public:
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

	uint GetRowSize(uint Kmer) const
		{
		assert(Kmer < m_DictSize);
		uint n = m_Finger[Kmer+1] - m_Finger[Kmer];
		assert(m_Finger[Kmer] + n <= m_Size);
		return n;
		}

	void Put(uint DataOffset, uint32_t SeqIdx, uint16_t SeqPos)
		{
		uint8_t *ptr = m_Data + DataOffset;
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
		const uint8_t *ptr = m_Data + DataOffset;
		SeqIdx = *(uint32_t *) ptr;
		SeqPos = *(uint16_t *) (ptr + 4);
		}
	};
