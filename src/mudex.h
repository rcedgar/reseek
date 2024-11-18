#pragma once

// Mu 5-mer index
class MuDex
	{
public:
	static const uint32_t m_DictSize = 36*36*36*36*36;
	uint32_t *m_RowSizes = 0;
	uint32_t *m_RowSizes2 = 0;
	uint32_t **m_SeqIdxRows = 0;
	uint16_t **m_SeqPosRows = 0;

public:
	const char *KmerToStr(uint Kmer, string &s) const;
	uint StrToKmer(const string &s) const;
	void Alloc_Pass1();
	void Alloc_Pass2();
	void AddSeq_Pass1(uint SeqIdx, const char *Label, const char *Seq, uint L);
	void AddSeq_Pass2(uint SeqIdx, const char *Label, const char *Seq, uint L);
	void LogStats() const;
	};
