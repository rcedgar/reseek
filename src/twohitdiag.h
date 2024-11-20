#pragma once

#include <list>

/***
Radix
	= seqidx%1024, diag%16
	= 10 bits + 14 bits, 2^14
	= 16k headers

Item 
	= Seq+Pos+Diag
	= 32+16+16 bytes

m_Data has fixed space for 16 entries per radix.

Memory per thread
	m_BusyRdxs = 4*16k bytes = 65536 bytes
	m_RdxSizes = 4*16k bytes = 65536 bytes
	m_Data = 2097152 (2 M)

Overflow if >16 entries in a radix:
	List of int8_t* is allocated for this radix
	Each list entry points to 16-entry data
	All overflow lists & data are freed and 
	  deleted in Reset()
***/

class TwoHitDiag
	{
public:
	typedef list<uint8_t *> *listptr_t;

	static const uint32_t m_SeqMod = 1024;
	static const uint32_t m_DiagMod = 16;
	static const uint32_t m_NrRdxs = m_SeqMod*m_DiagMod; // 16384
	static const uint32_t m_SeqMask = m_NrRdxs - 1;
	static const uint32_t m_DiagShift = 10;
	static const uint32_t m_FixedEntriesPerRdx = 16;
	static const uint32_t m_OverflowPointerBytes = sizeof(uint8_t *);
	static const uint32_t m_ItemSize = 8;
	static const uint32_t m_FixedBytes = 
		m_NrRdxs*m_FixedEntriesPerRdx*m_ItemSize; // 2097152 (2 M)

// Total number of items, for validation only
	uint m_Size = 0;

// Arrays of size m_NrRdxs with index Rdx
	list<uint8_t *> **m_Overflows = 0;
	uint32_t *m_Sizes = 0;

// Vector of size m_BusyRdxCount <= m_NrRdxs
	uint32_t *m_BusyRdxs = 0;
	uint32_t m_BusyCount = 0;

// Fixed buffer size m_FixedBytes
	uint8_t *m_Data = 0;

public:
	TwoHitDiag();
	~TwoHitDiag();

	uint32_t UnpackRdx(uint32_t Rdx, uint16_t &Diag) const
		{
		Diag = (Rdx >> m_DiagShift);
		uint32_t SeqIdx = (Rdx & m_SeqMask);
		return SeqIdx;
		}

	uint32_t GetRdx(uint32_t SeqIdx, uint16_t Diag) const
		{
		uint32_t Rdx = SeqIdx%m_SeqMod + ((Diag%m_DiagMod) << m_DiagShift);
		assert(Rdx < m_NrRdxs);
		return Rdx;
		}

	void PutRaw(uint8_t *ptrData, uint32_t SeqIdx, uint16_t Pos, uint16_t Diag)
		{
		uint64_t u = uint64_t(SeqIdx)
			| (uint64_t(Pos) << 32)
			| (uint64_t(Diag) << 48);
		uint64_t *ptrData64 = (uint64_t *) ptrData;
		*ptrData64 = u;
		}

	uint32_t GetRaw(const uint8_t *ptrData, uint16_t &Pos, uint16_t &Diag) const
		{
		uint64_t *ptrData64 = (uint64_t *) ptrData;
		uint64 u = *ptrData64;
		uint32_t SeqIdx = uint32_t(u);
		uint32_t v = (u >> 32);
		Pos = uint16_t(v);
		Diag = (v >> 16);
		return SeqIdx;
		}

	void Add(uint32_t SeqIdx, uint16_t Pos, uint16_t Diag);
	void Reset();
	void ValidateEmpty() const;
	void ValidateEmptyRdx(uint Rdx) const;
	void ValidateRdx(uint Rdx, uint MaxSeqIdx, uint MaxPos, uint MaxDiag) const;
	void Validate(uint MaxSeqIdx, uint MaxPos, uint MaxDiag) const;
	void LogRdx(uint Rdx) const;
	void LogStats() const;
	};
