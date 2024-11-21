#pragma once

#include <list>

/***
Radix
	= 10 lo bits of SeqIdx, 4 lo bits of Diag
	= SeqIdx%1024, Diag%16
	= 10 bits + 14 bits
	= 2^14 radixes = 16k

SeqIdx
	32 bits, SeqIdx/1024 is 32-10 = 22 bits
	Max hi bits = 2^22 = 4M, so max 4M SeqIdx's per radix
	Bit map is 4M/8 = 0.5Mb
	Could search for duplicates with two bitmaps = 1Mb per thread

Diag 
	14 bits (max 16k), Diag/16 = 14-4 = 10 bits

Item 
	= SeqIdx/1024,Diag/16
	= 22 hi bits of SeqIdx, 10 hi bits of Diag
	= 32 bits

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
	static const uint32_t m_Mask10 = 0b1111111111;
	static const uint32_t m_Mask14 = 0b11111111111111;
	static const uint32_t m_Mask22 = 0b1111111111111111111111;

	static const uint32_t m_DiagHi = m_Mask14 + 1; // MaxDiag + 1 = 2^14 = 16k
	static const uint32_t m_SeqMod = 1024;
	static const uint32_t m_DiagMod = 16;
	static const uint32_t m_NrRdxs = m_SeqMod*m_DiagMod; // 16384
	static const uint32_t m_DiagShiftRdx = 10;
	static const uint32_t m_DiagShiftItem = 22;
	static const uint32_t m_FixedItemsPerRdx = 16;
	static const uint32_t m_ItemSize = sizeof(uint32_t);
	static const uint32_t m_TotalFixedItems = m_NrRdxs*m_FixedItemsPerRdx;

	static const uint32_t m_SeqIdxMaskRdx =	m_Mask10;
	static const uint32_t m_SeqIdxMaskItem = m_Mask22;

// Total number of items, for validation only
	uint m_Size = 0;

// Arrays of size m_NrRdxs with index Rdx
	list<uint32_t *> **m_Overflows = 0;
	uint32_t *m_Sizes = 0;

// Vector of size m_BusyCount <= m_NrRdxs
	uint32_t *m_BusyRdxs = 0;
	uint32_t m_BusyCount = 0;

// Fixed buffer size m_FixedBytes
	uint32_t *m_Data = 0;

public:
	TwoHitDiag();
	~TwoHitDiag();

	uint32_t UnpackRdx(uint32_t Rdx, uint16_t &DiagLoBits) const
		{
		DiagLoBits = (Rdx >> m_DiagShiftRdx);
		uint32_t SeqIdxLoBits = (Rdx & m_SeqIdxMaskRdx);
		return SeqIdxLoBits;
		}

	uint32_t GetRdx(uint32_t SeqIdx, uint16_t Diag) const
		{
		asserta(Diag <= m_Mask14);

		uint32_t SeqIdxLoBits = SeqIdx%m_SeqMod;
		uint16_t DiagLoBits = Diag%m_DiagMod;
		uint32_t Rdx = (SeqIdxLoBits | (DiagLoBits << m_DiagShiftRdx));
#if DEBUG
		assert(Rdx < m_NrRdxs);
		uint16_t DiagLoBits2;
		uint32_t SeqIdxLoBits2 = UnpackRdx(Rdx, DiagLoBits2);
		assert(SeqIdxLoBits2 == SeqIdxLoBits);
		assert(DiagLoBits2 == DiagLoBits);
#endif
		return Rdx;
		}

	void PutRaw(uint32_t *ptrData, uint32_t SeqIdx, uint16_t Diag)
		{
		asserta(Diag <= m_Mask14);

		uint32_t Item = ((SeqIdx/m_SeqMod) | (Diag/m_DiagMod) << m_DiagShiftItem);
		*ptrData = Item;
		}

	uint32_t GetHiBits(const uint32_t *ptrData, uint16_t &DiagHiBits) const
		{
		uint32_t Item = *ptrData;
		uint32_t SeqIdxHiBits = (Item & m_SeqIdxMaskItem);
		DiagHiBits = uint16_t(Item >> m_DiagShiftItem);
		return SeqIdxHiBits;
		}

	uint32_t GetAllBits(uint Rdx, const uint32_t *ptrData, uint16_t &Diag) const
		{
		uint32_t Item = *ptrData;
		uint32_t SeqIdxHiBits = (Item & m_SeqIdxMaskItem);
		uint16_t DiagHiBits = uint16_t(Item >> m_DiagShiftItem);

		uint16_t DiagLoBits;
		uint32_t SeqIdxLoBits = UnpackRdx(Rdx, DiagLoBits);

		uint32_t SeqIdx = (SeqIdxLoBits | (SeqIdxHiBits*m_SeqMod));

		Diag = (DiagLoBits | (DiagHiBits*m_DiagMod));
		asserta(Diag <= m_Mask14);

		return SeqIdx;
		}

	void TestPutGet(uint32_t SeqIdx, uint16_t Diag) const;

	void Add(uint32_t SeqIdx, uint16_t Diag);
	void Reset();
	void ValidateEmpty() const;
	void ValidateEmptyRdx(uint Rdx) const;
	void ValidateRdx(uint Rdx, uint MaxSeqIdx, uint MaxDiag) const;
	void Validate(uint MaxSeqIdx, uint MaxDiag) const;
	void LogRdx(uint Rdx) const;
	void LogStats() const;
	void GetAll(set<pair<uint32_t, uint16_t> > &SeqIdxDiagPairs) const;
	void AppendAll(uint Rdx, set<pair<uint32_t, uint16_t> > &SeqIdxDiagPairs) const;
	};
