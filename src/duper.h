#pragma once

/***
Find duplicates in list of 32-bit integers != UINT_MAX
  provided serially through calls to Add().

Number of integers is known in advance.

Implemented as hash table with two vectors:
	Input integers where UINT_MAX indicates empty.
	Bitmap where 1=seen before
***/

class Duper
	{
public:
	uint32_t m_InputSize = 0;
	uint32_t m_TableSize = 0;
	uint32_t *m_Ints = 0;
	uint8_t *m_Bits = 0;
	uint32_t *m_Dupes = 0;
	uint32_t m_DupeCount = 0;

public:
	Duper(uint32_t InputSize);
	~Duper();

public:
	uint32_t Hash(uint32_t i) const
		{
		return i%m_TableSize;
		}

	void SetBit(uint Idx)
		{
		assert(Idx < m_TableSize);
		assert(GetBit(Idx) == 0);
		m_Bits[Idx/8] |= (1 << uint8_t(Idx%8));
		}

	int GetBit(uint Idx) const
		{
		return m_Bits[Idx/8] & (1 << uint8_t(Idx%8));
		}

	void FoundDupe(uint32_t i)
		{
		assert(m_DupeCount < m_InputSize);
		m_Dupes[m_DupeCount++] = i;
		}

	void Add(uint32_t i);
	};
