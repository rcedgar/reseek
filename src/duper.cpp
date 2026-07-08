#include "myutils.h"
#include "duper.h"

static void FreeDuperTables(Duper &D)
	{
	myfree(D.m_Ints);
	myfree(D.m_Bits);
	myfree(D.m_Dupes);
	D.m_Ints = 0;
	D.m_Bits = 0;
	D.m_Dupes = 0;
	D.m_TableSize = 0;
	D.m_AllocInputSize = 0;
	}

static void AllocDuperTables(Duper &D, uint InputSize)
	{
	asserta(InputSize < UINT_MAX/4);
	D.m_InputSize = InputSize;
	D.m_AllocInputSize = InputSize;

	uint FindPrime(uint Min, uint Max);
	uint N = InputSize + 7;
	uint Lo = N*3;
	uint Hi = Lo + Lo/4;
	D.m_TableSize = FindPrime(Lo, Hi);
	D.m_Ints = myalloc(uint32_t, D.m_TableSize);
	memset(D.m_Ints, 0xff, D.m_TableSize*sizeof(uint32_t));

	uint BitBytes = (D.m_TableSize + 7)/8;
	D.m_Bits = myalloc(uint8_t, BitBytes);
	memset(D.m_Bits, 0, BitBytes);

	D.m_Dupes = myalloc(uint32_t, InputSize);
	D.m_DupeCount = 0;
	}

Duper::Duper(uint InputSize)
	{
	AllocDuperTables(*this, InputSize);
	}

Duper::~Duper()
	{
	FreeDuperTables(*this);
	}

void Duper::Init(uint InputSize)
	{
	m_DupeCount = 0;
	m_InputSize = InputSize;
	if (InputSize <= m_AllocInputSize && m_Ints != 0)
		{
		memset(m_Ints, 0xff, m_TableSize*sizeof(uint32_t));
		uint BitBytes = (m_TableSize + 7)/8;
		memset(m_Bits, 0, BitBytes);
		return;
		}
	FreeDuperTables(*this);
	AllocDuperTables(*this, InputSize);
	}

void Duper::Add(uint32_t i)
	{
	uint32_t h = Hash(i);
	for (uint k = 0; k < m_TableSize; ++k)
		{
		uint32_t Idx = (h + k)%m_TableSize;
		if (m_Ints[Idx] == UINT_MAX)
			{
			m_Ints[Idx] = i;
			return;
			}
		else if (m_Ints[Idx] == i)
			{
			if (GetBit(Idx))
				return;
			SetBit(Idx);
			FoundDupe(i);
			return;
			}
		}
	Die("Duper::Add() overflow");
	}

static uint Test()
	{
	vector<uint32_t> v;
	map<uint32_t, uint32_t> countmap;
	uint N = 7 + randu32()%100024;
	uint DupeCount = 0;
	for (uint k = 0; k < N; ++k)
		{
		uint32_t i = randu32();
		if (i == UINT_MAX)
			i = 0;
		v.push_back(i);
		if (countmap.find(i) == countmap.end())
			countmap[i] = 1;
		else
			{
			if (countmap[i] == 1)
				++DupeCount;
			countmap[i] += 1;
			}
		}

	Duper D(N);
	for (uint i = 0; i < N; ++i)
		D.Add(v[i]);
	asserta(DupeCount == D.m_DupeCount);

	Duper D2(1);
	for (uint i = 0; i < N; ++i)
		{
		D2.Init(N);
		D2.m_DupeCount = 0;
		for (uint j = 0; j < N; ++j)
			D2.Add(v[j]);
		asserta(DupeCount == D2.m_DupeCount);
		}
	return DupeCount;
	}

void cmd_duper()
	{
	const uint Iters = 10000;
	uint n = 0;
	for (uint Try = 0; Try < Iters; ++Try)
		{
		ProgressStep(Try, Iters, "Duping %u", n);
		n = Test();
		}
	}
