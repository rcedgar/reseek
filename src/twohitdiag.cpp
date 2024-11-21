#include "myutils.h"
#include "twohitdiag.h"

TwoHitDiag::~TwoHitDiag()
	{
// Prefer to leak memory for faster exit
	}

TwoHitDiag::TwoHitDiag()
	{
	m_BusyRdxs = myalloc(uint32_t, m_NrRdxs);
	zero_array(m_BusyRdxs, m_NrRdxs);

	m_Sizes = myalloc(uint32_t, m_NrRdxs);
	zero_array(m_Sizes, m_NrRdxs);

	m_Overflows = myalloc(list<uint32_t *> *, m_NrRdxs);
	zero_array(m_Overflows, m_NrRdxs);

	m_Data = myalloc(uint32_t, m_TotalFixedItems);
	zero_array(m_Data, m_TotalFixedItems);
	}

void TwoHitDiag::Add(uint32_t SeqIdx, uint16_t Diag)
	{
	uint32_t Rdx = GetRdx(SeqIdx, Diag);
	assert(Rdx < m_NrRdxs);
	uint32_t Size = m_Sizes[Rdx];

	if (Size == 0)
		{
		assert(m_BusyCount < m_NrRdxs);
		m_BusyRdxs[m_BusyCount++] = Rdx;
		}

	if (Size < m_FixedItemsPerRdx)
		{
		// No overflow
		m_Data[Rdx*m_FixedItemsPerRdx + Size];
		uint32_t Offset = Rdx*m_FixedItemsPerRdx + Size;
		PutRaw(m_Data + Offset, SeqIdx, Diag);
		m_Sizes[Rdx] = Size + 1;
		m_Size += 1;
		}
	else if (Size%m_FixedItemsPerRdx == 0)
		{
		// Overflow, needs new data
		list<uint32_t *> *Overflow = m_Overflows[Rdx];
		if (Overflow == 0)
			{
			// First overflow
			Overflow = new list<uint32_t *>;
			m_Overflows[Rdx] = Overflow;
			}
		uint32_t *NewData = myalloc(uint32_t, m_FixedItemsPerRdx);
		Overflow->push_back(NewData);
		m_Overflows[Rdx] = Overflow;
		PutRaw(NewData, SeqIdx, Diag);
		m_Sizes[Rdx] = Size + 1;
		m_Size += 1;
		}
	else
		{
		// Previous overflow, fits into last data
		list<uint32_t *> *Overflow = m_Overflows[Rdx];
		assert(Overflow != 0);
		uint32_t *Data = Overflow->back();
		uint32_t Idx = Size%m_FixedItemsPerRdx;
		PutRaw(Data + Idx, SeqIdx, Diag);
		m_Sizes[Rdx] = Size + 1;
		m_Size += 1;
		}
	}

void TwoHitDiag::LogRdx(uint Rdx) const
	{
	uint Size = m_Sizes[Rdx];
	Log("LogRdx(%08x) size=%u", Rdx, m_Sizes[Rdx]);
	const uint32_t *BasePtr = m_Data + Rdx*m_FixedItemsPerRdx;
	uint n1 = min(Size, m_FixedItemsPerRdx);
	for (uint i = 0; i < n1; ++i)
		{
		uint16_t Diag;
		uint32_t SeqIdx = GetAllBits(Rdx, BasePtr + i, Diag);
		Log(" %u:%u", SeqIdx, Diag);
		}
	Log("\n");

	if (Size <= m_FixedItemsPerRdx)
		{
		asserta(m_Overflows[Rdx] == 0);
		return;
		}

	const list<uint32_t *> *Overflow = m_Overflows[Rdx];
	asserta(Overflow != 0);

	uint Sum = 0;
	uint BlockIdx = 0;
	uint Remainder = Size - m_FixedItemsPerRdx;
	for (list<uint32_t *>::const_iterator iter = Overflow->begin();
		 iter != Overflow->end(); ++iter)
		{
		Log("  Overflow [%u] ", BlockIdx);
		++BlockIdx;
		// Test that Data is readable
		const uint32_t *Data = *iter;
		uint n = Remainder;
		if (n > m_FixedItemsPerRdx)
			n = m_FixedItemsPerRdx;
		Remainder -= n;
		for (uint i = 0; i < n; ++i)
			{
			const uint32_t *ptr = Data + i;
			uint16_t Diag;
			uint32_t SeqIdx = GetAllBits(Rdx, ptr, Diag);
			Log(" %u:%u", SeqIdx, Diag);
			}
		Log("\n");
		}
	uint ExpectedBlockCount = (Size - 1)/m_FixedItemsPerRdx;
	asserta(BlockIdx == ExpectedBlockCount);
	}

void TwoHitDiag::Reset()
	{
	for (uint i = 0; i < m_BusyCount; ++i)
		{
		uint32_t Rdx = m_BusyRdxs[i];
		assert(m_Sizes[Rdx] > 0);
		list<uint32_t *> *Overflow = m_Overflows[Rdx];
		if (Overflow != 0)
			{
			for (list<uint32_t *>::iterator iter = Overflow->begin();
				iter != Overflow->end(); ++iter)
				{
				uint32_t *Data = *iter;
				myfree(Data);
				}
			delete Overflow;
			m_Overflows[Rdx] = 0;
			}
		m_Sizes[Rdx] = 0;
		m_BusyRdxs[i] = 0;
		}
	m_BusyRdxs = 0;
	m_BusyCount = 0;
	m_Size = 0;
	}

void TwoHitDiag::ValidateEmptyRdx(uint Rdx) const
	{
	asserta(m_Sizes[Rdx] == 0);
	asserta(m_Overflows[Rdx] == 0);
	}

void TwoHitDiag::ValidateRdx(uint Rdx, uint MaxSeqIdx, uint MaxDiag) const
	{
	uint Size = m_Sizes[Rdx];
	if (Size == 0)
		{
		ValidateEmptyRdx(Rdx);
		return;
		}

	if (Size <= m_FixedItemsPerRdx)
		{
		asserta(m_Overflows[Rdx] == 0);
		return;
		}

	const list<uint32_t *> *Overflow = m_Overflows[Rdx];
	asserta(Overflow != 0);

	uint Sum = 0;
	uint BlockIdx = 0;
	uint Remainder = Size - m_FixedItemsPerRdx;
	for (list<uint32_t *>::const_iterator iter = Overflow->begin();
		 iter != Overflow->end(); ++iter)
		{
		++BlockIdx;
		// Test that Data is readable
		const uint32_t *Data = *iter;
		uint n = Remainder;
		if (n > m_FixedItemsPerRdx)
			n = m_FixedItemsPerRdx;
		Remainder -= n;
		for (uint i = 0; i < n; ++i)
			{
			const uint32_t *ptr = Data + i;
			uint16_t Diag;
			uint32_t SeqIdx = GetAllBits(Rdx, ptr, Diag);
			asserta(SeqIdx <= MaxSeqIdx);
			asserta(Diag <= MaxDiag);
			}
		}
	uint ExpectedBlockCount = (Size - 1)/m_FixedItemsPerRdx;
	asserta(BlockIdx == ExpectedBlockCount);
	}

void TwoHitDiag::ValidateEmpty() const
	{
	for (uint Rdx = 0; Rdx < m_NrRdxs; ++Rdx)
		ValidateEmptyRdx(Rdx);
	asserta(m_BusyCount == 0);
	}

void TwoHitDiag::Validate(uint MaxSeqIndex, uint MaxDiag) const
	{
	for (uint Rdx = 0; Rdx < m_NrRdxs; ++Rdx)
		ValidateRdx(Rdx, MaxSeqIndex, MaxDiag);

	uint BusyCount = 0;
	uint TotalSize = 0;
	for (uint Rdx = 0; Rdx < m_NrRdxs; ++Rdx)
		{
		uint Size = m_Sizes[Rdx];
		if (Size > 0)
			{
			TotalSize += Size;
			++BusyCount;
			}
		}
	asserta(m_Size == TotalSize);
	asserta(m_BusyCount == BusyCount);
	}

void TwoHitDiag::LogStats() const
	{
	uint BusyCount = 0;
	uint OverflowCount = 0;
	for (uint Rdx = 0; Rdx < m_NrRdxs; ++Rdx)
		{
		if (m_Sizes[Rdx] > 0)
			{
			++BusyCount;
			if (m_Overflows[Rdx] != 0)
				++OverflowCount;
			}
		}
	Log("%u rdx, %u busy, %u overflow\n",
		m_NrRdxs, BusyCount, OverflowCount);
	asserta(m_BusyCount == BusyCount);
	}

void TwoHitDiag::TestPutGet(uint32_t SeqIdx, uint16_t Diag) const
	{
	asserta(Diag <= m_Mask14);

	Log("TracePutGet()\n");
	Log("%8x  SeqIdx\n", SeqIdx);
	Log("%8x  Diag\n", Diag);

	uint32_t SeqIdxLoBits = SeqIdx%m_SeqMod;
	uint16_t DiagLoBits = Diag%m_DiagMod;
	uint32_t Rdx = (SeqIdxLoBits | (DiagLoBits << m_DiagShiftRdx));

	Log("%8x  SeqIdxLoBits\n", SeqIdxLoBits);
	Log("%8x  DiagLoBits\n", DiagLoBits);
	Log("%8x  Rdx\n", Rdx);

	uint32_t DiagLoBits2 = (Rdx >> m_DiagShiftRdx);
	uint32_t SeqIdxLoBits2 = (Rdx & m_SeqIdxMaskRdx);

	Log("%8x  SeqIdxLoBits2\n", SeqIdxLoBits);
	Log("%8x  DiagLoBits2\n", DiagLoBits);

	asserta(DiagLoBits2 == DiagLoBits);
	asserta(SeqIdxLoBits2 == SeqIdxLoBits);

	uint32_t SeqIdxHiBits = SeqIdx/m_SeqMod;
	uint32_t DiagHiBits = Diag/m_DiagMod;

	asserta(SeqIdx == ((SeqIdxHiBits*m_SeqMod) + SeqIdxLoBits));
	asserta(SeqIdx == ((SeqIdxHiBits*m_SeqMod) | SeqIdxLoBits));
	
	asserta(Diag == ((DiagHiBits*m_DiagMod) + DiagLoBits));
	asserta(Diag == ((DiagHiBits*m_DiagMod) | DiagLoBits));

	uint32_t Item = (SeqIdxHiBits | (DiagHiBits << m_DiagShiftItem));

	Log("%8x  SeqIdxHiBits\n", SeqIdxHiBits);
	Log("%8x  DiagHiBits\n", DiagHiBits);
	Log("%8x  Item\n", Item);

	uint32_t SeqIdxHiBits2 = (Item & m_SeqIdxMaskItem);
	uint16_t DiagHiBits2 = uint16_t(Item >> m_DiagShiftItem);

	Log("%8x  SeqIdxHiBits2\n", SeqIdxHiBits2);
	Log("%8x  DiagHiBits2\n", DiagHiBits2);

	asserta(DiagHiBits2 == DiagHiBits);
	asserta(SeqIdxHiBits2 == SeqIdxHiBits);
	}

void TwoHitDiag::AppendAll(uint Rdx, set<pair<uint32_t, uint16_t> > &SeqIdxDiagPairs) const
	{
	uint Size = m_Sizes[Rdx];
	const uint32_t *BasePtr = m_Data + Rdx*m_FixedItemsPerRdx;
	uint n1 = min(Size, m_FixedItemsPerRdx);
	for (uint i = 0; i < n1; ++i)
		{
		uint16_t Diag;
		uint32_t SeqIdx = GetAllBits(Rdx, BasePtr + i, Diag);
		SeqIdxDiagPairs.insert(pair<uint32_t, uint16_t>(SeqIdx, Diag));
		}

	if (Size <= m_FixedItemsPerRdx)
		{
		asserta(m_Overflows[Rdx] == 0);
		return;
		}

	const list<uint32_t *> *Overflow = m_Overflows[Rdx];
	asserta(Overflow != 0);

	uint Sum = 0;
	uint Remainder = Size - m_FixedItemsPerRdx;
	for (list<uint32_t *>::const_iterator iter = Overflow->begin();
		 iter != Overflow->end(); ++iter)
		{
		const uint32_t *Data = *iter;
		uint n = Remainder;
		if (n > m_FixedItemsPerRdx)
			n = m_FixedItemsPerRdx;
		Remainder -= n;
		for (uint i = 0; i < n; ++i)
			{
			const uint32_t *ptr = Data + i;
			uint16_t Diag;
			uint32_t SeqIdx = GetAllBits(Rdx, ptr, Diag);
			SeqIdxDiagPairs.insert(pair<uint32_t, uint16_t>(SeqIdx, Diag));
			}
		}
	}

void TwoHitDiag::GetAll(set<pair<uint32_t, uint16_t> > &SeqIdxDiagPairs) const
	{
	SeqIdxDiagPairs.clear();
	for (uint Rdx = 0; Rdx < m_NrRdxs; ++Rdx)
		AppendAll(Rdx, SeqIdxDiagPairs);
	}

static void TestPutGet()
	{
	TwoHitDiag T;
	for (uint Try = 0; Try < 100; ++Try)
		{
		uint32_t SeqIdx = randu32();
		uint16_t Diag = uint16_t(randu32()%T.m_DiagHi);
		uint Offset = randu32()%(T.m_TotalFixedItems - 32);
		uint32_t *ptr = T.m_Data + Offset;
		uint Rdx = T.GetRdx(SeqIdx, Diag);
		T.PutRaw(ptr, SeqIdx, Diag);

		uint16_t Diag2;
		uint SeqIdx2 = T.GetAllBits(Rdx, ptr, Diag2);
		asserta(SeqIdx2 == SeqIdx);
		asserta(Diag2 == Diag);
		}
	ProgressLog("TestPutGet ok\n");
	}

static void TestPutGet1()
	{
	TwoHitDiag T;

	uint32_t SeqIdx = 0xb519e073;
	uint16_t Diag = 0xdcd;
	T.TestPutGet(SeqIdx, Diag);
	return;

	uint Offset = 0x0003265e;
	uint32_t *ptr = T.m_Data + Offset;
	uint Rdx = T.GetRdx(SeqIdx, Diag);
	T.PutRaw(ptr, SeqIdx, Diag);

	uint16_t Diag2;
	uint SeqIdx2 = T.GetAllBits(Rdx, ptr, Diag2);

	asserta(SeqIdx2 == SeqIdx);
	asserta(Diag2 == Diag);
	}

static void SimpleTests()
	{
	TestPutGet();

	TwoHitDiag T;
	T.ValidateEmpty();
	ProgressLog("ValidateEmpty #1 ok\n");

	T.Validate(0, 0);
	ProgressLog("Validate #1 ok\n");

	T.Add(1, 3);
	ProgressLog("Validate #2 ok\n");

	T.Reset();
	T.ValidateEmpty();
	ProgressLog("ValidateEmpty #2 ok\n");
	}

static void Test2()
	{
	TwoHitDiag T;
	uint32_t SeqIdx = 0;
	uint16_t Diag = 123;
	for (uint Pos = 0; Pos < 23; ++Pos)
		T.Add(SeqIdx, Diag);

	uint Rdx = T.GetRdx(SeqIdx, Diag);
	T.LogRdx(Rdx);
	T.ValidateRdx(Rdx, 0, 999);
	T.Validate(0, 9999);
	T.Reset();
	T.ValidateEmpty();
	}

static void Test3(uint MaxSeqIdx, uint16_t MaxDiag, uint Tries)
	{
	TwoHitDiag T;

	set<pair<uint32_t, uint16_t> > Pairs;
	for (uint i = 0; i < Tries; ++i)
		{
		uint SeqIdx = randu32()%(MaxSeqIdx + 1);
		uint16_t Diag = uint16_t(randu32()%MaxDiag + 1);
		uint n = 1 + randu32()%128;
		for (uint j = 0; j < n; ++j)
			{
			T.Add(SeqIdx, Diag);
			Pairs.insert(pair<uint32_t, uint16_t>(SeqIdx, Diag));
			}
		}
	T.Validate(MaxSeqIdx, MaxDiag);
//	ProgressLog("Validate #1 ok\n");

	set<pair<uint32_t, uint16_t> > Pairs2;
	T.GetAll(Pairs2);

	for (set<pair<uint32_t, uint16_t> >::const_iterator iter = Pairs.begin();
		 iter != Pairs.end(); ++iter)
		{
		const pair<uint32_t, uint16_t> &Pair = *iter;
		set<pair<uint32_t, uint16_t> >::const_iterator iter2 = Pairs2.find(Pair);
		asserta(iter2 != Pairs2.end());
		}
	asserta(SIZE(Pairs2) == SIZE(Pairs));

	T.Reset();
	T.Validate(0, 0);
//	ProgressLog("Validate #2 ok\n");

	T.ValidateEmpty();
//	ProgressLog("ValidateEmpty #3 ok\n");

	T.LogStats();
	ProgressLog("Test3(%u, %u) Ok\n", MaxSeqIdx, Tries);
	}

void cmd_twohit()
	{
	TestPutGet1();
	TestPutGet();
	SimpleTests();
	Test3(156789987, 555, 100);
	Test3(89987, 878, 1000);
	Test3(8339987, 999, 9999);
	}
