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

	m_Overflows = myalloc(list<uint8_t *> *, m_NrRdxs);
	zero_array(m_Overflows, m_NrRdxs);

	m_Data = myalloc(uint8_t, m_FixedBytes);
	zero_array(m_Data, m_FixedBytes);
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

	if (Size < m_FixedEntriesPerRdx)
		{
		// No overflow
		m_Data[Rdx*m_FixedEntriesPerRdx + Size];
		uint32_t Offset = Rdx*m_FixedEntriesPerRdx*m_ItemSize + Size*m_ItemSize;
		PutRaw(m_Data + Offset, SeqIdx, Diag);
		m_Sizes[Rdx] = Size + 1;
		m_Size += 1;
		}
	else if (Size%m_FixedEntriesPerRdx == 0)
		{
		// Overflow, needs new data
		list<uint8_t *> *Overflow = m_Overflows[Rdx];
		if (Overflow == 0)
			{
			// First overflow
			Overflow = new list<uint8_t *>;
			m_Overflows[Rdx] = Overflow;
			}
		uint8_t *NewData = myalloc(uint8_t, m_FixedEntriesPerRdx*m_ItemSize);
		Overflow->push_back(NewData);
		m_Overflows[Rdx] = Overflow;
		PutRaw(NewData, SeqIdx, Diag);
		m_Sizes[Rdx] = Size + 1;
		m_Size += 1;
		}
	else
		{
		// Previous overflow, fits into last data
		list<uint8_t *> *Overflow = m_Overflows[Rdx];
		assert(Overflow != 0);
		uint8_t *Data = Overflow->back();
		uint32_t Idx = Size%m_FixedEntriesPerRdx;
		PutRaw(Data + Idx*m_ItemSize, SeqIdx, Diag);
		m_Sizes[Rdx] = Size + 1;
		m_Size += 1;
		}
	}

void TwoHitDiag::LogRdx(uint Rdx) const
	{
	uint Size = m_Sizes[Rdx];
	Log("LogRdx(%08x) size=%u", Rdx, m_Sizes[Rdx]);
	const uint8_t *BasePtr = m_Data + Rdx*m_FixedEntriesPerRdx*m_ItemSize;
	uint n1 = min(Size, m_FixedEntriesPerRdx);
	for (uint i = 0; i < n1; ++i)
		{
		uint32_t Offset = i*m_ItemSize;
		asserta(Offset < m_FixedBytes);
		uint16_t Diag;
		uint32_t SeqIdx = GetRaw(BasePtr + Offset, Diag);
		Log(" %u:%u", SeqIdx, Diag);
		}
	Log("\n");

	if (Size <= m_FixedEntriesPerRdx)
		{
		asserta(m_Overflows[Rdx] == 0);
		return;
		}

	const list<uint8_t *> *Overflow = m_Overflows[Rdx];
	asserta(Overflow != 0);

	uint Sum = 0;
	uint BlockIdx = 0;
	uint Remainder = Size - m_FixedEntriesPerRdx;
	for (list<uint8_t *>::const_iterator iter = Overflow->begin();
		 iter != Overflow->end(); ++iter)
		{
		Log("  Overflow [%u] ", BlockIdx);
		++BlockIdx;
		// Test that Data is readable
		const uint8_t *Data = *iter;
		uint n = Remainder;
		if (n > m_FixedEntriesPerRdx)
			n = m_FixedEntriesPerRdx;
		Remainder -= n;
		for (uint i = 0; i < n; ++i)
			{
			const uint8_t *ptr = Data + i*m_ItemSize;
			uint16_t Diag;
			uint32_t SeqIdx = GetRaw(ptr, Diag);
			Log(" %u:%u,%u", SeqIdx, Diag);
			}
		Log("\n");
		}
	uint ExpectedBlockCount = (Size - 1)/m_FixedEntriesPerRdx;
	asserta(BlockIdx == ExpectedBlockCount);
	}

void TwoHitDiag::Reset()
	{
	for (uint i = 0; i < m_BusyCount; ++i)
		{
		uint32_t Rdx = m_BusyRdxs[i];
		assert(m_Sizes[Rdx] > 0);
		list<uint8_t *> *Overflow = m_Overflows[Rdx];
		if (Overflow != 0)
			{
			for (list<uint8_t *>::iterator iter = Overflow->begin();
				iter != Overflow->end(); ++iter)
				{
				uint8_t *Data = *iter;
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

	if (Size <= m_FixedEntriesPerRdx)
		{
		asserta(m_Overflows[Rdx] == 0);
		return;
		}

	const list<uint8_t *> *Overflow = m_Overflows[Rdx];
	asserta(Overflow != 0);

	uint Sum = 0;
	uint BlockIdx = 0;
	uint Remainder = Size - m_FixedEntriesPerRdx;
	for (list<uint8_t *>::const_iterator iter = Overflow->begin();
		 iter != Overflow->end(); ++iter)
		{
		++BlockIdx;
		// Test that Data is readable
		const uint8_t *Data = *iter;
		uint n = Remainder;
		if (n > m_FixedEntriesPerRdx)
			n = m_FixedEntriesPerRdx;
		Remainder -= n;
		for (uint i = 0; i < n; ++i)
			{
			const uint8_t *ptr = Data + i*m_ItemSize;
			uint16_t Diag;
			uint32_t SeqIdx = GetRaw(ptr, Diag);
			asserta(SeqIdx <= MaxSeqIdx);
			asserta(Diag <= MaxDiag);
			}
		}
	uint ExpectedBlockCount = (Size - 1)/m_FixedEntriesPerRdx;
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

static void TestPutGet()
	{
	TwoHitDiag T;
	for (uint Try = 0; Try < 100; ++Try)
		{
		uint32_t SeqIdx = randu32();
		uint16_t Diag = uint16_t(randu32());
		uint Offset = randu32()%(T.m_FixedBytes - 32);
		uint8_t *ptr = T.m_Data + Offset;
		T.PutRaw(ptr, SeqIdx, Diag);

		uint16_t Diag2;
		uint SeqIdx2 = T.GetRaw(ptr, Diag2);
		asserta(SeqIdx2 == SeqIdx);
		asserta(Diag2 == Diag);
		}
	ProgressLog("TestPutGet ok\n");
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

static void Test3(uint MaxSeqIdx, uint16_t MaxDiag, uint Tries)
	{
	TwoHitDiag T;

	for (uint i = 0; i < Tries; ++i)
		{
		uint SeqIdx = randu32()%(MaxSeqIdx + 1);
		uint16_t Diag = uint16_t(randu32()%MaxDiag + 1);
		uint n = 1 + randu32()%128;
		for (uint j = 0; j < n; ++j)
			T.Add(SeqIdx, Diag);
		}
	T.Validate(MaxSeqIdx, MaxDiag);
//	ProgressLog("Validate #1 ok\n");

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
	Test3(156789987, 555, 100);
	Test3(89987, 878, 1000);
	Test3(8339987, 999, 9999);
	}
