#include "myutils.h"
#include "seqdb.h"
#include "twohitdiag.h"
#include "duper.h"
#include "mudex.h"
#include "mermx.h"
#include "alpha.h"
#include "diaghsp.h"
#include "diag.h"

#define LOGDIAGALNS	0

void StrToMuLetters(const string &StrSeq, byte *Letters)
	{
	const uint L = SIZE(StrSeq);
	for (uint i = 0; i < L; ++i)
		{
		byte Letter = g_CharToLetterMu[StrSeq[i]];
		assert(Letter < 36);
		Letters[i] = Letter;
		}
	}

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

	m_MaxSeqIdxHiBits = myalloc(uint32_t, m_NrRdxs);
	zero_array(m_MaxSeqIdxHiBits, m_NrRdxs);

	m_Overflows = myalloc(list<uint32_t *> *, m_NrRdxs);
	zero_array(m_Overflows, m_NrRdxs);

	m_Data = myalloc(uint32_t, m_TotalFixedItems);
	zero_array(m_Data, m_TotalFixedItems);
	}

void TwoHitDiag::Add(uint32_t SeqIdx, uint16_t Diag)
	{
	uint32_t Rdx = GetRdx(SeqIdx, Diag);
	assert(Rdx < m_NrRdxs);

	uint32_t SeqIdxHiBits = SeqIdx/m_SeqMod;
	m_MaxSeqIdxHiBits[Rdx] = max(m_MaxSeqIdxHiBits[Rdx], SeqIdxHiBits);

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
		m_MaxSeqIdxHiBits[Rdx] = 0;
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
	uint32_t MaxSeqIdxHiBits = 0;
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

	uint n1 = min(Size, m_FixedItemsPerRdx);
	const uint32_t *BasePtr = m_Data + Rdx*m_FixedItemsPerRdx;
	for (uint i = 0; i < n1; ++i)
		{
		uint16_t Diag;
		uint32_t SeqIdx = GetAllBits(Rdx, BasePtr + i, Diag);
		uint32_t SeqIdxHiBits = SeqIdx/m_SeqMod;
		MaxSeqIdxHiBits = max(SeqIdxHiBits, MaxSeqIdxHiBits);
		}

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

			uint32_t SeqIdxHiBits = SeqIdx/m_SeqMod;
			MaxSeqIdxHiBits = max(SeqIdxHiBits, MaxSeqIdxHiBits);
			}
		}
	uint ExpectedBlockCount = (Size - 1)/m_FixedItemsPerRdx;
	asserta(BlockIdx == ExpectedBlockCount);
	asserta(m_MaxSeqIdxHiBits[Rdx] == MaxSeqIdxHiBits);
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

void TwoHitDiag::AppendAll(uint Rdx, vector<pair<uint32_t, uint16_t> > &SeqIdxDiagPairs) const
	{
	uint Size = m_Sizes[Rdx];
	const uint32_t *BasePtr = m_Data + Rdx*m_FixedItemsPerRdx;
	uint n1 = min(Size, m_FixedItemsPerRdx);
	for (uint i = 0; i < n1; ++i)
		{
		uint16_t Diag;
		uint32_t SeqIdx = GetAllBits(Rdx, BasePtr + i, Diag);
		SeqIdxDiagPairs.push_back(pair<uint32_t, uint16_t>(SeqIdx, Diag));
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
			brk(SeqIdx == 1 && Diag == 432);
			SeqIdxDiagPairs.push_back(pair<uint32_t, uint16_t>(SeqIdx, Diag));
			}
		}
	}

void TwoHitDiag::GetAll(vector<pair<uint32_t, uint16_t> > &SeqIdxDiagPairs) const
	{
	SeqIdxDiagPairs.clear();
	for (uint Rdx = 0; Rdx < m_NrRdxs; ++Rdx)
		AppendAll(Rdx, SeqIdxDiagPairs);
	}

void TwoHitDiag::AddItems(Duper &D, uint Rdx) const
	{
	uint Size = m_Sizes[Rdx];
	const uint32_t *Data = m_Data + Rdx*m_FixedItemsPerRdx;
	uint n1 = min(Size, m_FixedItemsPerRdx);
	for (uint i = 0; i < n1; ++i)
		D.Add(Data[i]);
	if (Size <= m_FixedItemsPerRdx)
		{
		assert(m_Overflows[Rdx] == 0);
		return;
		}

	const list<uint32_t *> *Overflow = m_Overflows[Rdx];
	assert(Overflow != 0);

	uint Sum = 0;
	uint BlockIdx = 0;
	uint Remainder = Size - m_FixedItemsPerRdx;
	for (list<uint32_t *>::const_iterator iter = Overflow->begin();
		iter != Overflow->end(); ++iter)
		{
		++BlockIdx;
		Data = *iter;
		uint n = Remainder;
		if (n > m_FixedItemsPerRdx)
			n = m_FixedItemsPerRdx;
		Remainder -= n;
		for (uint i = 0; i < n; ++i)
			D.Add(Data[i]);
		}
	uint ExpectedBlockCount = (Size - 1)/m_FixedItemsPerRdx;
	asserta(BlockIdx == ExpectedBlockCount);
	}

void TwoHitDiag::SetDupesRdx(uint Rdx)
	{
	uint Size = m_Sizes[Rdx];
	if (Size < 2)
		return;
	Duper &D = *new Duper(Size);
	AddItems(D, Rdx);
	for (uint j = 0; j < D.m_DupeCount; ++j)
		{
		uint32_t Item = D.m_Dupes[j];

		uint16_t Diag;
		uint32_t SeqIdx = CvtItem(Rdx, Item, Diag);

		m_DupeSeqIdxs[m_DupeCount] = SeqIdx;
		m_DupeDiags[m_DupeCount] = Diag;
		++m_DupeCount;
		}
	}

void TwoHitDiag::SetDupes()
	{
	m_DupeCount = 0;
	if (m_Size < 2)
		return;
	m_DupeSeqIdxs = myalloc(uint32_t, m_Size);
	m_DupeDiags = myalloc(uint16_t, m_Size);
	for (uint i = 0; i < m_BusyCount; ++i)
		SetDupesRdx(m_BusyRdxs[i]);
	}

void TwoHitDiag::ClearDupes()
	{
	if (m_DupeSeqIdxs != 0)
		{
		assert(m_DupeDiags != 0);
		myfree(m_DupeSeqIdxs);
		myfree(m_DupeDiags);
		m_DupeSeqIdxs = 0;
		m_DupeDiags = 0;
		}
	else
		assert(m_DupeDiags == 0);
	}

void TwoHitDiag::CheckDupes()
	{
	vector<pair<uint32_t, uint16_t> > Pairs;
	GetAll(Pairs);
	map<pair<uint32_t, uint16_t>, uint32_t> PairToCount;
	uint DupeCount = 0;
	for (uint i = 0; i < SIZE(Pairs); ++i)
		{
		const pair<uint32_t, uint16_t> &Pair = Pairs[i];
		if (PairToCount.find(Pair) == PairToCount.end())
			PairToCount[Pair] = 1;
		else
			{
			uint n = PairToCount[Pair];
			if (n == 1)
				++DupeCount;
			PairToCount[Pair] = n + 1;
			}
		}
	for (uint i = 0; i < m_DupeCount; ++i)
		{
		uint32_t SeqIdx = m_DupeSeqIdxs[i];
		uint16_t Diag = m_DupeDiags[i];
		pair<uint32_t, uint16_t> Pair(SeqIdx, Diag);
		map<pair<uint32_t, uint16_t>, uint32_t>::const_iterator iter =
			PairToCount.find(Pair);
		asserta(iter != PairToCount.end());
		asserta(iter->second > 1);
		}
	asserta(DupeCount == m_DupeCount);
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

static void VecToSet(const vector<pair<uint32_t, uint16_t> > &PairsVec,
					 set<pair<uint32_t, uint16_t>> &Pairs)
	{
	Pairs.clear();
	const uint N = SIZE(PairsVec);
	for (uint i = 0; i < N; ++i)
		Pairs.insert(PairsVec[i]);
	}

static void Test3(uint MaxSeqIdx, uint16_t MaxDiag, uint Tries)
	{
	TwoHitDiag T;

	vector<pair<uint32_t, uint16_t> > PairsVec;
	for (uint i = 0; i < Tries; ++i)
		{
		uint SeqIdx = randu32()%(MaxSeqIdx + 1);
		uint16_t Diag = uint16_t(randu32()%MaxDiag + 1);
		uint n = 1 + randu32()%128;
		for (uint j = 0; j < n; ++j)
			{
			T.Add(SeqIdx, Diag);
			PairsVec.push_back(pair<uint32_t, uint16_t>(SeqIdx, Diag));
			}
		}
	T.Validate(MaxSeqIdx, MaxDiag);
//	ProgressLog("Validate #1 ok\n");

	set<pair<uint32_t, uint16_t> > Pairs;
	VecToSet(PairsVec, Pairs);

	vector<pair<uint32_t, uint16_t> > Pairs2Vec;
	T.GetAll(Pairs2Vec);

	set<pair<uint32_t, uint16_t> > Pairs2;
	VecToSet(Pairs2Vec, Pairs2);

	for (set<pair<uint32_t, uint16_t> >::const_iterator iter = Pairs.begin();
		 iter != Pairs.end(); ++iter)
		{
		const pair<uint32_t, uint16_t> &Pair = *iter;
		set<pair<uint32_t, uint16_t> >::const_iterator iter2 = Pairs2.find(Pair);
		asserta(iter2 != Pairs2.end());
		}
	asserta(SIZE(Pairs2) == SIZE(Pairs));

	T.SetDupes();
	uint DupeCount = T.m_DupeCount;
	T.CheckDupes();
//	T.ClearDupes();

	T.Reset();
	T.Validate(0, 0);
//	ProgressLog("Validate #2 ok\n");

	T.ValidateEmpty();
//	ProgressLog("ValidateEmpty #3 ok\n");

	T.LogStats();
	ProgressLog("Test3(%u, %u, %u) Ok, %u dupes\n",
				MaxSeqIdx, MaxDiag, Tries, DupeCount);
	}

static void SavedTests()
	{
	TestPutGet1();
	TestPutGet();
	SimpleTests();

	for (uint Try = 0; Try < 100; ++Try)
		{
		uint MaxSeqIndex = 12345 + randu32()%8339987;
		uint16_t MaxDiag = (123 + randu32()%63001) & TwoHitDiag::m_Mask14;
		if (MaxDiag < 10)
			continue;
		uint Tries = 23 + randu32()%222;
		Test3(MaxSeqIndex, MaxDiag, Tries);
		}
	}

void LogDiagAln(const byte *MuSeqQ, uint LQ, const char *LabelQ,
				const byte *MuSeqT, uint LT, const char *LabelT,
				int Diag, int Lo, int Len)
	{
	extern const short * const *Mu_S_ij_short;
	diag dg(LQ, LT);
	int mini = dg.getmini(Diag) + Lo;
	int minj = dg.getminj(Diag) + Lo;
	Log("\n");
	Log("%5d  ", mini+1);
	for (int diagpos = 0; diagpos < Len; ++diagpos)
		{
		int posq = mini+diagpos;
		asserta(posq < int(LQ));
		char c = g_LetterToCharMu[MuSeqQ[posq]];
		Log("%c", c);
		}
	Log("  >%s (%u)\n", LabelQ, LQ);

	int d00 = dg.getd(0, 0);

	int Sum = 0;
	Log("       ");
	for (int diagpos = 0; diagpos < Len; ++diagpos)
		{
		int posq = mini+diagpos;
		int post = minj+diagpos;
		asserta(posq < int(LQ));
		asserta(post < int(LT));
		byte iq = MuSeqQ[posq];
		byte it = MuSeqT[post];
		int Score = Mu_S_ij_short[iq][it];
		Sum += Score;
		if (Score >= 5)
			Log("|");
		else if (Score > 0)
			Log(".");
		else
			Log(" ");
		}
	Log("   score %d (d=%d/%+d)\n", Sum, Diag, Diag-d00);

	Log("%5d  ", minj+1);
	for (int diagpos = 0; diagpos < Len; ++diagpos)
		{
		int post = minj+diagpos;
		asserta(post < int(LT));
		char c = g_LetterToCharMu[MuSeqT[post]];
		Log("%c", c);
		}
	Log("  >%s (%u)\n", LabelT, LT);
	}

void cmd_twohit()
	{
	const int HSMinScore = 35;
	const int MinDiagScore = 130;
	bool Self = false;
	FILE *fOut = CreateStdioFile(opt_output);

	SeqDB Input;
	Input.FromFasta(g_Arg1);
	Input.SetLabelToIndex();
	const uint SeqCount = Input.GetSeqCount();

	byte **MuSeqs = myalloc(byte *, SeqCount);
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		uint L = Input.GetSeqLength(SeqIdx);
		byte *MuSeq = myalloc(byte, L);
		StrToMuLetters(Input.GetSeq(SeqIdx), MuSeq);
		MuSeqs[SeqIdx] = MuSeq;
		}

	MuDex MD;
	MD.FromSeqDB(Input);
	const uint k = MD.m_k;
	const uint DictSize = MD.m_DictSize;

	MerMx MM;
	extern const short * const *Mu_S_ij_short;
	MM.Init(Mu_S_ij_short, 5, 36, 2);
	asserta(MM.m_AS_pow[5] == DictSize);

	DiagHSP DH;

// Buffer for high-scoring k-mers
	uint *HSKmers = myalloc(uint, DictSize);

	int *SeqIdxTToBestDiagScore = myalloc(int, SeqCount);
	int *SeqIdxTToBestDiag = myalloc(int, SeqCount);
	int *SeqIdxTToBestDiagLo = myalloc(int, SeqCount);
	int *SeqIdxTToBestDiagLen = myalloc(int, SeqCount);
	uint *SeqIdxTs = myalloc(uint, SeqCount);
	uint HitCount = 0;

	zero_array(SeqIdxTToBestDiagScore, SeqCount);
	zero_array(SeqIdxTToBestDiag, SeqCount);
	zero_array(SeqIdxTToBestDiagLo, SeqCount);
	zero_array(SeqIdxTToBestDiagLen, SeqCount);
	zero_array(SeqIdxTs, SeqCount);

	LeakCheck("Before_loop");

	TwoHitDiag TH;
	for (uint SeqIdxQ = 0; SeqIdxQ < SeqCount; ++SeqIdxQ)
		{
		ProgressStep(SeqIdxQ, SeqCount, "Searching");
		if (SeqIdxQ == 5)
			LeakCheck("SeqIdxQ == 5");
		else if (SeqIdxQ == 10)
			{
			LeakCheck("SeqIdxQ == 10");
			Die("TODO");
			}
		const char *SeqQ = Input.GetSeq(SeqIdxQ).c_str();
		const char *LabelQ = Input.GetLabel(SeqIdxQ).c_str();
		uint LQ32 = Input.GetSeqLength(SeqIdxQ);
		asserta(LQ32 < UINT16_MAX);
		uint16_t LQ = uint16_t(LQ32);
#if 0
		{
		const char *LabelQ = Input.GetLabel(SeqIdxQ).c_str();
		Log("Q>%s\n", LabelQ);
		}
#endif

	// Initialize k-mer scan of Query
		uint Kmer = 0;
		for (uint SeqPosQ = 0; SeqPosQ < k-1; ++SeqPosQ)
			{
			byte Letter = g_CharToLetterMu[SeqQ[SeqPosQ]];
			assert(Letter < 36);
			Kmer = Kmer*36 + Letter;
			}

	// k-mer scan of Query
		for (uint SeqPosQ = k-1; SeqPosQ < LQ; ++SeqPosQ)
			{
			byte Letter = g_CharToLetterMu[SeqQ[SeqPosQ]];
			assert(Letter < 36);
			Kmer = Kmer*36 + Letter;
			Kmer %= DictSize;
		// Construct high-scoring neighborhood of current k-mer (Kmer)
			const uint HSKmerCount =
				MM.GetHighScoring5mers(Kmer, HSMinScore, HSKmers);
#if 0
			{
			if (HSKmerCount > 0)
				{
				string Tmp;
				Log(" %5u  %s  HSKmerCount=%u\n",
					SeqPosQ-4, MM.KmerToStr(Kmer, k, Tmp), HSKmerCount);
				}
			}
#endif
			for (uint HSKmerIdx = 0; HSKmerIdx < HSKmerCount; ++HSKmerIdx)
				{
				uint HSKmer = HSKmers[HSKmerIdx];
			// Look up HSKmer in MuDex (Mu k-mer index of db)
				uint RowSize = MD.GetRowSize(HSKmer);
				uint DataOffset = MD.GetRowStart(HSKmer);
#if 0
				{
				if (RowSize > 1)
					{
					string Tmp;
					Log("   %s row=%u", MM.KmerToStr(HSKmer, k, Tmp), RowSize);
					}
				}
#endif
				for (uint ColIdx = 0; ColIdx < RowSize; ++ColIdx)
					{
					uint32_t SeqIdxT;
					uint16_t SeqPosT;
					MD.Get(DataOffset++, SeqIdxT, SeqPosT);
					if (SeqIdxT == SeqIdxQ)
						continue;
					if (Self && SeqIdxT < SeqIdxQ)
						continue;
					uint LT32 = Input.GetSeqLength(SeqIdxT);
					asserta(LT32 < UINT16_MAX);
					uint16_t LT = uint16_t(LT32);
					diag dg(LQ, LT);
					uint16_t Diag = dg.getd(SeqPosQ, SeqPosT);
					TH.Add(SeqIdxT, Diag);
#if 0
					{
					Log(" %u:%u", SeqIdxT, Diag);
					}
#endif
					}
#if 0
				{
				if (RowSize > 1)
					Log("\n");
				}
#endif
				}
			}
		TH.ClearDupes();
		TH.SetDupes();
#if DEBUG
		TH.Validate(SeqCount, INT16_MAX);
#endif
		vector<pair<uint32_t, uint16_t> > SeqIdxDiagPairs;
		uint DupeCount = TH.m_DupeCount;
		const byte *MuSeqQ = MuSeqs[SeqIdxQ];
		DH.SetQ(MuSeqQ, LQ);
#if 0
		Log("%u dupes\n", DupeCount);
		Log(" Idx   Diag  Score     Lo    Len\n");
		//  12345  12345  12345  12345  12345
#endif
		assert(HitCount == 0);
		for (uint i = 0; i < DupeCount; ++i)
			{
			uint32_t SeqIdxT = TH.m_DupeSeqIdxs[i];
			uint16_t Diag = TH.m_DupeDiags[i];
			const byte *MuSeqT = MuSeqs[SeqIdxT];
			const uint LT = Input.GetSeqLength(SeqIdxT);
			DH.SetT(MuSeqT, LT);
			int Lo, Len;
			int DiagScore = DH.Search(Diag, Lo, Len);
#if 0
			Log("%5u  %5u  %5d  %5d  %5d  >%s (%u)\n",
				SeqIdxT, Diag, DiagScore, Lo, Len,
				Input.GetLabel(SeqIdxT).c_str(), LT);
#endif
			if (DiagScore >= MinDiagScore)
				{
				int BestDiagScoreT = SeqIdxTToBestDiagScore[SeqIdxT];
				if (BestDiagScoreT == 0)
					{
					SeqIdxTs[HitCount++] = SeqIdxT;
					SeqIdxTToBestDiagScore[SeqIdxT] = DiagScore;
					SeqIdxTToBestDiag[SeqIdxT] = Diag;
					SeqIdxTToBestDiagLo[SeqIdxT] = Lo;
					SeqIdxTToBestDiagLen[SeqIdxT] = Len;
					}
				else if (DiagScore > SeqIdxTToBestDiagScore[SeqIdxT])
					{
					SeqIdxTToBestDiagScore[SeqIdxT] = DiagScore;
					SeqIdxTToBestDiag[SeqIdxT] = Diag;
					SeqIdxTToBestDiagLo[SeqIdxT] = Lo;
					SeqIdxTToBestDiagLen[SeqIdxT] = Len;
					}
				}
			}
#if 1
		{
		if (HitCount > 0)
			{
			for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
				{
				uint SeqIdxT = SeqIdxTs[HitIdx];
				int BestDiagScore = SeqIdxTToBestDiagScore[SeqIdxT];
				int BestDiag = SeqIdxTToBestDiag[SeqIdxT];
				int BestDiagLo = SeqIdxTToBestDiagLo[SeqIdxT];
				int BestDiagLen = SeqIdxTToBestDiagLen[SeqIdxT];
				const char *LabelT = Input.GetLabel(SeqIdxT).c_str();
				const byte *MuSeqT = MuSeqs[SeqIdxT];
				uint LT = Input.GetSeqLength(SeqIdxT);
				if (fOut != 0)
					fprintf(fOut, "%s\t%s\n", LabelQ, LabelT);
#if LOGDIAGALNS
				Log("%6d  %6d  >%s\n", BestDiagScore, BestDiag, LabelT);
				LogDiagAln(MuSeqQ, LQ, LabelQ,
						   MuSeqT, LT, LabelT,
						   BestDiag, BestDiagLo, BestDiagLen);
#endif

				}
			}
		}
#endif
		for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
			{
			uint SeqIdxT = SeqIdxTs[HitIdx];
			SeqIdxTToBestDiagScore[SeqIdxT] = 0;
			}
#if DEBUG
		{
		for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
			{
			assert(SeqIdxTToBestDiagScore[SeqIdx] == 0);
			}
		}
#endif
		HitCount = 0;
		}
	CloseStdioFile(fOut);
	}
