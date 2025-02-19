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

	m_Overflows = myalloc(uint32_t *, m_NrRdxs);
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
		}
	else if (Size%m_FixedItemsPerRdx == 0)
		{
		// Overflow, needs new data
		uint32_t *OldOverflow = m_Overflows[Rdx];
		uint32_t *NewOverflow = myalloc(uint32_t, Size + m_FixedItemsPerRdx);
		if (Size == m_FixedItemsPerRdx)
			{
			// First overflow
			assert(OldOverflow == 0);
			m_Overflows[Rdx] = OldOverflow;
			}
		else
			{
			// Expand overflow buffer
			memcpy(NewOverflow, OldOverflow, Size*sizeof(uint32_t));
			myfree(OldOverflow);
			}
		m_Overflows[Rdx] = NewOverflow;
		PutRaw(NewOverflow + Size - m_FixedItemsPerRdx, SeqIdx, Diag);
		}
	else
		{
		// Previous overflow, fits into current buffer
		uint32_t *Overflow = m_Overflows[Rdx];
		assert(Overflow != 0);
		PutRaw(Overflow + Size - m_FixedItemsPerRdx, SeqIdx, Diag);
		}

	m_Sizes[Rdx] = Size + 1;
	m_Size += 1;
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

	uint32_t *Overflow = m_Overflows[Rdx];
	asserta(Overflow != 0);

	uint Sum = 0;
	uint BlockIdx = 0;
	uint Remainder = Size - m_FixedItemsPerRdx;
	for (uint i = 0; i < Remainder; ++i)
		{
		const uint32_t *ptr = Overflow + i;
		uint16_t Diag;
		uint32_t SeqIdx = GetAllBits(Rdx, ptr, Diag);
		Log(" %u:%u", SeqIdx, Diag);
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
		uint Size = m_Sizes[Rdx];
		assert(m_Size > 0);
		if (Size > m_FixedItemsPerRdx)
			{
			uint32_t *Overflows = m_Overflows[Rdx];
			assert(Overflows != 0);
			myfree(m_Overflows[Rdx]);
			m_Overflows[Rdx] = 0;
			}
		else
			assert(m_Overflows[Rdx] == 0);
		m_Sizes[Rdx] = 0;
		m_MaxSeqIdxHiBits[Rdx] = 0;
		m_BusyRdxs[i] = 0;
		}
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

	const uint32_t *Overflow = m_Overflows[Rdx];
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
	uint Remainder = Size - m_FixedItemsPerRdx;
	for (uint i = 0; i < Remainder; ++i)
		{
		const uint32_t *ptr = Overflow + i;
		uint16_t Diag;
		uint32_t SeqIdx = GetAllBits(Rdx, ptr, Diag);
		asserta(SeqIdx <= MaxSeqIdx);
		asserta(Diag <= MaxDiag);

		uint32_t SeqIdxHiBits = SeqIdx/m_SeqMod;
		MaxSeqIdxHiBits = max(SeqIdxHiBits, MaxSeqIdxHiBits);
		}
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

	const uint32_t *Overflow = m_Overflows[Rdx];
	asserta(Overflow != 0);

	uint Sum = 0;
	uint Remainder = Size - m_FixedItemsPerRdx;
	for (uint i = 0; i < Remainder; ++i)
		{
		const uint32_t *Data = Overflow + i;
		uint16_t Diag;
		uint32_t SeqIdx = GetAllBits(Rdx, Data, Diag);
		SeqIdxDiagPairs.push_back(pair<uint32_t, uint16_t>(SeqIdx, Diag));
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

	uint32_t *Overflow = m_Overflows[Rdx];
	assert(Overflow != 0);
	uint Remainder = Size - m_FixedItemsPerRdx;
	for (uint i = 0; i < Remainder; ++i)
		D.Add(Overflow[i]);
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
	delete &D;
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
		uint16_t Diag = uint16_t(randu32()%m_DiagHi);
		uint Offset = randu32()%(m_TotalFixedItems - 32);
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
		uint16_t MaxDiag = (123 + randu32()%63001) & m_Mask14;
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

void cmd_twohittest()
	{
	SimpleTests();
	SavedTests();
	ProgressLog("Tests completed ok\n");
	}
