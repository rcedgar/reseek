#include "myutils.h"
#include "alpha.h"
#include "mudex.h"
#include "mermx.h"
#include "seqdb.h"
#include "quarts.h"
#include "prefiltermuparams.h"
#include "binner.h"

const MerMx &GetMuMerMx(uint k);

const uint8_t *MuDex::m_Offsets = s_Offsets;
const uint32_t MuDex::m_DictSize = 36*36*36*36*36; // 60466176 (60.4 M)
const uint32_t MuDex::m_k = k;
const uint32_t MuDex::m_K = K;

/***
32 bits 2^5, 64 bits 2^6
5 bits per aa letter
6 bits per Mu letter
8 bits per char
***/
const uint32_t MuDex::m_ItemSize = 6;	// 4 byte SeqIdx + 2 byte Pos

#define TRACE	0

void MuDex::SetSeq(uint SeqIdx, const char *Label, const byte *Seq, uint L)
	{
	m_SeqIdx = SeqIdx;
	m_Label = Label;
	m_Seq = Seq;
	m_L = L;
	GetKmers(Seq, L, m_Kmers);
#if TRACE
	LogSeq();
#endif
	}

uint MuDex::StrToKmer(const string &s) const
	{
	assert(SIZE(s) == m_k);
	return StrToKmer(s.c_str());
	}

uint MuDex::BytesToKmer(const byte *s) const
	{
	uint Kmer = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		byte Letter = s[m_Offsets[i]];
		Kmer = Kmer*36 + Letter;
		}
	return Kmer;
	}

uint MuDex::StrToKmer(const char *s) const
	{
	uint Kmer = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		Kmer *= 36;
		byte c = s[i];
		uint Letter = g_CharToLetterMu[c];
		asserta(Letter < 36);
		Kmer += Letter;
		}
	return Kmer;
	}

uint MuDex::GetKmerMaxLetterCount(uint Kmer)
	{
	uint8_t KmerLetterCounts[36];
	memset((void *) KmerLetterCounts, 0, 36);
	for (uint i = 0; i < m_k; ++i)
		{
		byte Letter = Kmer%36;
		KmerLetterCounts[Letter] += 1;
		Kmer /= 36;
		}

	uint8_t maxn = 1;
	for (uint Letter = 0; Letter < 36; ++Letter)
		maxn = max(maxn, KmerLetterCounts[Letter]);
	return maxn;
	}

const char *MuDex::KmerToStr(uint Kmer, string &s) const
	{
	s.clear();
	for (uint i = 0; i < m_k; ++i)
		{
		byte Letter = Kmer%36;
		s.push_back(g_LetterToCharMu[Letter]);
		Kmer /= 36;
		}
	reverse(s.begin(), s.end());
	return s.c_str();
	}

void MuDex::Alloc_Pass1()
	{
	if (m_AddNeighborhood && m_NeighborKmers == 0)
		m_NeighborKmers = myalloc(uint, MAX_HOOD_SIZE);

// Pass1 m_Finger[Kmer] = Count
	asserta(m_Finger == 0 && m_Data == 0);
	m_Finger = myalloc(uint32_t, m_DictSize + 2);
	zero_array(m_Finger, m_DictSize+2);
#if DEBUG_CHECKS
	m_KmerToCount1.resize(m_DictSize, 0);
#endif
	}

void MuDex::Alloc_Pass2()
	{
// 6 bytes for uint32_t:uint16_t (SeqIdx:Pos)
	const uint64_t Bytes = uint64(m_ItemSize)*m_Size;
	m_Data = myalloc64(uint8_t, Bytes);
#if DEBUG_CHECKS
	m_KmerToCount2.resize(m_DictSize, 0);
	memset(m_Data, 0xff, Bytes);
#endif
	}

void  MuDex::AddSeq_Pass1()
	{
#if TRACE
	string Tmp;
	Log("AddSeq_Pass1(%s) L=%u\n", m_Label, m_L);
#endif
	const uint KmerCount = SIZE(m_Kmers);
	for (uint SeqPos = 0; SeqPos < KmerCount; ++SeqPos)
		{
		uint Kmer = m_Kmers[SeqPos];
		if (Kmer == UINT_MAX)
			{
#if TRACE
			Log("[%4u] ***\n", SeqPos);
#endif
			continue;
			}

	// Pass 1, m_Finger[Kmer+1] is count
		asserta(m_Size < UINT_MAX);
		asserta(m_Finger[Kmer+1] < UINT_MAX);
		m_Finger[Kmer+1] += 1;
		++m_Size;
#if DEBUG_CHECKS
		m_KmerToCount1[Kmer] += 1;
#endif
#if TRACE
		Log("[%4u] %08x %s", SeqPos, Kmer, KmerToStr(Kmer, Tmp));
		if (m_KmerSelfScores != 0)
			Log(" self=%d", m_KmerSelfScores[Kmer]);
		Log("\n");
#endif

		if (m_AddNeighborhood)
			{
			uint n = m_ptrScoreMx->GetHighScoringKmers(Kmer, 
			  MIN_KMER_PAIR_SCORE, m_NeighborKmers);
			for (uint j = 0; j < n; ++j)
				{
				uint NeighborKmer = m_NeighborKmers[j];
				asserta(NeighborKmer < DICT_SIZE);
				asserta(m_Size < UINT_MAX);
				asserta(m_Finger[NeighborKmer+1] < UINT_MAX);
				m_Finger[NeighborKmer+1] += 1;
				++m_Size;
#if DEBUG_CHECKS
				m_KmerToCount1[NeighborKmer] += 1;
#endif
				}
			}
		}
	}

void MuDex::AddSeq_Pass2()
	{
#if TRACE
	string Tmp;
	Log("AddSeq_Pass2(%s) L=%u\n", m_Label, m_L);
#endif
	const uint KmerCount = SIZE(m_Kmers);
	for (uint SeqPos = 0; SeqPos < KmerCount; ++SeqPos)
		{
		uint Kmer = m_Kmers[SeqPos];
		if (Kmer == UINT_MAX)
			{
#if TRACE
			Log("[%4u] %08x %s --LOW\n", SeqPos, Kmer, KmerToStr(Kmer, Tmp));
#endif
			continue;
			}
		uint DataOffset = m_Finger[Kmer+1];
		Put(DataOffset, m_SeqIdx, SeqPos);
		asserta(m_Finger[Kmer+1] < UINT_MAX);
		m_Finger[Kmer+1] += 1;
#if DEBUG_CHECKS
		assert(m_KmerToDataStart[Kmer] + m_KmerToCount2[Kmer] == DataOffset);
		m_KmerToCount2[Kmer] += 1;
#endif
#if TRACE
		Log("[%4u] %08x %s DO=%u\n",
			SeqPos, Kmer, KmerToStr(Kmer, Tmp), DataOffset);
#endif
		if (m_AddNeighborhood)
			{
			uint n = m_ptrScoreMx->GetHighScoringKmers(Kmer, 
			  MIN_KMER_PAIR_SCORE, m_NeighborKmers);
			for (uint j = 0; j < n; ++j)
				{
				uint NeighborKmer = m_NeighborKmers[j];
				asserta(NeighborKmer < DICT_SIZE);
				uint DataOffset = m_Finger[NeighborKmer+1];
				Put(DataOffset, m_SeqIdx, SeqPos);
				asserta(m_Finger[NeighborKmer+1] < UINT_MAX);
				m_Finger[NeighborKmer+1] += 1;
#if DEBUG_CHECKS
				assert(m_KmerToDataStart[NeighborKmer] + 
					   m_KmerToCount2[NeighborKmer] == DataOffset);
				m_KmerToCount2[NeighborKmer] += 1;
#endif
				}
			}
		}
	}

void MuDex::LogStats() const
	{
	vector<uint> RowSizes;
	uint Sum = 0;
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint Size = GetRowSize(Kmer);
		Sum += Size;
		RowSizes.push_back(Size);
		}
	Quarts Q;
	GetQuarts(RowSizes, Q);
	Log("RowSizes: ");
	Q.LogMe();
	Log("Total = %u (%s)\n", Sum, IntToStr(Sum));
	}

#if DEBUG_CHECKS
void MuDex::CheckAfterPass1() const
	{
	uint Check_Size = 0;
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint n = m_Finger[Kmer+1];
		uint Check_n = m_KmerToCount1[Kmer];
		if (Check_n != n)
			{
			Log("Kmer %08x DictSize %08x Check_n %u n %u\n",
				Kmer, m_DictSize, Check_n, n);
			Die("CheckAfterPass1");
			}
		Check_Size += n;
		}
	assert(Check_Size == m_Size);
	ProgressLog("CheckAfterPass1 OK\n");
	}

void MuDex::CheckAfterAdjust() const
	{
	assert(m_Finger[m_DictSize+1] == m_Size);
	uint Check_Size = 0;
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint n = m_Finger[Kmer+2] - m_Finger[Kmer+1];
		uint Check_n = m_KmerToCount1[Kmer];
		if (Check_n != n)
			{
			Log("Kmer %08x DictSize %08x Check_n %u n %u\n",
				Kmer, m_DictSize, Check_n, n);
			Die("CheckAfterAdjust");
			}
		Check_Size += n;
		}
	assert(Check_Size == m_Size);
	ProgressLog("CheckAfterAdjust OK\n");
	}

void MuDex::CheckAfterPass2() const
	{
	uint Check_Size = 0;
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint n = m_Finger[Kmer+1] - m_Finger[Kmer];
		uint Check_n1 = m_KmerToCount1[Kmer];
		uint Check_n2 = m_KmerToCount2[Kmer];
		if (Check_n1 != n || Check_n2 != n)
			{
			Log("Kmer %08x DictSize %08x Check_n1 %u Check_n2 %u n %u\n",
				Kmer, m_DictSize, Check_n1, Check_n2, n);
			Die("CheckAfterPass2");
			}
		Check_Size += n;
		}
	assert(Check_Size == m_Size);
	ProgressLog("CheckAfterPass2 OK\n");
	}
#endif

void MuDex::AdjustFinger()
	{
	uint Sum = 0;
	for (uint Kmer = 0; Kmer <= m_DictSize; ++Kmer)
		{
#if DEBUG_CHECKS
		m_KmerToDataStart.push_back(Sum);
#endif
		uint Kmer_Size = m_Finger[Kmer+1];
		asserta(m_Finger[Kmer+1] < UINT_MAX);
		m_Finger[Kmer+1] = Sum;
		Sum += Kmer_Size;
		}
	asserta(Sum == m_Size);
	}

void MuDex::Validate() const
	{
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		ValidateKmer(Kmer);
	}

void MuDex::LogIndexKmer(uint Kmer) const
	{
	uint n = GetRowSize(Kmer);
	string Tmp;
	uint DataOffset = m_Finger[Kmer];
	Log("LogIndexKmer(%08x) %s size=%u DO=%u",
		Kmer, KmerToStr(Kmer, Tmp), n, DataOffset);
	for (uint i = 0; i < n; ++i)
		{
		uint32_t SeqIdx;
		uint16_t SeqPos;
		Get(DataOffset+i, SeqIdx, SeqPos);
		Log(" %u:%u", SeqIdx, SeqPos);
		//uint Check_Kmer = GetSeqKmer(SeqIdx, SeqPos);
		//asserta(Check_Kmer == Kmer);
		}
	Log("\n");
	}

void MuDex::ValidateKmer(uint Kmer) const
	{
	const uint QSeqCount = m_SeqDB->GetSeqCount();
	uint n = GetRowSize(Kmer);
	uint DataOffset = m_Finger[Kmer];
	asserta(DataOffset <= m_Size);
	for (uint i = 0; i < n; ++i)
		{
		uint32_t SeqIdx;
		uint16_t SeqPos;
		Get(DataOffset, SeqIdx, SeqPos);
		if (SeqIdx >= QSeqCount)
			{
			Log("m_Size = %s\n", Int64ToStr(m_Size));
			Log("m_Finger[0x%x] = %u\n", Kmer, m_Finger[Kmer]);
			Log("i=%u n=%u\n", i, n);
			Log("QSeqIdx=%u, SeqPos=%u\n", SeqIdx, SeqPos);
			Die("MuDex::ValidateKmer(Kmer=0x%x)", Kmer);
			}
		uint QL = m_SeqDB->GetSeqLength(SeqIdx);
		asserta(SeqPos < QL);
		if (!m_AddNeighborhood)
			{
			const byte *Seq = m_SeqDB->GetByteSeq(SeqIdx);
			uint Check_Kmer = GetSeqKmer(Seq, SeqPos, false);
			asserta(Check_Kmer == Kmer);
			}
		}
	}

uint MuDex::GetSeqKmer(const byte *Seq, uint SeqPos, bool SelfScoreMask) const
	{
	uint Kmer = BytesToKmer(Seq + SeqPos);
	if (SelfScoreMask && m_KmerSelfScores[Kmer] < m_MinKmerSelfScore)
		Kmer = UINT_MAX;
	return Kmer;
	}

void MuDex::FromSeqDB(const SeqDB &Input)
	{
	m_SeqDB = &Input;
	const uint SeqCount = Input.GetSeqCount();
	if (m_AddNeighborhood && m_ptrScoreMx == 0)
		m_ptrScoreMx = &GetMuMerMx(m_k);

	Alloc_Pass1();
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, SeqCount, "MuDex pass 1");
		const char *Label = m_SeqDB->GetLabel(SeqIdx).c_str();
		const byte *Seq = Input.GetByteSeq(SeqIdx);
		const uint L = Input.GetSeqLength(SeqIdx);
		SetSeq(SeqIdx, Label, Seq, L);
		AddSeq_Pass1();
		}
#if DEBUG_CHECKS
	CheckAfterPass1();
#endif

	AdjustFinger();
#if DEBUG_CHECKS
	CheckAfterAdjust();
#endif

	Alloc_Pass2();
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, SeqCount, "MuDex pass 2");
		const char *Label = m_SeqDB->GetLabel(SeqIdx).c_str();
		const byte *Seq = Input.GetByteSeq(SeqIdx);
		const uint L = Input.GetSeqLength(SeqIdx);
		SetSeq(SeqIdx, Label, Seq, L);
		AddSeq_Pass2();
		}
	SetRowSizes();
#if DEBUG_CHECKS
	CheckAfterPass2();
#endif

#if DEBUG_CHECKS
	{
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint RowSize = GetRowSize(Kmer);
		uint Check_RowSize = m_KmerToCount1[Kmer];
		asserta(Check_RowSize == RowSize);
		if (RowSize == 0)
			continue;
		uint Offset = m_Finger[Kmer];
		uint Check_Offset = m_KmerToDataStart[Kmer];
		}
	}
	Validate();
#endif
	}

void MuDex::SetRowSizes()
	{
	//m_RowSizes = myalloc(uint16_t, m_DictSize);
	m_RowSizes = myalloc(uint32_t, m_DictSize);
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint32_t RowSize32 = m_Finger[Kmer+1] - m_Finger[Kmer];
		//uint16_t RowSize = uint16_t(RowSize32);
		//assert(uint32_t(RowSize) == RowSize32);
		//m_RowSizes[Kmer] = RowSize;
		m_RowSizes[Kmer] = RowSize32;
		}
	}

void MuDex::Put(uint DataOffset, uint32_t SeqIdx, uint16_t SeqPos)
	{
	assert(DataOffset < m_Size);
	uint64 Bytes64 = uint64(m_ItemSize)*uint64(DataOffset);
	uint8_t *ptr = m_Data + Bytes64;
	*(uint32_t *) ptr = SeqIdx;
	*(uint16_t *) (ptr + 4) = SeqPos;
#if DEBUG_CHECKS
	{
	uint32_t Check_SeqIdx;
	uint16_t Check_SeqPos;
	Get(DataOffset, Check_SeqIdx, Check_SeqPos);
	assert(Check_SeqIdx == SeqIdx);
	assert(Check_SeqPos == SeqPos);
	}
#endif
	}

void MuDex::Get(uint DataOffset, uint32_t &SeqIdx, uint16_t &SeqPos) const
	{
	const uint8_t *ptr = m_Data + m_ItemSize*uint64(DataOffset);
	SeqIdx = *(uint32_t *) ptr;
	SeqPos = *(uint16_t *) (ptr + 4);
	}

void MuDex::GetKmersAndSizes(const byte *Seq, uint L,
							 vector<uint> &Kmers, vector<uint> &Sizes) const
	{
	Kmers.reserve(L);
	Sizes.reserve(L);
	Kmers.clear();
	Sizes.clear();
	for (uint KmerStartPos = 0; KmerStartPos + m_K <= L; ++KmerStartPos)
		{
		uint Kmer = 0;
		for (uint i = 0; i < m_k; ++i)
			{
			byte Letter = Seq[KmerStartPos + m_Offsets[i]];
			Kmer = Kmer*36 + Letter;
			}
#if DEBUG_CHECKS
		uint CheckKmer = GetSeqKmer(Seq, KmerStartPos, false);
		asserta(CheckKmer == Kmer);
#endif
		if (m_KmerSelfScores != 0 && m_KmerSelfScores[Kmer] < m_MinKmerSelfScore)
			{
			Kmers.push_back(UINT_MAX);
			Sizes.push_back(UINT_MAX);
			}
		else
			{
			Kmer = Kmer%m_DictSize;
			uint Size = GetRowSize(Kmer);
			Kmers.push_back(Kmer);
			Sizes.push_back(Kmer);
			}
		}
	}

void MuDex::GetKmers(const byte *Seq, uint L, vector<uint> &Kmers) const
	{
	Kmers.reserve(L);
	Kmers.clear();
	for (uint KmerStartPos = 0; KmerStartPos + m_K <= L; ++KmerStartPos)
		{
		uint Kmer = 0;
		for (uint i = 0; i < m_k; ++i)
			{
			byte Letter = Seq[KmerStartPos + m_Offsets[i]];
			Kmer = Kmer*36 + Letter;
			}
#if DEBUG_CHECKS
		uint CheckKmer = GetSeqKmer(Seq, KmerStartPos, false);
		asserta(CheckKmer == Kmer);
#endif
		if (m_KmerSelfScores != 0 && m_KmerSelfScores[Kmer] < m_MinKmerSelfScore)
			Kmers.push_back(UINT_MAX);
		else
			Kmers.push_back(Kmer%m_DictSize);
		}
	}

void cmd_mudex()
	{
	SeqDB Input;
	Input.FromFasta(g_Arg1);
	Input.ToLetters(g_CharToLetterMu);

	const MerMx &ScoreMx = GetMuMerMx(k);
	asserta(ScoreMx.m_k == k);

	MuDex MD;
	MD.FromSeqDB(Input);
	MD.m_KmerSelfScores = ScoreMx.BuildSelfScores_Kmers();
	MD.LogStats();
	MD.Validate();
	ProgressLog("Validate OK\n");

	vector<uint> v;
	for (uint i = 0; i < MuDex::m_DictSize; ++i)
		{
		uint score = MD.m_KmerSelfScores[i];
		v.push_back(score);
		}
	Quarts Q;
	GetQuarts(v, Q);
// k=5 SelfScores: 
// N=60466176, Min=20, LoQ=43, Med=47, HiQ=51, Max=75, Avg=47.3611
	Log("SelfScores: ");
	Q.LogMe();

	vector<uint> Counts(6);
	uint Total = 0;
	for (uint Kmer = 0; Kmer < MuDex::m_DictSize; ++Kmer)
		{
		uint n = MD.GetKmerMaxLetterCount(Kmer);
		uint Size = MD.GetRowSize(Kmer);
		Total += Size;
		Counts[n] += Size;
		}

// SCOP40
//Max letters [1] = 237643 (12.6%)
//Max letters [2] = 1023722 (54.4%)
//Max letters [3] = 463046 (24.6%)
//Max letters [4] = 112057 (6.0%)
//Max letters [5] = 44543 (2.4%)

// Dictionary
//Max letters [1] = 45239040 (74.8%)
//Max letters [2] = 14779800 (24.4%)
//Max letters [3] = 441000 (0.7%)
//Max letters [4] = 6300 (0.0%)
//Max letters [5] = 36 (0.0%)
	for (uint i = 1; i <= 5; ++i)
		{
		uint n = Counts[i];
		double Pct = GetPct(n, Total);
		ProgressLog("Max letters [%u] = %u (%.1f%%)\n",
					i, n, Pct);
		}
	}
