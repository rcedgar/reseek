#include "myutils.h"
#include "alpha.h"
#include "kappa_dex.h"
#include "kappa_mermx.h"
#include "seqdb.h"
#include "quarts.h"
#include "kappa_prefilter_params.h"
#include "dssparams.h"
#include "binner.h"

uint8_t *kappa_dex::m_Offsets;
uint32_t kappa_dex::m_DictSize;
uint32_t kappa_dex::m_k;
uint32_t kappa_dex::m_K;

void kappa_dex::Init()
	{
	m_Offsets = DSSParams::m_PrefilterKappaKmerOnesOffsets;
	m_DictSize = DSSParams::m_PrefilterKappaDictSize;
	m_k = DSSParams::m_PrefilterKappaKmerNrOnes;
	m_K = DSSParams::m_PrefilterKappaKmerWidth;
	}

const uint32_t kappa_dex::m_ItemSize = 6;	// 4 byte SeqIdx + 2 byte Pos

#define TRACE	0

void kappa_dex::SetSeq(uint SeqIdx, const char *Label, const byte *Seq, uint L)
	{
#if KAPPA_DEBUG_CHECKS
	{
	for (uint i = 0; i < L; ++i)
		{
		uint Letter = Seq[i];
		if (Letter >= KAPPA_AS)
			Die("kappa_dex::SetSeq(SeqIdx=%u, %s) Pos=%u, L=%u, Mu letter=%u",
				SeqIdx, Label, i, L, Letter);
		}
	}
#endif
	m_SeqIdx = SeqIdx;
	m_Label = Label;
	m_Seq = Seq;
	m_L = L;
	GetKmers(Seq, L, m_Kmers);
#if TRACE
	LogSeq();
#endif
	}

uint kappa_dex::StrToKmer(const string &s) const
	{
	assert(SIZE(s) == m_k);
	return StrToKmer(s.c_str());
	}

uint kappa_dex::BytesToKmer(const byte *s) const
	{
	uint Kmer = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		byte Letter = s[m_Offsets[i]];
		Kmer = Kmer*KAPPA_AS + Letter;
		}
	return Kmer;
	}

uint kappa_dex::StrToKmer(const char *s) const
	{
	uint Kmer = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		Kmer *= KAPPA_AS;
		byte c = s[i];
		uint Letter = g_CharToLetterMu[c];
		asserta(Letter < KAPPA_AS);
		Kmer += Letter;
		}
	return Kmer;
	}

uint kappa_dex::GetKmerMaxLetterCount(uint Kmer)
	{
	uint8_t KmerLetterCounts[KAPPA_AS];
	memset((void *) KmerLetterCounts, 0, KAPPA_AS);
	for (uint i = 0; i < m_k; ++i)
		{
		byte Letter = Kmer%KAPPA_AS;
		KmerLetterCounts[Letter] += 1;
		Kmer /= KAPPA_AS;
		}

	uint8_t maxn = 1;
	for (uint Letter = 0; Letter < KAPPA_AS; ++Letter)
		maxn = max(maxn, KmerLetterCounts[Letter]);
	return maxn;
	}

const char *kappa_dex::KmerToStr(uint Kmer, string &s) const
	{
	s.clear();
	for (uint i = 0; i < m_k; ++i)
		{
		byte Letter = Kmer%KAPPA_AS;
		s.push_back(g_LetterToCharMu[Letter]);
		Kmer /= KAPPA_AS;
		}
	reverse(s.begin(), s.end());
	return s.c_str();
	}

void kappa_dex::Alloc_Pass1()
	{
	if (m_AddNeighborhood && m_NeighborKmers == 0)
		m_NeighborKmers = myalloc(uint, DSSParams::m_PrefilterKappaDictSize);

// Pass1 m_Finger[Kmer] = Count
	asserta(m_Finger == 0 && m_Data == 0);
	m_Finger = myalloc(uint32_t, m_DictSize + 2);
	zero_array(m_Finger, m_DictSize+2);
#if KAPPA_DEBUG_CHECKS
	m_KmerToCount1.resize(m_DictSize, 0);
#endif
	}

void kappa_dex::Alloc_Pass2()
	{
// 6 bytes for uint32_t:uint16_t (SeqIdx:Pos)
	const uint64_t Bytes = uint64(m_ItemSize)*m_Size;
	m_Data = myalloc64(uint8_t, Bytes);
#if KAPPA_DEBUG_CHECKS
	m_KmerToCount2.resize(m_DictSize, 0);
	memset(m_Data, 0xff, Bytes);
#endif
	}

void  kappa_dex::AddSeq_Pass1()
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
#if KAPPA_DEBUG_CHECKS
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
			   DSSParams::m_PrefilterMinMuKmerPairScore, m_NeighborKmers);
			for (uint j = 0; j < n; ++j)
				{
				uint NeighborKmer = m_NeighborKmers[j];
				asserta(NeighborKmer < DSSParams::m_PrefilterKappaDictSize);
				asserta(m_Size < UINT_MAX);
				asserta(m_Finger[NeighborKmer+1] < UINT_MAX);
				m_Finger[NeighborKmer+1] += 1;
				++m_Size;
#if KAPPA_DEBUG_CHECKS
				m_KmerToCount1[NeighborKmer] += 1;
#endif
				}
			}
		}
	}

void kappa_dex::AddSeq_Pass2()
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
#if KAPPA_DEBUG_CHECKS
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
			   DSSParams::m_PrefilterMinMuKmerPairScore, m_NeighborKmers);
			for (uint j = 0; j < n; ++j)
				{
				uint NeighborKmer = m_NeighborKmers[j];
				asserta(NeighborKmer < DSSParams::m_PrefilterKappaDictSize);
				uint DataOffset = m_Finger[NeighborKmer+1];
				Put(DataOffset, m_SeqIdx, SeqPos);
				asserta(m_Finger[NeighborKmer+1] < UINT_MAX);
				m_Finger[NeighborKmer+1] += 1;
#if KAPPA_DEBUG_CHECKS
				assert(m_KmerToDataStart[NeighborKmer] + 
					   m_KmerToCount2[NeighborKmer] == DataOffset);
				m_KmerToCount2[NeighborKmer] += 1;
#endif
				}
			}
		}
	}

void kappa_dex::LogStats() const
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

#if KAPPA_DEBUG_CHECKS
void kappa_dex::CheckAfterPass1() const
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

void kappa_dex::CheckAfterAdjust() const
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

void kappa_dex::CheckAfterPass2() const
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

void kappa_dex::AdjustFinger()
	{
	uint Sum = 0;
	for (uint Kmer = 0; Kmer <= m_DictSize; ++Kmer)
		{
#if KAPPA_DEBUG_CHECKS
		m_KmerToDataStart.push_back(Sum);
#endif
		uint Kmer_Size = m_Finger[Kmer+1];
		asserta(m_Finger[Kmer+1] < UINT_MAX);
		m_Finger[Kmer+1] = Sum;
		Sum += Kmer_Size;
		}
	asserta(Sum == m_Size);
	}

void kappa_dex::Validate() const
	{
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		ValidateKmer(Kmer);
	}

void kappa_dex::LogIndexKmer(uint Kmer) const
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

void kappa_dex::ValidateKmer(uint Kmer) const
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
			Die("kappa_dex::ValidateKmer(Kmer=0x%x)", Kmer);
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

uint kappa_dex::GetSeqKmer(const byte *Seq, uint SeqPos, bool SelfScoreMask) const
	{
	uint Kmer = BytesToKmer(Seq + SeqPos);
	if (SelfScoreMask && m_KmerSelfScores[Kmer] < m_MinKmerSelfScore)
		Kmer = UINT_MAX;
	return Kmer;
	}

void kappa_dex::FromSeqDB(const SeqDB &Input)//TODO FromBags already have Mu k-mers
	{
	m_SeqDB = &Input;
	const uint SeqCount = Input.GetSeqCount();
	if (m_AddNeighborhood && m_ptrScoreMx == 0)
		m_ptrScoreMx = &GetMuMerMx(m_k);

	Alloc_Pass1();
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, SeqCount, "kappa_dex pass 1");
		const char *Label = m_SeqDB->GetLabel(SeqIdx).c_str();
		const byte *Seq = Input.GetByteSeq(SeqIdx);
		const uint L = Input.GetSeqLength(SeqIdx);
		SetSeq(SeqIdx, Label, Seq, L);
		AddSeq_Pass1();
		}
#if KAPPA_DEBUG_CHECKS
	CheckAfterPass1();
#endif

	AdjustFinger();
#if KAPPA_DEBUG_CHECKS
	CheckAfterAdjust();
#endif

	Alloc_Pass2();
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, SeqCount, "kappa_dex pass 2");
		const char *Label = m_SeqDB->GetLabel(SeqIdx).c_str();
		const byte *Seq = Input.GetByteSeq(SeqIdx);
		const uint L = Input.GetSeqLength(SeqIdx);
		SetSeq(SeqIdx, Label, Seq, L);
		AddSeq_Pass2();
		}
	SetRowSizes();
#if KAPPA_DEBUG_CHECKS
	CheckAfterPass2();
#endif

#if KAPPA_DEBUG_CHECKS
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

void kappa_dex::SetRowSizes()
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

void kappa_dex::Put(uint DataOffset, uint32_t SeqIdx, uint16_t SeqPos)
	{
	assert(DataOffset < m_Size);
	uint64 Bytes64 = uint64(m_ItemSize)*uint64(DataOffset);
	uint8_t *ptr = m_Data + Bytes64;
	*(uint32_t *) ptr = SeqIdx;
	*(uint16_t *) (ptr + 4) = SeqPos;
#if KAPPA_DEBUG_CHECKS
	{
	uint32_t Check_SeqIdx;
	uint16_t Check_SeqPos;
	Get(DataOffset, Check_SeqIdx, Check_SeqPos);
	assert(Check_SeqIdx == SeqIdx);
	assert(Check_SeqPos == SeqPos);
	}
#endif
	}

void kappa_dex::Get(uint DataOffset, uint32_t &SeqIdx, uint16_t &SeqPos) const
	{
	const uint8_t *ptr = m_Data + m_ItemSize*uint64(DataOffset);
	SeqIdx = *(uint32_t *) ptr;
	SeqPos = *(uint16_t *) (ptr + 4);
	}

void kappa_dex::GetKmersAndSizes(const byte *Seq, uint L,
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
			Kmer = Kmer*KAPPA_AS + Letter;
			}
#if KAPPA_DEBUG_CHECKS
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

void kappa_dex::GetKmers(const byte *Seq, uint L, vector<uint> &Kmers) const
	{
#if DEBUG
	{
	for (uint i = 0; i < L; ++i)
		assert(Seq[i] < KAPPA_AS);
	}
#endif
	Kmers.reserve(L);
	Kmers.clear();
	for (uint KmerStartPos = 0; KmerStartPos + m_K <= L; ++KmerStartPos)
		{
		uint Kmer = 0;
		for (uint i = 0; i < m_k; ++i)
			{
			uint off = m_Offsets[i];
			assert(KmerStartPos + off < L);
			byte Letter = Seq[KmerStartPos + off];
			assert(Letter < KAPPA_AS);
			Kmer = Kmer*KAPPA_AS + Letter;
			}
#if KAPPA_DEBUG_CHECKS
		uint CheckKmer = GetSeqKmer(Seq, KmerStartPos, false);
		asserta(CheckKmer == Kmer);
#endif
		assert(Kmer < DSSParams::m_PrefilterKappaDictSize);
		if (m_KmerSelfScores != 0 && m_KmerSelfScores[Kmer] < m_MinKmerSelfScore)
			Kmers.push_back(UINT_MAX);
		else
			Kmers.push_back(Kmer%m_DictSize);
		}
	}
