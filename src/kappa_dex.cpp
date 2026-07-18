#include "myutils.h"
#include "alpha.h"
#include "kappa_dex.h"
#include "kappa_mermx.h"
#include "seqdb.h"
#include "quarts.h"
#include "flat_params.h"
#include "binner.h"
#include <atomic>
#include <mutex>
#include <thread>
#include <vector>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

static inline uint64_t AtomicFetchAddU64(uint64_t *p, uint64_t add)
	{
#if defined(__GNUC__) || defined(__clang__)
	return __atomic_fetch_add(p, add, __ATOMIC_RELAXED);
#elif defined(_MSC_VER)
	return (uint64_t)_InterlockedExchangeAdd64(
		reinterpret_cast<volatile long long *>(p),
		(long long)add);
#else
	return std::atomic_ref<uint64_t>(*p).fetch_add(add, std::memory_order_relaxed);
#endif
	}

uint8_t *kappa_dex::m_Offsets;
uint32_t kappa_dex::m_DictSize;
uint32_t kappa_dex::m_k;
uint32_t kappa_dex::m_K;

void kappa_dex::Init()
	{
	m_Offsets = flat_params::m_kappa_kmer_onesoffsets;
	m_DictSize = flat_params::m_kappa_dict_size;
	m_k = flat_params::m_kappa_kmer_nrones;
	m_K = flat_params::m_kappa_kmer_width;
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
		m_NeighborKmers = myalloc(uint, flat_params::m_kappa_dict_size);

// Pass1 m_Finger[Kmer] = Count
	asserta(m_Finger == 0 && m_Data == 0);
	m_Finger = myalloc(uint64_t, m_DictSize + 2);
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
			   flat_params::m_kappa_min_kmerpairscore, m_NeighborKmers);
			for (uint j = 0; j < n; ++j)
				{
				uint NeighborKmer = m_NeighborKmers[j];
				asserta(NeighborKmer < flat_params::m_kappa_dict_size);
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
		uint64_t DataOffset = m_Finger[Kmer+1];
		Put(DataOffset, m_SeqIdx, SeqPos);
		m_Finger[Kmer+1] += 1;
#if KAPPA_DEBUG_CHECKS
		assert(m_KmerToDataStart[Kmer] + m_KmerToCount2[Kmer] == DataOffset);
		m_KmerToCount2[Kmer] += 1;
#endif
#if TRACE
		Log("[%4u] %08x %s DO=%s\n",
			SeqPos, Kmer, KmerToStr(Kmer, Tmp), Int64ToStr(DataOffset));
#endif
		if (m_AddNeighborhood)
			{
			uint n = m_ptrScoreMx->GetHighScoringKmers(Kmer, 
			   flat_params::m_kappa_min_kmerpairscore, m_NeighborKmers);
			for (uint j = 0; j < n; ++j)
				{
				uint NeighborKmer = m_NeighborKmers[j];
				asserta(NeighborKmer < flat_params::m_kappa_dict_size);
				uint64_t DataOffset = m_Finger[NeighborKmer+1];
				Put(DataOffset, m_SeqIdx, SeqPos);
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
	uint64_t Sum = 0;
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint64_t Size = GetRowSize(Kmer);
		Sum += Size;
		RowSizes.push_back(Size > UINT_MAX ? UINT_MAX : uint(Size));
		}
	Quarts Q;
	GetQuarts(RowSizes, Q);
	Log("RowSizes: ");
	Q.LogMe();
	Log("Total = %s\n", Int64ToStr(Sum));
	}

#if KAPPA_DEBUG_CHECKS
void kappa_dex::CheckAfterPass1() const
	{
	uint64_t Check_Size = 0;
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint64_t n = m_Finger[Kmer+1];
		uint64_t Check_n = m_KmerToCount1[Kmer];
		if (Check_n != n)
			{
			Log("Kmer %08x DictSize %08x Check_n %" PRIu64 " n %" PRIu64 "\n",
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
	uint64_t Check_Size = 0;
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint64_t n = m_Finger[Kmer+2] - m_Finger[Kmer+1];
		uint64_t Check_n = m_KmerToCount1[Kmer];
		if (Check_n != n)
			{
			Log("Kmer %08x DictSize %08x Check_n %" PRIu64 " n %" PRIu64 "\n",
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
	uint64_t Check_Size = 0;
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint64_t n = m_Finger[Kmer+1] - m_Finger[Kmer];
		uint64_t Check_n1 = m_KmerToCount1[Kmer];
		uint64_t Check_n2 = m_KmerToCount2[Kmer];
		if (Check_n1 != n || Check_n2 != n)
			{
			Log("Kmer %08x DictSize %08x Check_n1 %" PRIu64 " Check_n2 %" PRIu64 " n %" PRIu64 "\n",
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
	uint64_t Sum = 0;
	for (uint Kmer = 0; Kmer <= m_DictSize; ++Kmer)
		{
#if KAPPA_DEBUG_CHECKS
		m_KmerToDataStart.push_back(Sum);
#endif
		uint64_t Kmer_Size = m_Finger[Kmer+1];
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
	uint64_t n = GetRowSize(Kmer);
	string Tmp;
	uint64_t DataOffset = m_Finger[Kmer];
	Log("LogIndexKmer(%08x) %s size=%s DO=%s",
		Kmer, KmerToStr(Kmer, Tmp), Int64ToStr(n), Int64ToStr(DataOffset));
	for (uint64_t i = 0; i < n; ++i)
		{
		uint32_t SeqIdx;
		uint16_t SeqPos;
		Get(DataOffset+i, SeqIdx, SeqPos);
		Log(" %u:%u", SeqIdx, SeqPos);
		}
	Log("\n");
	}

void kappa_dex::ValidateKmer(uint Kmer) const
	{
	const uint QSeqCount = m_nseq;
	asserta(QSeqCount > 0);
	uint64_t n = GetRowSize(Kmer);
	uint64_t DataOffset = m_Finger[Kmer];
	asserta(DataOffset <= m_Size);
	for (uint64_t i = 0; i < n; ++i)
		{
		uint32_t SeqIdx;
		uint16_t SeqPos;
		Get(DataOffset, SeqIdx, SeqPos);
		if (SeqIdx >= QSeqCount)
			{
			Log("m_Size = %s\n", Int64ToStr(m_Size));
			Log("m_Finger[0x%x] = %s\n", Kmer, Int64ToStr(m_Finger[Kmer]));
			Log("i=%s n=%s\n", Int64ToStr(i), Int64ToStr(n));
			Log("QSeqIdx=%u, SeqPos=%u\n", SeqIdx, SeqPos);
			Die("kappa_dex::ValidateKmer(Kmer=0x%x)", Kmer);
			}
		uint QL = get_seq_length(SeqIdx);
		asserta(SeqPos < QL);
		if (!m_AddNeighborhood)
			{
			const byte *Seq = get_byte_seq(SeqIdx);
			uint Check_Kmer = GetSeqKmer(Seq, SeqPos, false);
			asserta(Check_Kmer == Kmer);
			}
		++DataOffset;
		}
	}

uint kappa_dex::GetSeqKmer(const byte *Seq, uint SeqPos, bool SelfScoreMask) const
	{
	uint Kmer = BytesToKmer(Seq + SeqPos);
	if (SelfScoreMask && m_KmerSelfScores[Kmer] < m_MinKmerSelfScore)
		Kmer = UINT_MAX;
	return Kmer;
	}

// Extract exact kmers (skip low self-score). Returns count written to Out[].
// Out must hold at least max(0, L-K+1) slots.
static uint ExtractExactKmers(
	const byte *Seq, uint L,
	const uint8_t *Offsets, uint k, uint K,
	const int16_t *SelfScores, int MinSelf,
	uint *Out)
	{
	if (L < K)
		return 0;
	uint n = 0;
	for (uint pos = 0; pos + K <= L; ++pos)
		{
		if (pos > UINT16_MAX)
			break;
		uint Kmer = 0;
		for (uint i = 0; i < k; ++i)
			Kmer = Kmer*KAPPA_AS + Seq[pos + Offsets[i]];
		if (SelfScores != 0 && SelfScores[Kmer] < MinSelf)
			continue;
		Out[n++] = Kmer;
		}
	return n;
	}

void kappa_dex::from_codeseqs(
	uint8_t **kappa_codeseqs,
	const uint *lengths,
	const vector<string> &labels,
	uint nseq)
	{
	m_kappa_codeseqs = kappa_codeseqs;
	m_seq_lengths = lengths;
	m_nseq = nseq;
	if (m_AddNeighborhood)
		asserta(m_ptrScoreMx != 0);

	if (nseq == 0)
		{
		Alloc_Pass1();
		AdjustFinger();
		Alloc_Pass2();
		SetRowSizes();
		return;
		}

	Alloc_Pass1();

	const uint ThreadCount = GetRequestedThreadCount();
	const bool do_parallel = (!m_AddNeighborhood && ThreadCount > 1 && nseq >= 64);

	if (!do_parallel)
		{
		for (uint SeqIdx = 0; SeqIdx < m_nseq; ++SeqIdx)
			{
			ProgressStep(SeqIdx, m_nseq, "kappa_dex pass 1");
			const char *Label = 0;
			const byte *Seq = kappa_codeseqs[SeqIdx];
			const uint L = lengths[SeqIdx];
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
		for (uint SeqIdx = 0; SeqIdx < m_nseq; ++SeqIdx)
			{
			ProgressStep(SeqIdx, m_nseq, "kappa_dex pass 2");
			const char *Label = 0;
			const byte *Seq = kappa_codeseqs[SeqIdx];
			const uint L = lengths[SeqIdx];
			SetSeq(SeqIdx, Label, Seq, L);
			AddSeq_Pass2();
			}
		SetRowSizes();
#if KAPPA_DEBUG_CHECKS
		CheckAfterPass2();
		Validate();
#endif
		return;
		}

	ProgressLog("kappa_dex parallel build  threads=%u  nseq=%u\n",
		ThreadCount, nseq);

	const uint DictSize = m_DictSize;
	const uint k = m_k;
	const uint K = m_K;
	const uint8_t *Offsets = m_Offsets;
	const int16_t *SelfScores = m_KmerSelfScores;
	const int MinSelf = m_MinKmerSelfScore;

	uint64_t **tls_counts = myalloc(uint64_t *, ThreadCount);
	for (uint t = 0; t < ThreadCount; ++t)
		{
		tls_counts[t] = myalloc(uint64_t, DictSize);
		zero_array(tls_counts[t], DictSize);
		}

	{
	std::atomic<uint> next_seq{0};
	std::mutex progress_lock;
	vector<thread *> ts;
	ProgressStep(0, nseq, "kappa_dex pass 1");
	for (uint tid = 0; tid < ThreadCount; ++tid)
		{
		ts.push_back(new thread(
			[&, tid]()
			{
			uint64_t *counts = tls_counts[tid];
			vector<uint> kmers;
			kmers.reserve(512);
			for (;;)
				{
				const uint SeqIdx = next_seq.fetch_add(1, std::memory_order_relaxed);
				if (SeqIdx >= nseq)
					break;
				if ((SeqIdx & 0x3ff) == 0)
					{
					lock_guard<mutex> lock(progress_lock);
					if (SeqIdx < nseq)
						ProgressStep(SeqIdx, nseq, "kappa_dex pass 1");
					}
				const byte *Seq = kappa_codeseqs[SeqIdx];
				const uint L = lengths[SeqIdx];
				if (L < K)
					continue;
				const uint maxn = L - K + 1;
				if (kmers.size() < maxn)
					kmers.resize(maxn);
				const uint nk = ExtractExactKmers(Seq, L, Offsets, k, K,
					SelfScores, MinSelf, kmers.data());
				for (uint i = 0; i < nk; ++i)
					counts[kmers[i]] += 1;
				}
			}));
		}
	for (uint tid = 0; tid < ThreadCount; ++tid)
		{
		ts[tid]->join();
		delete ts[tid];
		}
	ProgressStep(nseq - 1, nseq, "kappa_dex pass 1");
	}

	m_Size = 0;
	for (uint Kmer = 0; Kmer < DictSize; ++Kmer)
		{
		uint64_t c = 0;
		for (uint t = 0; t < ThreadCount; ++t)
			c += tls_counts[t][Kmer];
		m_Finger[Kmer + 1] = c;
		m_Size += c;
		}
	for (uint t = 0; t < ThreadCount; ++t)
		myfree(tls_counts[t]);
	myfree(tls_counts);

#if KAPPA_DEBUG_CHECKS
	CheckAfterPass1();
#endif
	AdjustFinger();
#if KAPPA_DEBUG_CHECKS
	CheckAfterAdjust();
#endif
	Alloc_Pass2();

	{
	std::atomic<uint> next_seq{0};
	std::mutex progress_lock;
	vector<thread *> ts;
	ProgressStep(0, nseq, "kappa_dex pass 2");
	for (uint tid = 0; tid < ThreadCount; ++tid)
		{
		ts.push_back(new thread(
			[&]()
			{
			for (;;)
				{
				const uint SeqIdx = next_seq.fetch_add(1, std::memory_order_relaxed);
				if (SeqIdx >= nseq)
					break;
				if ((SeqIdx & 0x3ff) == 0)
					{
					lock_guard<mutex> lock(progress_lock);
					if (SeqIdx < nseq)
						ProgressStep(SeqIdx, nseq, "kappa_dex pass 2");
					}
				const byte *Seq = kappa_codeseqs[SeqIdx];
				const uint L = lengths[SeqIdx];
				if (L < K)
					continue;
				for (uint pos = 0; pos + K <= L; ++pos)
					{
					if (pos > UINT16_MAX)
						break;
					uint Kmer = 0;
					for (uint j = 0; j < k; ++j)
						Kmer = Kmer*KAPPA_AS + Seq[pos + Offsets[j]];
					if (SelfScores != 0 && SelfScores[Kmer] < MinSelf)
						continue;
					const uint64_t DataOffset =
						AtomicFetchAddU64(&m_Finger[Kmer + 1], 1);
					Put(DataOffset, SeqIdx, uint16_t(pos));
					}
				}
			}));
		}
	for (uint tid = 0; tid < ThreadCount; ++tid)
		{
		ts[tid]->join();
		delete ts[tid];
		}
	ProgressStep(nseq - 1, nseq, "kappa_dex pass 2");
	}

	SetRowSizes();
#if KAPPA_DEBUG_CHECKS
	CheckAfterPass2();
	Validate();
#endif
	}

void kappa_dex::FromSeqDB(const SeqDB &Input)//TODO FromBags already have Mu k-mers
	{
	//m_SeqDB = &Input;
	const uint SeqCount = Input.GetSeqCount();
	m_labels = &Input.m_Labels;
	uint *seq_lengths = myalloc(uint, SeqCount);
	m_seq_lengths = seq_lengths;
	m_nseq = SeqCount;
	if (m_AddNeighborhood)
		asserta(m_ptrScoreMx != 0);
	//if (m_AddNeighborhood && m_ptrScoreMx == 0)
	//	m_ptrScoreMx = &GetMuMerMx(m_k);

	Alloc_Pass1();
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, SeqCount, "kappa_dex pass 1");
		const char *Label = get_label(SeqIdx);
		const byte *Seq = Input.GetByteSeq(SeqIdx);
		const uint L = Input.GetSeqLength(SeqIdx);
		seq_lengths[SeqIdx] = L;
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
		const char *Label = get_label(SeqIdx);
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
	m_RowSizes = myalloc(uint64_t, m_DictSize);
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		m_RowSizes[Kmer] = m_Finger[Kmer+1] - m_Finger[Kmer];
	}

void kappa_dex::Put(uint64_t DataOffset, uint32_t SeqIdx, uint16_t SeqPos)
	{
	asserta(DataOffset < m_Size);
	asserta(m_Data != 0);
	uint64 Bytes64 = uint64(m_ItemSize)*DataOffset;
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

void kappa_dex::Get(uint64_t DataOffset, uint32_t &SeqIdx, uint16_t &SeqPos) const
	{
	const uint8_t *ptr = m_Data + m_ItemSize*DataOffset;
	SeqIdx = *(uint32_t *) ptr;
	SeqPos = *(uint16_t *) (ptr + 4);
	}

/***
Binary kappa_dex layout (little-endian):
	uint32 MAGIC ('KDEX')
	uint32 VERSION (2)
	uint32 KAPPA_AS (must be 32)
	uint32 ItemSize (must be 6)
	uint32 k
	uint32 K
	uint32 DictSize
	uint32 nseq
	uint64 Size          // posting count (v1 used uint32)
	uint32 MinKmerSelfScore
	uint8  Offsets[k]
	uint64 Finger[DictSize+2]     // v1: uint32
	uint64 RowSizes[DictSize]     // v1: uint32
	uint8  Data[Size*ItemSize]
	uint32 MAGIC
***/
void kappa_dex::ToFile(const string &FN) const
	{
	if (FN == "") return;
	asserta(FN != "");
	asserta(m_Finger != 0);
	asserta(m_RowSizes != 0);
	asserta(m_Offsets != 0);
	asserta(m_Size == 0 || m_Data != 0);
	asserta(m_k > 0 && m_k <= 32);
	asserta(m_K >= m_k);
	asserta(m_DictSize > 0);
	asserta(m_Finger[0] == 0);
	asserta(m_Finger[m_DictSize + 1] == m_Size);

	FILE *f = CreateStdioFile(FN);
	uint32_t Magic = KAPPA_DEX_MAGIC;
	uint32_t Version = KAPPA_DEX_VERSION;
	uint32_t AS = KAPPA_AS;
	uint32_t ItemSize = m_ItemSize;
	uint32_t k = m_k;
	uint32_t K = m_K;
	uint32_t DictSize = m_DictSize;
	uint32_t nseq = m_nseq;
	uint64_t Size = m_Size;
	uint32_t MinSelf = uint32_t(m_MinKmerSelfScore);

	WriteStdioFile(f, &Magic, sizeof(Magic));
	WriteStdioFile(f, &Version, sizeof(Version));
	WriteStdioFile(f, &AS, sizeof(AS));
	WriteStdioFile(f, &ItemSize, sizeof(ItemSize));
	WriteStdioFile(f, &k, sizeof(k));
	WriteStdioFile(f, &K, sizeof(K));
	WriteStdioFile(f, &DictSize, sizeof(DictSize));
	WriteStdioFile(f, &nseq, sizeof(nseq));
	WriteStdioFile(f, &Size, sizeof(Size));
	WriteStdioFile(f, &MinSelf, sizeof(MinSelf));
	WriteStdioFile(f, m_Offsets, k);

	const uint64 FingerBytes = uint64(DictSize + 2) * sizeof(uint64_t);
	const uint64 RowBytes = uint64(DictSize) * sizeof(uint64_t);
	const uint64 DataBytes = Size * uint64(ItemSize);
	asserta(FingerBytes <= UINT_MAX);
	asserta(RowBytes <= UINT_MAX);
	WriteStdioFile(f, m_Finger, uint32(FingerBytes));
	WriteStdioFile(f, m_RowSizes, uint32(RowBytes));
	if (DataBytes > 0)
		WriteStdioFile64(f, m_Data, DataBytes);
	WriteStdioFile(f, &Magic, sizeof(Magic));
	CloseStdioFile(f);

	ProgressLog("Wrote kappa_dex %s  nseq=%u  postings=%s (%s)  dict=%u\n",
		FN.c_str(), nseq, Int64ToStr(Size), MemBytesToStr(double(DataBytes)),
		DictSize);
	}

void kappa_dex::FromFile(const string &FN)
	{
	asserta(FN != "");
	asserta(m_Finger == 0);
	asserta(m_Data == 0);
	asserta(m_RowSizes == 0);

	FILE *f = OpenStdioFile(FN);
	const uint64 FileSize = GetStdioFileSize64(f);
	asserta(FileSize > 40);

	uint32_t Magic = 0;
	uint32_t Version = 0;
	uint32_t AS = 0;
	uint32_t ItemSize = 0;
	uint32_t k = 0;
	uint32_t K = 0;
	uint32_t DictSize = 0;
	uint32_t nseq = 0;
	uint32_t MinSelf = 0;

	ReadStdioFile(f, &Magic, sizeof(Magic));
	asserta(Magic == KAPPA_DEX_MAGIC);
	ReadStdioFile(f, &Version, sizeof(Version));
	asserta(Version == KAPPA_DEX_VERSION || Version == KAPPA_DEX_VERSION_V1);
	ReadStdioFile(f, &AS, sizeof(AS));
	asserta(AS == KAPPA_AS);
	ReadStdioFile(f, &ItemSize, sizeof(ItemSize));
	asserta(ItemSize == m_ItemSize);
	ReadStdioFile(f, &k, sizeof(k));
	ReadStdioFile(f, &K, sizeof(K));
	ReadStdioFile(f, &DictSize, sizeof(DictSize));
	ReadStdioFile(f, &nseq, sizeof(nseq));

	uint64_t Size = 0;
	if (Version == KAPPA_DEX_VERSION)
		ReadStdioFile(f, &Size, sizeof(Size));
	else
		{
		uint32_t Size32 = 0;
		ReadStdioFile(f, &Size32, sizeof(Size32));
		Size = Size32;
		}

	ReadStdioFile(f, &MinSelf, sizeof(MinSelf));
	asserta(k > 0 && k <= 32);
	asserta(K >= k);
	asserta(DictSize > 0);
	asserta(DictSize == myipow(AS, k));

	m_k = k;
	m_K = K;
	m_DictSize = DictSize;
	m_nseq = nseq;
	m_Size = Size;
	m_MinKmerSelfScore = int(MinSelf);
	m_AddNeighborhood = false;
	m_ptrScoreMx = 0;
	m_KmerSelfScores = 0;

	m_Offsets = myalloc(uint8_t, k);
	ReadStdioFile(f, m_Offsets, k);
	for (uint i = 0; i < k; ++i)
		asserta(m_Offsets[i] < K);

	const uint64 DataBytes = Size * uint64(ItemSize);
	m_Finger = myalloc(uint64_t, DictSize + 2);
	m_RowSizes = myalloc(uint64_t, DictSize);

	if (Version == KAPPA_DEX_VERSION)
		{
		const uint64 FingerBytes = uint64(DictSize + 2) * sizeof(uint64_t);
		const uint64 RowBytes = uint64(DictSize) * sizeof(uint64_t);
		asserta(FingerBytes <= UINT_MAX);
		asserta(RowBytes <= UINT_MAX);
		ReadStdioFile(f, m_Finger, uint32(FingerBytes));
		ReadStdioFile(f, m_RowSizes, uint32(RowBytes));
		}
	else
		{
		const uint64 FingerBytes = uint64(DictSize + 2) * sizeof(uint32_t);
		const uint64 RowBytes = uint64(DictSize) * sizeof(uint32_t);
		asserta(FingerBytes <= UINT_MAX);
		asserta(RowBytes <= UINT_MAX);
		uint32_t *Finger32 = myalloc(uint32_t, DictSize + 2);
		uint32_t *RowSizes32 = myalloc(uint32_t, DictSize);
		ReadStdioFile(f, Finger32, uint32(FingerBytes));
		ReadStdioFile(f, RowSizes32, uint32(RowBytes));
		for (uint i = 0; i < DictSize + 2; ++i)
			m_Finger[i] = Finger32[i];
		for (uint i = 0; i < DictSize; ++i)
			m_RowSizes[i] = RowSizes32[i];
		myfree(Finger32);
		myfree(RowSizes32);
		}

	asserta(m_Finger[0] == 0);
	asserta(m_Finger[DictSize + 1] == Size);

	if (DataBytes > 0)
		{
		m_Data = myalloc64(uint8_t, DataBytes);
		ReadStdioFile64NoPos(f, m_Data, DataBytes);
		}
	else
		m_Data = 0;

	for (uint Kmer = 0; Kmer < DictSize; ++Kmer)
		{
		asserta(m_Finger[Kmer] <= m_Finger[Kmer + 1]);
		asserta(m_Finger[Kmer + 1] - m_Finger[Kmer] == m_RowSizes[Kmer]);
		}

	uint32_t Magic2 = 0;
	ReadStdioFile(f, &Magic2, sizeof(Magic2));
	asserta(Magic2 == KAPPA_DEX_MAGIC);
	asserta(GetStdioFilePos64(f) == FileSize);
	CloseStdioFile(f);

	ProgressLog("Read kappa_dex %s  nseq=%u  postings=%s (%s)  dict=%u\n",
		FN.c_str(), nseq, Int64ToStr(Size), MemBytesToStr(double(DataBytes)),
		DictSize);
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
			uint64_t Size = GetRowSize(Kmer);
			Kmers.push_back(Kmer);
			Sizes.push_back(Size > UINT_MAX ? UINT_MAX : uint(Size));
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
		assert(Kmer < flat_params::m_kappa_dict_size);
		if (m_KmerSelfScores != 0 && m_KmerSelfScores[Kmer] < m_MinKmerSelfScore)
			Kmers.push_back(UINT_MAX);
		else
			Kmers.push_back(Kmer%m_DictSize);
		}
	}
