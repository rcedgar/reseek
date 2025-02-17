#include "myutils.h"
#include "alpha.h"
#include "threedex.h"
#include "seqdb.h"
#include "quarts.h"

#define TRACE	0

/***
32 bits 2^5, 64 bits 2^6
5 bits per aa letter
6 bits per Mu letter
8 bits per char
***/

const char ThreeDex::m_Pattern[11] = "1101010011";
const byte ThreeDex::m_Offsets[6] = { 0, 1, 3, 5, 8, 9 };
const uint32_t ThreeDex::m_k = 6;	// number of 1s in pattern
const uint32_t ThreeDex::m_K = 10;	// pattern length
const uint32_t ThreeDex::m_DictSize = 20*20*20*20*20*20;
const uint32_t ThreeDex::m_ItemSize = 6;	// 4 byte SeqIdx + 2 byte Pos

static bool ValidatePattern()
	{
	asserta(strlen(ThreeDex::m_Pattern) == ThreeDex::m_K);
	uint ones = 0;
	uint k = 0;
	for (uint i = 0; i < ThreeDex::m_K; ++i)
		if (ThreeDex::m_Pattern[i] == '1')
			{
			++ones;
			asserta(ThreeDex::m_Offsets[k++] == i);
			}
	asserta(ThreeDex::m_k == ones);
	asserta(ThreeDex::m_DictSize == myipow(20, 6));
	return true;
	}

static bool s_ValidatePatternDone = ValidatePattern();

uint ThreeDex::StrToKmer(const string &s) const
	{
	assert(SIZE(s) == m_k);
	return StrToKmer(s.c_str());
	}

uint ThreeDex::StrToKmer(const char *s) const
	{
	uint Kmer = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		Kmer *= 20;
		byte c = s[i];
		uint Letter = g_CharToLetterAmino[c];
		asserta(Letter < 20);
		Kmer += Letter;
		}
	return Kmer;
	}

uint ThreeDex::BytesToKmer(const byte *s) const
	{
	uint Kmer = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		Kmer *= 20;
		uint Letter = s[i];
		asserta(Letter < 20);
		Kmer += Letter;
		}
	return Kmer;
	}

const char *ThreeDex::KmerToStr(uint Kmer, string &s) const
	{
	if (Kmer == UINT_MAX)
		{
		s.clear();
		for (uint i = 0; i < m_K; ++i)
			s += '*';
		return s.c_str();
		}

	s.clear();
	for (uint i = 0; i < m_K; ++i)
		{
		if (m_Pattern[m_K-i-1] == '0')
			{
			s.push_back('_');
			continue;
			}
		byte Letter = Kmer%20;
		s.push_back(g_LetterToCharAmino[Letter]);
		Kmer /= 20;
		}
	reverse(s.begin(), s.end());
	return s.c_str();
	}

void ThreeDex::Alloc_Pass1()
	{
// Pass1 m_Finger[Kmer] = Count
	asserta(m_Finger == 0 && m_Data == 0);
	m_Finger = myalloc(uint32_t, m_DictSize + 2);
	zero_array(m_Finger, m_DictSize+2);
#if DEBUG
	m_KmerToCount1.resize(m_DictSize, 0);
#endif
	}

void ThreeDex::Alloc_Pass2()
	{
// 6 bytes for uint32_t:uint16_t (SeqIdx:Pos)
	const uint Bytes = m_ItemSize*m_Size;
	m_Data = myalloc(uint8_t, Bytes);
#if DEBUG
	m_KmerToCount2.resize(m_DictSize, 0);
	memset(m_Data, 0xff, Bytes);
#endif
	}

void ThreeDex::SetSeq(uint SeqIdx, const char *Label, const byte *Seq, uint L)
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

void ThreeDex::AddSeq_Pass1()
	{
#if TRACE
	string Tmp;
	Log("AddSeq_Pass1(%s) L=%u\n", m_Label, m_L);
#endif
	const uint KmerCount = SIZE(m_Kmers);
	for (uint i = 0; i < KmerCount; ++i)
		{
		uint Kmer = m_Kmers[i];
		if (Kmer == UINT_MAX)
			{
#if TRACE
			Log("[%4u] ***\n", i);
#endif
			continue;
			}

	// Pass 1, m_Finger[Kmer+1] is count
		m_Finger[Kmer+1] += 1;
		++m_Size;
#if DEBUG
		m_KmerToCount1[Kmer] += 1;
#endif
#if TRACE
		Log("[%4u] %08x %s", i, Kmer, KmerToStr(Kmer, Tmp));
		if (m_KmerSelfScores != 0)
			Log(" self=%d", m_KmerSelfScores[Kmer]);
		Log("\n");
#endif
		}
	}

void ThreeDex::AddSeq_Pass2()
	{
#if TRACE
	string Tmp;
	Log("AddSeq_Pass2(%s) L=%u\n", m_Label, m_L);
#endif
	const uint KmerCount = SIZE(m_Kmers);
	for (uint i = 0; i < KmerCount; ++i)
		{
		uint Kmer = m_Kmers[i];
		if (Kmer == UINT_MAX)
			{
#if TRACE
			Log("[%4u] %08x %s --LOW\n", i, Kmer, KmerToStr(Kmer, Tmp));
#endif
			continue;
			}
		uint DataOffset = m_Finger[Kmer+1];
		Put(DataOffset, m_SeqIdx, i);
		m_Finger[Kmer+1] += 1;
#if DEBUG
		assert(m_KmerToDataStart[Kmer] + m_KmerToCount2[Kmer] == DataOffset);
		m_KmerToCount2[Kmer] += 1;
#endif
#if TRACE
		Log("[%4u] %08x %s DO=%u\n",
			i, Kmer, KmerToStr(Kmer, Tmp), DataOffset);
#endif
		}
	}

void ThreeDex::LogStats() const
	{
	vector<uint> RowSizes;
	uint Sum = 0;
	uint MaxSize = 0;
	uint MaxKmer = UINT_MAX;
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint Size = GetRowSize(Kmer);
		if (Size > MaxSize)
			{
			MaxSize = Size;
			MaxKmer = Kmer;
			}
		Sum += Size;
		RowSizes.push_back(Size);
		}
	Quarts Q;
	GetQuarts(RowSizes, Q);
	Log("RowSizes: ");
	Q.LogMe();
	Log("Total = %u (%s)\n", Sum, IntToStr(Sum));
	string s;
	Log("Max kmer %s", KmerToStr(MaxKmer, s));
	if (m_KmerSelfScores != 0)
		Log(" (self=%d)", m_KmerSelfScores[MaxKmer]);
	Log("\n");
	}

#if DEBUG
void ThreeDex::CheckAfterPass1() const
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

void ThreeDex::CheckAfterAdjust() const
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

void ThreeDex::CheckAfterPass2() const
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

void ThreeDex::AdjustFinger()
	{
	uint Sum = 0;
	for (uint Kmer = 0; Kmer <= m_DictSize; ++Kmer)
		{
#if DEBUG
		m_KmerToDataStart.push_back(Sum);
#endif
		uint Kmer_Size = m_Finger[Kmer+1];
		m_Finger[Kmer+1] = Sum;
		Sum += Kmer_Size;
		}
	asserta(Sum == m_Size);
	}

void ThreeDex::Validate() const
	{
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		ValidateKmer(Kmer);
	}

void ThreeDex::LogSeq() const
	{
	Log(">%s(%u)\n", m_Label, m_L);
	for (uint i = 0; i < m_L; ++i)
		{
		if (i > 0 && i%80 == 0)
			Log("\n");
		byte c = g_LetterToCharAmino[m_Seq[i]];
		Log("%c", c);
		}
	Log("\n");
	}

void ThreeDex::LogMe() const
	{
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		if (GetRowSize(Kmer) == 0)
			continue;
		LogIndexKmer(Kmer);
		}
	}

void ThreeDex::LogIndexKmer(uint Kmer) const
	{
	uint n = GetRowSize(Kmer);
	string Tmp;
	uint DataOffset = m_Finger[Kmer];
	int Self = -1;
	if (m_KmerSelfScores != 0)
		Self = m_KmerSelfScores[Kmer];
	Log("LogIndexKmer(%08x) %s size=%u self=%5d | ",
		Kmer, KmerToStr(Kmer, Tmp), n, Self);
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

void ThreeDex::ValidateKmer(uint Kmer) const
	{
	uint n = GetRowSize(Kmer);
	uint DataOffset = m_Finger[Kmer];
	for (uint i = 0; i < n; ++i)
		{
		uint32_t SeqIdx;
		uint16_t SeqPos;
		Get(DataOffset, SeqIdx, SeqPos);
		const byte *Seq = m_SeqDB->GetByteSeq(SeqIdx);
		uint Check_Kmer = GetSeqKmer(Seq, SeqPos, false);
		if (Check_Kmer != Kmer)
			{
			string Tmp;
			Log("ValidateKmer(%u=%s)\n", Kmer, KmerToStr(Kmer, Tmp));
			Log("GetSeqKmer(SeqIdx=%u, SeqPos=%u) %u=%s\n",
				SeqIdx, SeqPos, Check_Kmer, KmerToStr(Check_Kmer, Tmp));
			Die("ValidateKmer");
			}
		}
	}

uint ThreeDex::GetSeqKmer(const byte *Seq, uint SeqPos, bool SelfScoreMask) const
	{
	uint Kmer = SeqBytesToKmer(Seq, SeqPos, SelfScoreMask);
	return Kmer;
	}

void ThreeDex::FromSeqDB(const SeqDB &Input)
	{
	m_SeqDB = &Input;
	const uint SeqCount = Input.GetSeqCount();

	Alloc_Pass1();
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, SeqCount, "ThreeDex pass 1");
		const char *Label = m_SeqDB->GetLabel(SeqIdx).c_str();
		const byte *Seq = Input.GetByteSeq(SeqIdx);
		const uint L = Input.GetSeqLength(SeqIdx);
		SetSeq(SeqIdx, Label, Seq, L);
		AddSeq_Pass1();
		}
#if DEBUG
	CheckAfterPass1();
#endif

	AdjustFinger();
#if DEBUG
	CheckAfterAdjust();
#endif

	Alloc_Pass2();
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, SeqCount, "ThreeDex pass 2");
		const char *Label = m_SeqDB->GetLabel(SeqIdx).c_str();
		const byte *Seq = Input.GetByteSeq(SeqIdx);
		const uint L = Input.GetSeqLength(SeqIdx);
		SetSeq(SeqIdx, Label, Seq, L);
		AddSeq_Pass2();
		}
#if DEBUG
	CheckAfterPass2();
#endif

#if DEBUG
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
#endif
	}

uint ThreeDex::SeqBytesToKmer(const byte *Seq, uint Pos, bool SelfScoreMask) const
	{
	uint Kmer = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		byte Letter = Seq[Pos + m_Offsets[i]];
		Kmer = Kmer*20 + Letter;
		}
	if (SelfScoreMask && m_KmerSelfScores[Kmer] < m_MinKmerSelfScore)
		Kmer = UINT_MAX;
	return Kmer;
	}

void ThreeDex::GetKmers(const byte *Seq, uint L, vector<uint> &Kmers) const
	{
	Kmers.reserve(L);
	Kmers.clear();
	for (uint KmerStartPos = 0; KmerStartPos + m_K <= L; ++KmerStartPos)
		{
		uint Kmer = 0;
		for (uint i = 0; i < m_k; ++i)
			{
			byte Letter = Seq[KmerStartPos + m_Offsets[i]];
			Kmer = Kmer*20 + Letter;
			}
		asserta(Kmer == GetSeqKmer(Seq, KmerStartPos, false));
		if (m_KmerSelfScores != 0 && m_KmerSelfScores[Kmer] < m_MinKmerSelfScore)
			Kmers.push_back(UINT_MAX);
		else
			Kmers.push_back(Kmer%m_DictSize);
		}
	}

void cmd_threedex()
	{
	SeqDB Input;
	Input.FromFasta(g_Arg1);

	ThreeDex TD;
	TD.FromSeqDB(Input);
	//TD.LogIndexKmer(0);
	TD.LogIndexKmer(4800);
	TD.LogStats();
	TD.Validate();
	ProgressLog("Validate OK\n");
	}
