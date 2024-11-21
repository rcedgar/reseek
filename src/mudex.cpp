#include "myutils.h"
#include "alpha.h"
#include "mudex.h"
#include "seqdb.h"
#include "quarts.h"

#define TRACE	0

uint MuDex::StrToKmer(const string &s) const
	{
	assert(SIZE(s) == m_k);
	return StrToKmer(s.c_str());
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

const char *MuDex::KmerToStr(uint Kmer, string &s) const
	{
	s.clear();
	for (uint i = 0; i < 5; ++i)
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
// Pass1 m_Finger[Kmer] = Count
	asserta(m_Finger == 0 && m_Data == 0);
	m_Finger = myalloc(uint32_t, m_DictSize + 2);
	zero_array(m_Finger, m_DictSize+2);
#if DEBUG
	m_KmerToCount1.resize(m_DictSize, 0);
#endif
	}

void MuDex::Alloc_Pass2()
	{
// 6 bytes for uint32_t:uint16_t (SeqIdx:Pos)
	const uint Bytes = m_ItemSize*m_Size;
	m_Data = myalloc(uint8_t, Bytes);
#if DEBUG
	m_KmerToCount2.resize(m_DictSize, 0);
	memset(m_Data, 0xff, Bytes);
#endif
	}

void MuDex::AddSeq_Pass1(uint SeqIdx, const char *Label, const char *Seq, uint L)
	{
#if TRACE
	Log("AddSeq_Pass1(%s) L=%u\n", Label, L);
	SeqToFasta(g_fLog, Label, Seq, L);
	string Tmp;
#endif
	if (L < m_k)
		return;
	uint Kmer = 0;
	for (uint SeqPos = 0; SeqPos < m_k-1; ++SeqPos)
		{
		byte Letter = g_CharToLetterMu[Seq[SeqPos]];
		assert(Letter < 36);
		Kmer = Kmer*36 + Letter;
		}

	for (uint SeqPos = m_k-1; SeqPos < L; ++SeqPos)
		{
		byte Letter = g_CharToLetterMu[Seq[SeqPos]];
		assert(Letter < 36);
		Kmer = Kmer*36 + Letter;
		Kmer %= m_DictSize;
		assert(Kmer < m_DictSize);

	// Pass 1, m_Finger[Kmer+1] is count
		m_Finger[Kmer+1] += 1;
#if DEBUG
		m_KmerToCount1[Kmer] += 1;
#endif
#if TRACE
		Log("[%4u] %08x %s\n", SeqPos-4, Kmer, KmerToStr(Kmer, Tmp));
#endif
		}
	m_Size += L - (m_k-1);
	}

void MuDex::AddSeq_Pass2(uint SeqIdx, const char *Label, const char *Seq, uint L)
	{
#if TRACE
	Log("AddSeq_Pass2(%s) L=%u\n", Label, L);
	SeqToFasta(g_fLog, Label, Seq, L);
	string Tmp;
#endif
	if (L < m_k)
		return;
	uint Kmer = 0;
	for (uint SeqPos = 0; SeqPos < m_k-1; ++SeqPos)
		{
		byte Letter = g_CharToLetterMu[Seq[SeqPos]];
		assert(Letter < 36);
		Kmer = Kmer*36 + Letter;
		}

	for (uint SeqPos = m_k-1; SeqPos < L; ++SeqPos)
		{
		byte Letter = g_CharToLetterMu[Seq[SeqPos]];
		assert(Letter < 36);
		Kmer = Kmer*36 + Letter;
		Kmer %= m_DictSize;
		assert(Kmer < m_DictSize);
		uint DataOffset = m_Finger[Kmer+1];
		uint KmerStartPos = SeqPos - (m_k-1);
		Put(DataOffset, SeqIdx, KmerStartPos);
		m_Finger[Kmer+1] += 1;
#if DEBUG
		assert(m_KmerToDataStart[Kmer] + m_KmerToCount2[Kmer] == DataOffset);
		m_KmerToCount2[Kmer] += 1;
#endif
#if TRACE
		Log("[%4u] %08x %s DO=%u\n",
			KmerStartPos, Kmer, KmerToStr(Kmer, Tmp), DataOffset);
#endif
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

#if DEBUG
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
	assert(m_Finger[m_DictSize] == m_Size);
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
#if DEBUG
		m_KmerToDataStart.push_back(Sum);
#endif
		uint Kmer_Size = m_Finger[Kmer+1];
		m_Finger[Kmer+1] = Sum;
		Sum += Kmer_Size;
		}
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
	uint n = GetRowSize(Kmer);
	uint DataOffset = m_Finger[Kmer];
	for (uint i = 0; i < n; ++i)
		{
		uint32_t SeqIdx;
		uint16_t SeqPos;
		Get(DataOffset, SeqIdx, SeqPos);
		uint Check_Kmer = GetSeqKmer(SeqIdx, SeqPos);
		asserta(Check_Kmer == Kmer);
		}
	}

uint MuDex::GetSeqKmer(uint SeqIdx, uint SeqPos) const
	{
	const string &strSeq = m_SeqDB->GetSeq(SeqIdx);
#if DEBUG
	{
	const uint L = SIZE(strSeq);
	asserta(SeqPos + m_k <= L);
	}
#endif
	uint Kmer = StrToKmer(strSeq.c_str() + SeqPos);
	return Kmer;
	}

void MuDex::FromSeqDB(const SeqDB &Input)
	{
	m_SeqDB = &Input;
	const uint SeqCount = Input.GetSeqCount();

	Alloc_Pass1();
	uint SumL4 = 0;
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, SeqCount, "MuDex pass 1");
		const string &Seq = Input.GetSeq(SeqIdx);
		const uint L = SIZE(Seq);
		AddSeq_Pass1(SeqIdx, Input.GetLabel(SeqIdx).c_str(), Seq.c_str(), L);
		if (L >= 5)
			SumL4 += L - 4;
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
		ProgressStep(SeqIdx, SeqCount, "MuDex pass 2");
		const string &Seq = Input.GetSeq(SeqIdx);
		const uint L = SIZE(Seq);
		AddSeq_Pass2(SeqIdx, Input.GetLabel(SeqIdx).c_str(), Seq.c_str(), L);
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

void cmd_mudex()
	{
	SeqDB Input;
	MuDex MD;
	MD.FromSeqDB(Input);
	MD.LogIndexKmer(0);
	MD.LogIndexKmer(0x01733125);
	MD.LogIndexKmer(0x01bdf141);
	MD.LogIndexKmer(0x01710931);
	MD.LogStats();
	MD.Validate();
	ProgressLog("Validate OK\n");
	}
