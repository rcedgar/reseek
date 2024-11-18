#include "myutils.h"
#include "alpha.h"
#include "mudex.h"
#include "seqdb.h"
#include "quarts.h"

uint MuDex::StrToKmer(const string &s) const
	{
	asserta(SIZE(s) == 5);
	uint Kmer = 0;
	for (uint i = 0; i < 5; ++i)
		{
		Kmer *= 36;
		byte c = s[i];
		uint Letter = g_CharToLetterMu[c];
		asserta(Letter < 36);
		Kmer += Letter;
		}
	string Tmp;
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
	asserta(m_RowSizes == 0);
	m_SeqIdxRows = myalloc(uint32_t *, m_DictSize);
	m_SeqPosRows = myalloc(uint16_t *, m_DictSize);
	m_RowSizes = myalloc(uint32_t, m_DictSize);
	zero_array(m_RowSizes, m_DictSize);
	}

void MuDex::Alloc_Pass2()
	{
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint32_t Size = m_RowSizes[Kmer];
		if (Size == 0)
			continue;
		m_SeqIdxRows[Kmer] = myalloc(uint32_t, Size);
		m_SeqPosRows[Kmer] = myalloc(uint16_t, Size);
		}
	m_RowSizes2 = myalloc(uint32_t, m_DictSize);
	zero_array(m_RowSizes2, m_DictSize);
	}

void MuDex::AddSeq_Pass1(uint SeqIdx, const char *Label, const char *Seq, uint L)
	{
#if 0
	SeqToFasta(g_fLog, Label, Seq, L);
#endif
	if (L < 4)
		return;
	uint Kmer = 0;
	for (uint Pos = 0; Pos < 4; ++Pos)
		{
		byte Letter = g_CharToLetterMu[Seq[Pos]];
		assert(Letter < 36);
		Kmer = Kmer*36 + Letter;
		}

#if 0
	string Tmp;
#endif
	for (uint Pos = 4; Pos < L; ++Pos)
		{
		byte Letter = g_CharToLetterMu[Seq[Pos]];
		assert(Letter < 36);
		Kmer = Kmer*36 + Letter;
		Kmer %= m_DictSize;
		assert(Kmer < m_DictSize);
		m_RowSizes[Kmer] += 1;
#if 0
		Log("[%4u] %s\n", Pos-4, KmerToStr(Kmer, Tmp));
#endif
		}
	}

void MuDex::AddSeq_Pass2(uint SeqIdx, const char *Label, const char *Seq, uint L)
	{
#if 0
	SeqToFasta(g_fLog, Label, Seq, L);
#endif
	if (L < 4)
		return;
	uint Kmer = 0;
	for (uint Pos = 0; Pos < 4; ++Pos)
		{
		byte Letter = g_CharToLetterMu[Seq[Pos]];
		assert(Letter < 36);
		Kmer = Kmer*36 + Letter;
		}

#if 0
	string Tmp;
#endif
	for (uint Pos = 4; Pos < L; ++Pos)
		{
		byte Letter = g_CharToLetterMu[Seq[Pos]];
		assert(Letter < 36);
		Kmer = Kmer*36 + Letter;
		Kmer %= m_DictSize;
		assert(Kmer < m_DictSize);
		uint n = m_RowSizes2[Kmer];
		m_SeqIdxRows[Kmer][n] = SeqIdx;
		m_SeqPosRows[Kmer][n] = Pos-5;
		m_RowSizes2[Kmer] = n;
#if 0
		Log("[%4u] %s\n", Pos-4, KmerToStr(Kmer, Tmp));
#endif
		}
	}

void MuDex::LogStats() const
	{
	vector<uint> RowSizes;
	uint Sum = 0;
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		uint Size = m_RowSizes[Kmer];
		Sum += Size;
		RowSizes.push_back(Size);
		}
	Quarts Q;
	GetQuarts(RowSizes, Q);
	Log("RowSizes: ");
	Q.LogMe();
	Log("Total = %u (%s)\n", Sum, IntToStr(Sum));
	}

void cmd_mudex()
	{
	MuDex MD;
	MD.Alloc_Pass1();
	SeqDB Input;
	Input.FromFasta(g_Arg1);
	const uint SeqCount = Input.GetSeqCount();
	uint SumL = 0;
	uint SumL4 = 0;
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		ProgressStep(SeqIdx, SeqCount, "Pass 1");
		const string &Seq = Input.GetSeq(SeqIdx);
		const uint L = SIZE(Seq);
		MD.AddSeq_Pass1(SeqIdx, Input.GetLabel(SeqIdx).c_str(), Seq.c_str(), L);
		SumL += L;
		if (L >= 5)
			SumL4 += L - 4;
		}
	MD.LogStats();
	Log("SumL = %u %u\n", SumL, SumL4);
	}
