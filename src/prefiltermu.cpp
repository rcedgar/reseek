#include "myutils.h"
#include "prefiltermu.h"
#include "sort.h"

mutex PrefilterMu::m_Lock;
RankedScoresBag PrefilterMu::m_RSB;

//////////////////////////////////////////////
// 	FindHSP searches for the highest-scoring
// 	ungapped alignment on a given diagonal.
//////////////////////////////////////////////
int PrefilterMu::FindHSP(uint QSeqIdx, int Diag) const
	{
	const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	const uint QL = m_QDB->GetSeqLength(QSeqIdx);

	asserta(Diag >= 0);
#if 0 // TRACE
	if (DoTrace(QSeqIdx)) Log("FindHSP(QL=%u, Diag=%d)\n", QL, Diag);
#endif

	diag dg(QL, m_TL);
	int i = dg.getmini(Diag);
	int j = dg.getminj(Diag);
	int n = dg.getlen(Diag);

	int B = 0;
	int F = 0;
	for (int k = 0; k < n; ++k)
		{
		assert(i < int(QL));
		assert(j < int(m_TL));
		byte q = QSeq[i++];
		byte t = m_TSeq[j++];
		assert(q < ALPHABET_SIZE);
		assert(t < ALPHABET_SIZE);
		short Score = Mu_S_ij_i8[q][t];
		F += Score;
#if 0 // TRACE
		if (DoTrace(QSeqIdx)) Log(" i=%u j=%u F=%d B=%d score=%d\n", i, j, F, B, Score);
#endif
		if (F > B)
			B = F;
		else if (F < 0)
			F = 0;
		}
	return B;
	}

//////////////////////////////////////////////
// 	FindHSP plus "traceback", i.e. returns
// 	start position and length of HSP.
//////////////////////////////////////////////
int PrefilterMu::FindHSP2(uint QSeqIdx,
						int Diag, int &Lo, int &Len) const
	{
	const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	const uint QL = m_QDB->GetSeqLength(QSeqIdx);

	diag dg(QL, m_TL);
	int i = dg.getmini(Diag);
	int j = dg.getminj(Diag);
	int n = dg.getlen(Diag);

	int B = 0;
	int F = 0;
	int CurrLen = 0;
	Lo = 0;
	Len = 0;
	int SuffixLo = 0;
	for (int k = 0; k < n; ++k)
		{
		assert(i < int(QL));
		assert(j < int(m_TL));
		byte q = QSeq[i++];
		byte t = m_TSeq[j++];
		assert(q < ALPHABET_SIZE);
		assert(t < ALPHABET_SIZE);
		short Score = Mu_S_ij_i8[q][t];
		F += Score;
		if (F > B)
			{
			B = F;
			Lo = SuffixLo;
			Len = ++CurrLen;
			}
		else if (F > 0)
			++CurrLen;
		else
			{
			F = 0;
			SuffixLo = k+1;
			CurrLen = 0;
			}
		}
	return B;
	}

void PrefilterMu::SetQDB(const SeqDB &QDB)
	{
	m_QDB = &QDB;
	m_QSeqCount = QDB.GetSeqCount();

	m_QSeqIdxToBestDiagScore = myalloc(uint16_t, m_QSeqCount);
	m_QSeqIdxsWithTwoHitDiag = myalloc(uint16_t, m_QSeqCount);

	for (uint i = 0; i < m_QSeqCount; ++i)
		{
		m_QSeqIdxToBestDiagScore[i] = 0;
		m_QSeqIdxsWithTwoHitDiag[i] = UINT16_MAX;
		}

	bool TargetNeighborhood = !g_QueryNeighborhood;
	if (TargetNeighborhood)
		m_NeighborKmers = myalloc(uint, MAX_HOOD_SIZE);
	else
		m_NeighborKmers = 0;
	m_NrQueriesWithTwoHitDiag = 0;
	}

void PrefilterMu::Search_TargetKmers()
	{
#if KMER_SORT
	m_QKmerIndex->GetKmersAndSizes(m_TSeq, m_TL, m_TKmers, m_TKmerSizes);
	const uint NK = SIZE(m_TKmers);
	m_TKmerSizeOrder.resize(NK);
	QuickSortOrder(m_TKmerSizes.data(), NK, m_TKmerSizeOrder.data());
	const uint QueryCount = m_QDB->GetSeqCount();
	uint MaxTotalSize = QUERY_COUNT_MULTIPLIER*QueryCount;
	uint SumSize = 0;
	for (uint k = 0; k < NK; ++k)
		{
		uint i = m_TKmerSizeOrder[k];
		SumSize += m_TKmerSizes[i];
		if (SumSize > MaxTotalSize)
			m_TKmers[i] = UINT_MAX;
		}
#else
	m_QKmerIndex->GetKmers(m_TSeq, m_TL, m_TKmers);
	const uint NK = SIZE(m_TKmers);
#endif
	if (g_QueryNeighborhood)
		{
		for (uint TPos = 0; TPos < NK; ++TPos)
			{
			uint TKmer = m_TKmers[TPos];
			if (TKmer != UINT_MAX)
				{
#if TRACE
				m_TBaseKmer = TKmer;
#endif
				Search_TargetKmer(TKmer, TPos);
				}
			}
		}
	else
		{
		for (uint TPos = 0; TPos < NK; ++TPos)
			{
			uint Kmer = m_TKmers[TPos];
			if (Kmer != UINT_MAX)
				Search_TargetKmerNeighborhood(Kmer, TPos);
			}
		}
	}

void PrefilterMu::Search_TargetSeq(uint TSeqIdx, const string &TLabel,
				   const byte *TSeq, uint TL)
	{
	m_TSeqIdx = TSeqIdx;
	m_TLabel = TLabel;
	m_TSeq = TSeq;
	m_TL = TL;

	Reset();
	Search_TargetKmers();
	FindTwoHitDiags();
	ExtendTwoHitDiagsToHSPs();
	}

void PrefilterMu::Search_TargetKmerNeighborhood(uint Kmer, uint TPos)
	{
	if (Kmer == UINT_MAX)
		return;
#if TRACE
	m_TBaseKmer = Kmer;
#endif
	assert(Kmer < DICT_SIZE);
	assert(m_KmerSelfScores[Kmer] >= MIN_KMER_PAIR_SCORE);
	short MinKmerScore = MIN_KMER_PAIR_SCORE;

// Construct high-scoring neighborhood
	const uint HSKmerCount =
		m_ScoreMx->GetHighScoringKmers(Kmer, MinKmerScore, m_NeighborKmers);

#if TRACE
	string Tmp;
	Log("Search_TargetKmerNeighborhood TPos=%u Kmer=%s minscore=%d |nbrs|=%d\n",
		TPos, KmerToStr(Kmer, Tmp), MinKmerScore, HSKmerCount);
#endif
	for (uint HSKmerIdx = 0; HSKmerIdx < HSKmerCount; ++HSKmerIdx)
		{
		uint HSKmer = m_NeighborKmers[HSKmerIdx];
		Search_TargetKmer(HSKmer, TPos);
		}
	}

void PrefilterMu::Search_TargetKmer(uint TKmer, uint TPos)
	{
	uint RowSize = m_QKmerIndex->GetRowSize(TKmer);
#if TRACE
	{
	string KmerStr;
	m_QKmerIndex->KmerToStr(m_TBaseKmer, KmerStr);
	Log("Search_TargetKmer(TPos=%u, TKmer=%s) RowSize=%u\n",
		TPos, KmerStr.c_str(), RowSize);
	}
#endif
	if (RowSize == 0)
		return;
	uint DataOffset = m_QKmerIndex->GetRowStart(TKmer);
	for (uint ColIdx = 0; ColIdx < RowSize; ++ColIdx)
		{
		uint32_t QSeqIdx;
		uint16_t QSeqPos;
		m_QKmerIndex->Get(DataOffset++, QSeqIdx, QSeqPos);
		asserta(QSeqIdx < m_QSeqCount);
		uint QL32 = m_QDB->GetSeqLength(QSeqIdx);
		asserta(QL32 < UINT16_MAX);
		uint16_t QL = uint16_t(QL32);
		diag dg(QL, m_TL);
		uint16_t Diag = dg.getd(QSeqPos, TPos);
#if TRACE
		{
		string TKmerStr;
		string QKmerStr;
		uint QKmer = GetQKmer(QSeqIdx, QSeqPos);
		m_QKmerIndex->KmerToStr(m_TBaseKmer, TKmerStr);
		m_QKmerIndex->KmerToStr(QKmer, QKmerStr);
		const MerMx &MM = GetMuMerMx(k);
		int KmerPairScore = MM.GetScoreKmerPair(TKmer, QKmer);

		Log("@K@  [%4u] %5s  [%4u] %5s  /%5u/  %+3d\n",
			QSeqPos, QKmerStr.c_str(),
			TPos, TKmerStr.c_str(),
			Diag, KmerPairScore);
		}
#endif
		if (Diag > m_Mask14)
			continue;
		m_DiagBag.Add(QSeqIdx, Diag);
		}
	}
	
void PrefilterMu::FindTwoHitDiags()
	{
	m_DiagBag.ClearDupes();
	m_DiagBag.SetDupes();
#if DEBUG
	m_DiagBag.Validate(m_QSeqCount, INT16_MAX);
#endif
	}

void PrefilterMu::GetResults(vector<uint> &QSeqIdxs,
						   vector<uint16_t> &DiagScores) const
	{
	QSeqIdxs.clear();
	DiagScores.clear();
	DiagScores.reserve(m_NrQueriesWithTwoHitDiag);
 	for (uint i = 0; i < m_NrQueriesWithTwoHitDiag; ++i)
		{
		uint QSeqIdx = m_QSeqIdxsWithTwoHitDiag[i];
		uint16_t DiagScore = m_QSeqIdxToBestDiagScore[QSeqIdx];

		QSeqIdxs.push_back(QSeqIdx);
		DiagScores.push_back(DiagScore);
		}
	}

void PrefilterMu::AddTwoHitDiag(uint QSeqIdx, uint16_t Diag, int DiagScore)
	{
	if (DiagScore <= 0)
		return;
	asserta(QSeqIdx < UINT16_MAX);
	if (DiagScore >= UINT16_MAX)
		DiagScore = UINT16_MAX-1;
	uint16_t BestDiagScoreT = m_QSeqIdxToBestDiagScore[QSeqIdx];
	if (BestDiagScoreT == 0)
		{
		m_QSeqIdxsWithTwoHitDiag[m_NrQueriesWithTwoHitDiag++] = QSeqIdx;
		m_QSeqIdxToBestDiagScore[QSeqIdx] = DiagScore;
#if TRACE
		if (DoTrace(QSeqIdx)) Log("AddTwoHitDiag(TSeqIdx=%u, QSeqIdx=%u, Diag=%u, DiagScore=%d) (first)\n",
								  m_TSeqIdx, QSeqIdx, Diag, DiagScore);
#endif
		}
	else if (DiagScore > m_QSeqIdxToBestDiagScore[QSeqIdx])
		{
		m_QSeqIdxToBestDiagScore[QSeqIdx] = DiagScore;
#if TRACE
		if (DoTrace(QSeqIdx)) Log("AddTwoHitDiag(TSeqIdx=%u, QSeqIdx=%u, Diag=%u, DiagScore=%d) (better)\n",
								  m_TSeqIdx, QSeqIdx, Diag, DiagScore);
#endif
		}
	}

void PrefilterMu::ExtendTwoHitDiagsToHSPs()
	{
	uint DupeCount = m_DiagBag.m_DupeCount;
	m_NrQueriesWithTwoHitDiag = 0;
	for (uint i = 0; i < DupeCount; ++i)
		{
		uint32_t QSeqIdx = m_DiagBag.m_DupeSeqIdxs[i];
		uint16_t Diag = m_DiagBag.m_DupeDiags[i];
		int DiagScore = ExtendTwoHitDiagToHSP(QSeqIdx, Diag);
		AddTwoHitDiag(QSeqIdx, Diag, DiagScore);
		}
	}

int PrefilterMu::ExtendTwoHitDiagToHSP(uint32_t QSeqIdx, uint16_t Diag)
	{
	const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	const uint QL = m_QDB->GetSeqLength(QSeqIdx);
	int DiagScore = FindHSP(QSeqIdx, Diag);
#if TRACE
	LogDiag(QSeqIdx, Diag);
#endif
	return DiagScore;
	}

void PrefilterMu::Reset()
	{
	for (uint HitIdx = 0; HitIdx < m_NrQueriesWithTwoHitDiag; ++HitIdx)
		{
		uint QSeqIdx = m_QSeqIdxsWithTwoHitDiag[HitIdx];
		m_QSeqIdxToBestDiagScore[QSeqIdx] = 0;
		}
#if DEBUG
	{
	for (uint SeqIdx = 0; SeqIdx < m_QSeqCount; ++SeqIdx)
		{
		assert(m_QSeqIdxToBestDiagScore[SeqIdx] == 0);
		}
	}
#endif
	m_NrQueriesWithTwoHitDiag = 0;
	m_DiagBag.Reset();
	}

void PrefilterMu::LogDiag(uint QSeqIdx, uint16_t Diag) const
	{
	const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	uint QL = m_QDB->GetSeqLength(QSeqIdx);
	string QSeq_ascii;
	for (uint i = 0; i < QL; ++i)
		QSeq_ascii += g_LetterToCharMu[QSeq[i]];
	int Score = FindHSP(QSeqIdx, Diag);
	int Lo, Len;
	int Score2 = FindHSP2(QSeqIdx, Diag, Lo, Len);
	const string &QLabel = m_QDB->GetLabel(QSeqIdx);
	Log("LogDiag(%s, %u) lo %d, len %d, score %d\n",
		QLabel.c_str(), Diag, Lo, Len, Score);
	diag dg(QL, m_TL);
	int ilo = dg.getmini(Diag) + Lo;
	int jlo = dg.getminj(Diag) + Lo;
	string TSeq_ascii;
	for (uint i = 0; i < m_TL; ++i)
		TSeq_ascii += g_LetterToCharMu[m_TSeq[i]];
	Log(" Q %*.*s\n", Len, Len, QSeq_ascii.c_str() + ilo);
	Log(" T %*.*s\n", Len, Len, TSeq_ascii.c_str() + jlo);
	asserta(Score2 == Score);
	}

void PrefilterMu::Search(uint TSeqIdx, const string &TLabel,
				const byte *TSeq, uint TL)
	{
	Search_TargetSeq(TSeqIdx, TLabel, TSeq, TL);

	for (uint i = 0; i < m_NrQueriesWithTwoHitDiag; ++i)
		{
		uint QSeqIdx = m_QSeqIdxsWithTwoHitDiag[i];
		uint16_t DiagScore = m_QSeqIdxToBestDiagScore[QSeqIdx];
		m_RSB.AddScore(QSeqIdx, m_TSeqIdx, DiagScore);
		}
	}

uint PrefilterMu::GetQKmer(uint QSeqIdx, uint QPos) const
	{
	const byte *Q = m_QDB->GetByteSeq(QSeqIdx);
	uint Kmer = m_QKmerIndex->BytesToKmer(Q + QPos);
	return Kmer;
	}

void PrefilterMu::LogQueryKmers(uint QSeqIdx) const
	{
	const byte *Q = m_QDB->GetByteSeq(QSeqIdx);
	const uint QL = m_QDB->GetSeqLength(QSeqIdx);
	Log("\n");
	Log("PrefilterMu::LogQueryKmers() QL=%u >%s\n", 
		QL, m_QDB->GetLabel(QSeqIdx).c_str());
	for (uint PosQ = 0; PosQ + k <= QL; ++PosQ)
		{
		uint Kmer = m_QKmerIndex->BytesToKmer(Q + PosQ);
		string tmp;
		const char *KmerStr = m_QKmerIndex->KmerToStr(Kmer, tmp);
		Log("[%4u]  %08x  %s\n", PosQ, Kmer, KmerStr);
		}
	}

void PrefilterMu::LogTargetKmers() const
	{
	Log("\n");
	Log("PrefilterMu::LogTargetKmers() TL=%u >%s\n", 
		m_TL, m_TLabel);
	for (uint PosT = 0; PosT + k <= m_TL; ++PosT)
		{
		uint Kmer = m_QKmerIndex->BytesToKmer(m_TSeq + PosT);
		string tmp;
		const char *KmerStr = m_QKmerIndex->KmerToStr(Kmer, tmp);
		Log("[%4u]  %08x  %s\n", PosT, Kmer, KmerStr);
		}
	}
