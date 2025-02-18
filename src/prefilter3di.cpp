#include "myutils.h"
#include "prefilter3di.h"

mutex Prefilter3Di::m_Lock;
RankedScoresBag Prefilter3Di::m_RSB;

//////////////////////////////////////////////
// 	FindHSP searches for the highest-scoring
// 	ungapped alignment on a given diagonal.
//////////////////////////////////////////////
int Prefilter3Di::FindHSP(uint QSeqIdx, int Diag) const
	{
	const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	const uint QL = m_QDB->GetSeqLength(QSeqIdx);

	asserta(Diag >= 0);
#if TRACE
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
		short Score = threedi_substmx2[q][t];
		F += Score;
#if TRACE
		if (DoTrace(QSeqIdx)) Log(" i=%u j=%u F=%d B=%d score=%d\n", i, j, F, B, Score);
#endif
		if (F > B)
			B = F;
		else if (F < 0)
			F = 0;
		}
	return B;
	}

int Prefilter3Di::FindHSP_Biased(uint QSeqIdx, const vector<int8_t> &BiasVec8,
							  uint Qlo, uint Tlo) const
	{
	const byte *TSeq = m_TSeq + Tlo;
	const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx) + Qlo;
	const int8_t *QBias = BiasVec8.data() + Qlo;
	const uint QL = m_QDB->GetSeqLength(QSeqIdx);

	asserta(SIZE(BiasVec8) == QL);
#if TRACE
	if (DoTrace(QSeqIdx)) Log("FindHSP_Biased(Qlo=%u, Tlo=%u)\n", Qlo, Tlo);
#endif

	int nq = QL - Qlo;
	int nt = m_TL - Tlo;
	int n = min(nq, nt);
	uint i = Qlo;
	uint j = Tlo;

	int B = 0;
	int F = 0;
	for (int k = 0; k < n; ++k)
		{
		assert(i < int(QL));
		assert(j < int(m_TL));
		int8_t Bias = *QBias++;
		byte q = *QSeq++;
		byte t = *TSeq++;
		assert(q < ALPHABET_SIZE);
		assert(t < ALPHABET_SIZE);
		short Score = threedi_substmx2[q][t] + Bias;
		F += Score;
#if TRACE
		if (DoTrace(QSeqIdx)) Log(" k=%u F=%d B=%d score=%d\n", k, F, B, Score);
#endif
		if (F > B)
			B = F;
		else if (F < 0)
			F = 0;
		}
#if TRACE
	if (DoTrace(QSeqIdx)) Log("FindHSP_Biased(Qlo=%u, Tlo=%u)=%d\n", Qlo, Tlo, B);
#endif
	return B;
	}

//////////////////////////////////////////////
// 	FindHSP plus "traceback", i.e. returns
// 	start position and length of HSP.
//////////////////////////////////////////////
int Prefilter3Di::FindHSP2(uint QSeqIdx,
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
		short Score = threedi_substmx2[q][t];
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

void Prefilter3Di::SetQDB(const SeqDB &QDB)
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

	m_NeighborKmers = myalloc(uint, DICT_SIZE);
	m_NrQueriesWithTwoHitDiag = 0;
	}

void Prefilter3Di::Search_TargetKmers()
	{
	if (m_TL < 2*k)
		return;

	m_QKmerIndex->GetKmers(m_TSeq, m_TL, m_TKmers);
	const uint NK = SIZE(m_TKmers);
	for (uint i = 0; i < NK; ++i)
		{
		uint Kmer = m_TKmers[i];
		Search_TargetKmerNeighborhood(Kmer, i);
		}
	}

void Prefilter3Di::Search_TargetSeq()
	{
	Reset();
	Search_TargetKmers();
	FindTwoHitDiags();
	ExtendTwoHitDiagsToHSPs();
	}

void Prefilter3Di::Search_TargetKmerNeighborhood(uint Kmer, uint TPos)
	{
	if (Kmer == UINT_MAX)
		return;
#if TRACE
	m_TBaseKmer = Kmer;
#endif
	short Bias = GetTargetBiasCorrection(TPos);
	assert(Kmer < DICT_SIZE);
	assert(m_KmerSelfScores[Kmer] >= MIN_KMER_PAIR_SCORE);

// Construct high-scoring neighborhood
	short MinKmerScore = MIN_KMER_PAIR_SCORE - Bias;
	const uint HSKmerCount =
		m_ScoreMx->GetHighScoring6mers(Kmer, MinKmerScore, m_NeighborKmers);
#if TRACE
	string Tmp;
	Log("Search_TargetKmerNeighborhood TPos=%u Kmer=%s bias=%d minscore=%d |nbrs|=%d\n",
		TPos, KmerToStr(Kmer, Tmp), Bias, MinKmerScore, HSKmerCount);
#endif
	for (uint HSKmerIdx = 0; HSKmerIdx < HSKmerCount; ++HSKmerIdx)
		{
		uint HSKmer = m_NeighborKmers[HSKmerIdx];
		Search_Kmer(HSKmer, TPos);
		}
	}

void Prefilter3Di::Search_Kmer(uint Kmer, uint TPos)
	{
	uint RowSize = m_QKmerIndex->GetRowSize(Kmer);
	if (RowSize == 0)
		return;
	uint DataOffset = m_QKmerIndex->GetRowStart(Kmer);
	for (uint ColIdx = 0; ColIdx < RowSize; ++ColIdx)
		{
		uint32_t QSeqIdx;
		uint16_t QSeqPos;
		m_QKmerIndex->Get(DataOffset++, QSeqIdx, QSeqPos);
		uint QL32 = m_QDB->GetSeqLength(QSeqIdx);
		asserta(QL32 < UINT16_MAX);
		uint16_t QL = uint16_t(QL32);
		diag dg(QL, m_TL);
		uint16_t Diag = dg.getd(QSeqPos, TPos);
#if TRACE
		{
		string TKmerStr;
		string QKmerStr;
		m_QKmerIndex->KmerToStr(m_TBaseKmer, TKmerStr);
		m_QKmerIndex->KmerToStr(Kmer, QKmerStr);
		if (DoTrace(QSeqIdx)) Log("m_DiagBag(QSeqIdx=%u, Diag=%u) Q%u=%s T%u=%s\n",
								  QSeqIdx, Diag, QSeqPos, QKmerStr.c_str(), TPos, TKmerStr.c_str());
		}
#endif
		m_DiagBag.Add(QSeqIdx, Diag);
		}
	}

void Prefilter3Di::FindTwoHitDiags()
	{
	m_DiagBag.ClearDupes();
	m_DiagBag.SetDupes();
#if DEBUG
	m_DiagBag.Validate(m_QSeqCount, INT16_MAX);
#endif
	}

void Prefilter3Di::GetResults(vector<uint> &QSeqIdxs,
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

void Prefilter3Di::AddTwoHitDiag(uint QSeqIdx, uint16_t Diag, int DiagScore)
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

void Prefilter3Di::ExtendTwoHitDiagsToHSPs()
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

int Prefilter3Di::ExtendTwoHitDiagToHSP(uint32_t QSeqIdx, uint16_t Diag)
	{
	const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	const uint QL = m_QDB->GetSeqLength(QSeqIdx);
	assert(m_BiasVecs8 != 0);
	assert(QSeqIdx < SIZE(*m_BiasVecs8));
	const vector<int8_t> &BiasVec8 = (*m_BiasVecs8)[QSeqIdx];
	diag dg(QL, m_TL);
	uint Qlo = dg.getmini(Diag);
	uint Tlo = dg.getminj(Diag);
	int DiagScore = FindHSP_Biased(QSeqIdx, BiasVec8, Qlo, Tlo);
	return DiagScore;
	}

void Prefilter3Di::Reset()
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

// lib/mmseqs/src/prefiltering/QueryMatcher.cpp:257
short Prefilter3Di::GetTargetBiasCorrection(uint TPos) const
	{
	const byte *Offsets = ThreeDex::m_Offsets;
	const uint k = ThreeDex::m_k;
	float FloatBias = 0;
	for (uint i = 0; i < k; ++i)
		FloatBias += m_TBiasVec[TPos + Offsets[i]];
	short ShortBias = short(FloatBias < 0 ? FloatBias - 0.5 : FloatBias + 0.5);
	return ShortBias;
	}

void Prefilter3Di::SetTarget(uint TSeqIdx, const string &TLabel,
	const byte *TSeq, uint TL)
	{
	m_TSeqIdx = TSeqIdx;
	m_TLabel = TLabel;
	m_TSeq = TSeq;
	m_TL = TL;
	CalcLocalBiasCorrection_3Di(TSeq, TL, BIAS_WINDOW, TBIAS_SCALE, m_TBiasVec, m_TBiasVec8);

#if TRACE
	Log("\n");
	Log("SetTarget(%u) >%s L=%u\n", TSeqIdx, TLabel.c_str(), TL);
	for (uint i = 0; i < TL; ++i)
		{
		if (i > 0 && i%80 == 0)
			Log("\n");
		Log("%c", g_LetterToCharAmino[TSeq[i]]);
		}
	Log("\n");
	Log("Bias ");
	for (uint i = 0; i < TL; ++i)
		{
		if (i > 0 && i%80 == 0)
			Log("\n");
		Log(" %.3g", m_TBiasVec[i]);
		}
	Log("\n");
#endif
	}

void Prefilter3Di::LogDiag(uint QSeqIdx, uint16_t Diag) const
	{
	const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	uint QL = m_QDB->GetSeqLength(QSeqIdx);
	int Score = FindHSP(QSeqIdx, Diag);
	int Lo, Len;
	int Score2 = FindHSP2(QSeqIdx, Diag, Lo, Len);
	const string &QLabel = m_QDB->GetLabel(QSeqIdx);
	Log("LogDiag(%s, %u) lo %d, len %d, score %d\n",
		QLabel.c_str(), Diag, Lo, Len, Score);
	asserta(Score2 == Score);
	}

void Prefilter3Di::Search(FILE *fTsv, uint TSeqIdx, const string &TLabel,
				const byte *TSeq, uint TL)
	{
	SetTarget(TSeqIdx, TLabel, TSeq, TL);
	Search_TargetSeq();

	for (uint i = 0; i < m_NrQueriesWithTwoHitDiag; ++i)
		{
		uint QSeqIdx = m_QSeqIdxsWithTwoHitDiag[i];
		uint16_t DiagScore = m_QSeqIdxToBestDiagScore[QSeqIdx];
		m_RSB.AddScore(QSeqIdx, m_TSeqIdx, DiagScore);
		}

	if (fTsv == 0)
		return;

 	for (uint i = 0; i < m_NrQueriesWithTwoHitDiag; ++i)
		{
		uint QSeqIdx = m_QSeqIdxsWithTwoHitDiag[i];
		uint16_t DiagScore = m_QSeqIdxToBestDiagScore[QSeqIdx];
		const string &QLabel = m_QDB->GetLabel(QSeqIdx);

		m_Lock.lock();
		fprintf(fTsv, "%s", m_TLabel.c_str());
		fprintf(fTsv, "\t%s", QLabel.c_str());
		fprintf(fTsv, "\t%d", DiagScore);
		fprintf(fTsv, "\n");
		m_Lock.unlock();
		}
	}
