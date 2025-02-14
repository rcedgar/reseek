#include "myutils.h"
#include "prefilter.h"

const MerMx &Get3DiMerMx();

extern int8_t threedi_substmx[20][20];

static const uint k = 6;
static const uint ALPHABET_SIZE = 20;
static const uint DICT_SIZE = 64000000;	// 20^6
static const int MIN_KMER_SELF_SCORE = 78;

//////////////////////////////////////////////
// 	FindHSP searches for the highest-scoring
// 	ungapped alignment on a given diagonal.
//////////////////////////////////////////////
int Prefilter::FindHSP(const byte *QSeq, uint QL, int Diag) const
	{
	asserta(Diag >= 0);

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
		short Score = threedi_substmx[q][t];
		F += Score;
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
int Prefilter::FindHSP2(const byte *QSeq, uint QL,
						int Diag, int &Lo, int &Len) const
	{
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
		short Score = threedi_substmx[q][t];
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

void Prefilter::SetQDB(const SeqDB &QDB)
	{
	m_QDB = &QDB;
	m_QSeqCount = QDB.GetSeqCount();

	m_QSeqIdxToBestDiagScore = myalloc(int, m_QSeqCount);
	m_QSeqIdxsWithTwoHitDiag = myalloc(uint, m_QSeqCount);

	for (uint i = 0; i < m_QSeqCount; ++i)
		{
		m_QSeqIdxToBestDiagScore[i] = INT_MIN;
		m_QSeqIdxsWithTwoHitDiag[i] = UINT_MAX;
		}

	m_NeighborKmers = myalloc(uint, DICT_SIZE);
	m_NrQueriesWithTwoHitDiag = 0;
	}

void Prefilter::Search_TargetKmers()
	{
	if (m_TL < 2*k)
		return;

// Initialize k-mer scan of TSeq
	uint Kmer = 0;
	for (m_TSeqPos = 0; m_TSeqPos < k-1; ++m_TSeqPos)
		{
		byte Letter = m_TSeq[m_TSeqPos];
		assert(Letter < ALPHABET_SIZE);
		Kmer = Kmer*ALPHABET_SIZE + Letter;
		}

// Iterate through k-mers in TSeq
	for (m_TSeqPos = k-1; m_TSeqPos < m_TL; ++m_TSeqPos)
		{
		byte Letter = m_TSeq[m_TSeqPos];
		assert(Letter < ALPHABET_SIZE);
		Kmer = Kmer*ALPHABET_SIZE + Letter;
		Kmer %= DICT_SIZE;
		Search_TargetKmerNeighborhood(Kmer);
		}
	}

void Prefilter::Search_TargetSeq(const byte *TSeq, uint TL,
								 vector<uint> &QSeqIdxs,
								 vector<int> &DiagScores)

	{
	QSeqIdxs.clear();
	DiagScores.clear();

	m_TSeq = TSeq;
	m_TL = TL;

	Search_TargetKmers();
	FindTwoHitDiags();
	ExtendTwoHitDiagsToHSPs();
	GetResults(QSeqIdxs, DiagScores);
	Reset();
	}

void Prefilter::Search_TargetKmerNeighborhood(uint Kmer)
	{
// Construct high-scoring neighborhood
	const uint HSKmerCount =
		m_ScoreMx->GetHighScoring6mers(Kmer, m_KmerNeighborMinScore,
									   m_NeighborKmers);

	for (uint HSKmerIdx = 0; HSKmerIdx < HSKmerCount; ++HSKmerIdx)
		{
		uint HSKmer = m_NeighborKmers[HSKmerIdx];
		Search_Kmer(HSKmer);
		}
	}

void Prefilter::Search_Kmer(uint Kmer)
	{
	uint RowSize = m_QKmerIndex->GetRowSize(Kmer);
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
		uint16_t Diag = dg.getd(QSeqPos, m_TSeqPos);
		m_DiagBag.Add(QSeqIdx, Diag);

		LogDiag(QSeqIdx, Diag);//@@
		}
	}

void Prefilter::FindTwoHitDiags()
	{
	m_DiagBag.ClearDupes();
	m_DiagBag.SetDupes();
#if DEBUG
	m_DiagBag.Validate(m_QSeqCount, INT16_MAX);
#endif
	}

void Prefilter::GetResults(vector<uint> &QSeqIdxs,
						   vector<int> &DiagScores) const
	{
	QSeqIdxs.clear();
	DiagScores.clear();
	for (uint i = 0; i < m_NrQueriesWithTwoHitDiag; ++i)
		{
		uint QSeqIdx = m_QSeqIdxsWithTwoHitDiag[i];
		int DiagScore = m_QSeqIdxToBestDiagScore[i];

		QSeqIdxs.push_back(QSeqIdx);
		DiagScores.push_back(DiagScore);
		}
	}

void Prefilter::AddTwoHitDiag(uint QSeqIdx, uint16_t Diag, int DiagScore)
	{
	int BestDiagScoreT = m_QSeqIdxToBestDiagScore[QSeqIdx];
	if (BestDiagScoreT == INT_MIN)
		{
		m_QSeqIdxsWithTwoHitDiag[m_NrQueriesWithTwoHitDiag++] = QSeqIdx;
		m_QSeqIdxToBestDiagScore[QSeqIdx] = DiagScore;
		}
	else if (DiagScore > m_QSeqIdxToBestDiagScore[QSeqIdx])
		m_QSeqIdxToBestDiagScore[QSeqIdx] = DiagScore;
	}

void Prefilter::ExtendTwoHitDiagsToHSPs()
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

int Prefilter::ExtendTwoHitDiagToHSP(uint32_t QSeqIdx, uint16_t Diag)
	{
	const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	const uint QL = m_QDB->GetSeqLength(QSeqIdx);
	int DiagScore = FindHSP(QSeq, QL, Diag);
	return DiagScore;
	}

void Prefilter::Reset()
	{
	for (uint HitIdx = 0; HitIdx < m_NrQueriesWithTwoHitDiag; ++HitIdx)
		{
		uint QSeqIdx = m_QSeqIdxsWithTwoHitDiag[HitIdx];
		m_QSeqIdxToBestDiagScore[QSeqIdx] = INT_MIN;
		}
#if DEBUG
	{
	for (uint SeqIdx = 0; SeqIdx < m_QSeqCount; ++SeqIdx)
		{
		assert(m_QSeqIdxToBestDiagScore[SeqIdx] == INT_MIN);
		}
	}
#endif
	m_NrQueriesWithTwoHitDiag = 0;
	m_DiagBag.Reset();
	}

void Prefilter::LogDiag(uint QSeqIdx, uint16_t Diag) const
	{
	const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	uint QL = m_QDB->GetSeqLength(QSeqIdx);
	int Score = FindHSP(QSeq, QL, Diag);
	int Lo, Len;
	int Score2 = FindHSP2(QSeq, QL, Diag, Lo, Len);
	const string &QLabel = m_QDB->GetLabel(QSeqIdx);
	Log("LogDiag(%s, %u) lo %d, len %d, score %d\n",
		QLabel.c_str(), Diag, Lo, Len, Score);
	asserta(Score2 == Score);
	}

void cmd_prefilter_3di()
	{
	const string &Query3Di_FN = g_Arg1;
	const string &DB3Di_FN = opt_db;

	FILE *fOut = CreateStdioFile(opt_output);
	
	SeqDB QDB;
	SeqDB TDB;

	QDB.FromFasta(Query3Di_FN);
	TDB.FromFasta(DB3Di_FN);

	QDB.ToLetters(g_CharToLetterAmino);
	TDB.ToLetters(g_CharToLetterAmino);
	const uint QSeqCount = QDB.GetSeqCount();
	const uint TSeqCount = TDB.GetSeqCount();

	const MerMx &ScoreMx = Get3DiMerMx();
	asserta(ScoreMx.m_k == k);

	ThreeDex QKmerIndex;
	QKmerIndex.m_KmerSelfScores = ScoreMx.BuildSelfScores_6mers();
	QKmerIndex.m_MinKmerSelfScore = MIN_KMER_SELF_SCORE;
	QKmerIndex.FromSeqDB(QDB);
	QKmerIndex.Validate();
	QKmerIndex.LogStats();
	asserta(QKmerIndex.m_k == k);
	asserta(QKmerIndex.m_DictSize == DICT_SIZE);
	asserta(ScoreMx.m_AS_pow[k] == QKmerIndex.m_DictSize);

	Prefilter Pref;
	Pref.m_ScoreMx = &ScoreMx;
	Pref.m_QKmerIndex = &QKmerIndex;
	Pref.SetQDB(QDB);

	vector<uint> QSeqIdxs;
	vector<int> DiagScores;
	for (uint TSeqIdx = 0; TSeqIdx < TSeqCount; ++TSeqIdx)
		{
		const byte *TSeq = TDB.GetByteSeq(TSeqIdx);
		const string &TLabel = TDB.GetLabel(TSeqIdx);
		uint TL = TDB.GetSeqLength(TSeqIdx);
		Pref.Search_TargetSeq(TSeq, TL, QSeqIdxs, DiagScores);

		if (fOut != 0)
			{
			const uint n = SIZE(QSeqIdxs);
			asserta(SIZE(DiagScores) == n);
			for (uint i = 0; i < n; ++i)
				{
				uint QSeqIdx = QSeqIdxs[i];
				int DiagScore = DiagScores[i];

				const string &QLabel = QDB.GetLabel(QSeqIdx);

				fprintf(fOut, "%s", TLabel.c_str());
				fprintf(fOut, "\t%d", DiagScore);
				fprintf(fOut, "\t%s", QLabel.c_str());
				fprintf(fOut, "\n");
				}
			}
		}

	CloseStdioFile(fOut);
	}
