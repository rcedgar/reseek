#include "myutils.h"
#include "prefilter.h"

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
	QKmerIndex.m_MinKmerSelfScore = MIN_KMER_PAIR_SCORE;
	QKmerIndex.FromSeqDB(QDB);
#if DEBUG
	QKmerIndex.Validate();
#endif
	asserta(QKmerIndex.m_k == k);
	asserta(QKmerIndex.m_DictSize == DICT_SIZE);
	asserta(ScoreMx.m_AS_pow[k] == QKmerIndex.m_DictSize);

	vector<vector<int8_t> > BiasVecs8(QSeqCount);
	vector<float> BiasVec;
	const float Scale = 0.15f;
	for (uint QSeqIdx = 0; QSeqIdx < QSeqCount; ++QSeqIdx)
		{
		vector<int8_t> &BiasVec8 = BiasVecs8[QSeqIdx];
		const byte *QSeq = QDB.GetByteSeq(QSeqIdx);
		uint QL = QDB.GetSeqLength(QSeqIdx);
		ScoreMx.CalcLocalBiasCorrection(QSeq, QL, Scale, BiasVec, BiasVec8);
		}

	Prefilter Pref;
	Pref.m_ScoreMx = &ScoreMx;
	Pref.m_QKmerIndex = &QKmerIndex;
	Pref.m_KmerSelfScores = QKmerIndex.m_KmerSelfScores;
	Pref.SetQDB(QDB);
	Pref.m_BiasVecs8 = &BiasVecs8;

	vector<uint> QSeqIdxs;
	vector<int> DiagScores;
	for (uint TSeqIdx = 0; TSeqIdx < TSeqCount; ++TSeqIdx)
		{
		Pref.m_TSeqIdx = TSeqIdx;
		const byte *TSeq = TDB.GetByteSeq(TSeqIdx);
		const string &TLabel = TDB.GetLabel(TSeqIdx);
		uint TL = TDB.GetSeqLength(TSeqIdx);
		Pref.Search(fOut, TSeqIdx, TLabel, TSeq, TL);
		}

	CloseStdioFile(fOut);
	}
