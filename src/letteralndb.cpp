#include "myutils.h"
#include "letteralndb.h"
#include "pwalndb.h"
#include "bytevecdb.h"

void LetterAlnDB::Init(const ByteVecDB &BVDB, const PWAlnDB &PADB)
	{
	m_AlphaSize = BVDB.GetAlphaSize();
	m_ByteVecDB = &BVDB;
	m_PWAlnDB = &PADB;
	const uint SeqCount = m_ByteVecDB->GetSeqCount();
	m_LetterPairCount = m_PWAlnDB->GetTotalColCount();
	m_Letters = myalloc(uint8_t, 2*m_LetterPairCount);
	const uint AlnCount = PADB.GetAlnCount();
	uint Offset = 0;
	for (uint AlnIdx = 0; AlnIdx < AlnCount; ++AlnIdx)
		{
		uint Idx1 = PADB.GetLabelIdx1(AlnIdx);
		uint Idx2 = PADB.GetLabelIdx2(AlnIdx);
		const vector<uint8_t> &BV1 = BVDB.GetByteVec(Idx1);
		const vector<uint8_t> &BV2 = BVDB.GetByteVec(Idx2);
		const uint L1 = SIZE(BV1);
		const uint L2 = SIZE(BV2);
		const vector<uint> &ColToPos1 = PADB.GetColToPosVec1(AlnIdx);
		const vector<uint> &ColToPos2 = PADB.GetColToPosVec2(AlnIdx);
		const uint ColCount = SIZE(ColToPos1);
		asserta(SIZE(ColToPos2) == ColCount);
		for (uint ColIdx = 0; ColIdx < ColCount; ++ColIdx)
			{
			uint Pos1 = ColToPos1[ColIdx];
			uint Pos2 = ColToPos2[ColIdx];
			asserta(Pos1 < L1);
			asserta(Pos2 < L2);
			uint8_t Letter1 = BV1[Pos1];
			uint8_t Letter2 = BV2[Pos2];
			if (Letter1 == 0xff || Letter2 == 0xff)
				continue;
			asserta(Letter1 < m_AlphaSize);
			asserta(Letter2 < m_AlphaSize);
			m_Letters[Offset++] = Letter1;
			m_Letters[Offset++] = Letter2;
			}
		}
	asserta(Offset == 2*m_LetterPairCount);
	}

uint LetterAlnDB::GetCount(uint8_t Letter) const
	{
	asserta(m_CountsPtr != 0);
	asserta(Letter < m_AlphaSize);
	return m_CountsPtr[Letter];
	}

float LetterAlnDB::GetFreq(uint8_t Letter) const
	{
	asserta(m_FreqsPtr != 0);
	asserta(Letter < m_AlphaSize);
	return m_FreqsPtr[Letter];
	}

void LetterAlnDB::SetCountsAndFreqs(uint PsuedoCount)
	{
	asserta(m_CountsPtr == 0);
	asserta(m_FreqsPtr == 0);
	asserta(m_JointCountsMxPtr == 0);
	asserta(m_JointFreqsMxPtr == 0);
	uint AS = GetAlphaSize();
	m_CountsPtr = myalloc(uint, AS);
	m_FreqsPtr = myalloc(float, AS);
	m_JointCountsMxPtr = myalloc(uint, AS*AS);
	m_JointFreqsMxPtr = myalloc(float, AS*AS);
	zero_array(m_CountsPtr, AS);
	zero_array(m_FreqsPtr, AS);
	zero_array(m_JointCountsMxPtr, AS*AS);
	zero_array(m_JointFreqsMxPtr, AS*AS);

	for (uint i = 0; i < AS; ++i)
		m_CountsPtr[i] = PsuedoCount;

	for (uint i = 0; i < AS; ++i)
		for (uint j = 0; j < AS; ++j)
			m_JointCountsMxPtr[AS*i + j] = PsuedoCount;

	for (uint i = 0; i < m_LetterPairCount; ++i)
		{
		uint8_t Letter1 = m_Letters[2*i];
		uint8_t Letter2 = m_Letters[2*i+1];
		if (Letter1 == 0xff || Letter2 == 0xff)
			continue;
		asserta(Letter1 < AS);
		asserta(Letter2 < AS);
		m_CountsPtr[Letter1] += 1;
		m_CountsPtr[Letter2] += 1;
		m_JointCountsMxPtr[AS*Letter1 + Letter2] += 1;
		m_JointCountsMxPtr[AS*Letter2 + Letter1] += 1;
		}

	uint N2 = 2*m_LetterPairCount + m_AlphaSize*m_AlphaSize*PsuedoCount;
	float SumFreq = 0;
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		float Freq = float(m_CountsPtr[i])/N2;
		m_FreqsPtr[i] = Freq;
		SumFreq += Freq;
		}
	asserta(feq(SumFreq, 1));

	uint N1 = 2*m_LetterPairCount + m_AlphaSize*PsuedoCount;
	float SumFreq2 = 0;
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		for (uint j = 0; j < m_AlphaSize; ++j)
			{
			uint n = m_JointCountsMxPtr[AS*i + j];
			float Freq = float(n)/N1;
			m_JointFreqsMxPtr[AS*i + j] = Freq;
			SumFreq2 += Freq;
			}
		}
	asserta(feq(SumFreq2, 1));
	}

uint LetterAlnDB::GetJointCount(uint8_t Letter1, uint8_t Letter2) const
	{
	asserta(m_JointCountsMxPtr != 0);
	asserta(Letter1 < m_AlphaSize);
	asserta(Letter2 < m_AlphaSize);
	return m_JointCountsMxPtr[m_AlphaSize*Letter1 + Letter2];
	}

float LetterAlnDB::GetJointFreq(uint8_t Letter1, uint8_t Letter2) const
	{
	asserta(m_JointFreqsMxPtr != 0);
	asserta(Letter1 < m_AlphaSize);
	asserta(Letter2 < m_AlphaSize);
	return m_JointFreqsMxPtr[m_AlphaSize*Letter1 + Letter2];
	}

float LetterAlnDB::GetExpectedFreq(const float *LetterFreqsPtr,
	uint8_t Letter1, uint8_t Letter2) const
	{
	float ObsFreq1 = LetterFreqsPtr[Letter1];
	float ObsFreq2 = LetterFreqsPtr[Letter2];
	float ExpectedFreq = ObsFreq1*ObsFreq2;
	return ExpectedFreq;
	}

float LetterAlnDB::GetExpectedCount(const float *LetterFreqsPtr,
	uint8_t Letter1, uint8_t Letter2) const
	{
	float ExpectedFreq = GetExpectedFreq(LetterFreqsPtr, Letter1, Letter2);
	float ExpectedCount = ExpectedFreq*m_LetterPairCount;
	return ExpectedCount;
	}

float LetterAlnDB::GetLogOddsScore(uint8_t Letter1, uint8_t Letter2,
	const float *LetterFreqsPtr) const
	{
	uint ObsPairCount = GetJointCount(Letter1, Letter2);
	float ObsFreq = float(ObsPairCount)/m_LetterPairCount;
	float ObsFreq1 = LetterFreqsPtr[Letter1];
	float ObsFreq2 = LetterFreqsPtr[Letter2];
	float ExpectedFreq = GetExpectedFreq(LetterFreqsPtr, Letter1, Letter2);
	float Ratio = ObsFreq/ExpectedFreq;
	float LogOddsScore = log(Ratio);
	return LogOddsScore;
	}

float *LetterAlnDB::GetLogOddsMxPtr(
	const float *LetterFreqsPtr) const
	{
	float *Mx = myalloc(float, m_AlphaSize*m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		for (uint j = i; j < m_AlphaSize; ++j)
			{
			float Score = GetLogOddsScore(i, j, LetterFreqsPtr);
			Mx[m_AlphaSize*i + j] = Score;
			if (i != j)
				Mx[m_AlphaSize*j + i] = Score;
			}
	return Mx;
	}

/***
Expected Score.
Kullback-Leibler divergence between joint probabilities for homologs
and joint probabilities for non-homologs (random pairs) is precisely
the expected log-odds score when you draw pairs from the homolog distribution.
***/
float LetterAlnDB::GetES(const float *LetterFreqsPtr) const
	{
	float Sum = 0;
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		for (uint j = 0; j < m_AlphaSize; ++j)
			{
			float P_ij = GetJointFreq(i, j);
			float Score = GetLogOddsScore(i, j, LetterFreqsPtr);
			Sum += P_ij*Score;
			}
		}
	return Sum;
	}

float LetterAlnDB::GetES2(
	const uint AlphaSize,
	const float *LogOddsMxPtr,
	const float *JointFreqsPtr)
	{
	float Sum = 0;
	for (uint i = 0; i < AlphaSize; ++i)
		{
		for (uint j = 0; j < AlphaSize; ++j)
			{
			float P_ij = JointFreqsPtr[AlphaSize*i + j];
			float Score = LogOddsMxPtr[AlphaSize*i + j];
			Sum += P_ij*Score;
			}
		}
	return Sum;
	}

float LetterAlnDB::GetLogOddsScore2(
	uint8_t Letter1, uint8_t Letter2,
	const LetterAlnDB &DB_TP,
	const LetterAlnDB &DB_FP,
	const float *LetterFreqsPtr)
	{
	float Freq_TP = DB_TP.GetJointFreq(Letter1, Letter2);
	float Freq_FP = DB_FP.GetJointFreq(Letter1, Letter2);
	float Ratio = Freq_TP/Freq_FP;
	float LogOddsScore = log(Ratio);
	return LogOddsScore;
	}

float *LetterAlnDB::GetLogOddsMxPtr2(
	uint AlphaSize,
	const LetterAlnDB &LA_TP,
	const LetterAlnDB &LA_FP,
	const float *LetterFreqsPtr)
	{
	float *Mx = myalloc(float, AlphaSize*AlphaSize);
	for (uint i = 0; i < AlphaSize; ++i)
		for (uint j = i; j < AlphaSize; ++j)
			{
			float Score = GetLogOddsScore2(i, j, LA_TP, LA_FP, LetterFreqsPtr);
			Mx[AlphaSize*i + j] = Score;
			if (i != j)
				Mx[AlphaSize*j + i] = Score;
			}
	return Mx;
	}
