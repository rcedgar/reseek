#pragma once

#include "bytevecdb.h"
#include "pwalndb.h"

class LetterAlnDB
	{
public:
	const ByteVecDB *m_ByteVecDB = 0;
	const PWAlnDB *m_PWAlnDB = 0;
	uint8_t *m_Letters = 0;
	uint m_LetterPairCount = 0;
	uint m_AlphaSize = 0;
	uint *m_CountsPtr = 0;
	float *m_FreqsPtr = 0;
	uint *m_JointCountsMxPtr = 0;
	float *m_JointFreqsMxPtr = 0;

public:
	void Init(const ByteVecDB &BVDB, const PWAlnDB &PADB);
	uint GetAlphaSize() const { return m_ByteVecDB->GetAlphaSize(); }
	void SetCountsAndFreqs(uint PseudoCount);
	uint GetCount(uint8_t Letter) const;
	float GetFreq(uint8_t Letter) const;
	uint GetJointCount(uint8_t Letter1, uint8_t Letter2) const;
	float GetJointFreq(uint8_t Letter1, uint8_t Letter2) const;
	float GetExpectedCount(const float *LetterFreqsPtr,
		uint8_t Letter1, uint8_t Letter2) const;
	float GetExpectedFreq(const float *LetterFreqsPtr,
		uint8_t Letter1, uint8_t Letter2) const;
	float GetLogOddsScore(uint8_t Letter1, uint8_t Letter2,
		const float *LetterFreqsPtr) const;
	float *GetLogOddsMxPtr(
		const float *LetterFreqsPtr) const;
	float GetES(const float *LetterFreqsPtr) const;

public:
	static float GetES2(
		const uint AlphaSize,
		const float *LogOddsMxPtr,
		const float *JointFreqsPtr);

	static float GetLogOddsScore2(
		uint8_t Letter1, uint8_t Letter2,
		const LetterAlnDB &DB_TP,
		const LetterAlnDB &DB_FP,
		const float *LetterFreqsPtr);

	static float *GetLogOddsMxPtr2(
		uint AlphaSize,
		const LetterAlnDB &LA_TP,
		const LetterAlnDB &LA_FP,
		const float *LetterFreqsPtr);
	};