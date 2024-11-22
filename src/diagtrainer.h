#pragma once

#include "trainer.h"
#include "dss.h"
#include "diaghsp.h"

class DiagTrainer : public Trainer
	{
public:
	const DSSParams *m_Params = 0;
	DSS m_D;

	vector<float> m_PosScores;
	vector<float> m_NegScores;

	uint m_NT = 0;
	uint m_NF = 0;

	DiagHSP m_DH;

public:
	void Init(const DSSParams &Params, const string &PairAlnFN,
		const string &ChainsFN);
	void OnPair(uint ChainIdxQ, uint ChainIdxR, bool IsT);
	const byte *GetMuSeq(uint ChainIdx, uint &L);
	};

bool Scop40_IsTP_SF(const string &Label1, const string &Label2);
