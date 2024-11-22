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

	uint m_SeedCountT = 0;
	uint m_SeedCountF = 0;

	uint m_k = 5;
	int m_MinKmerScore = 30;
	uint m_MinSeedHits = 2;
	int m_MinDiagScore = 100;

public:
	void Init(const DSSParams &Params, const string &PairAlnFN,
		const string &ChainsFN);
	void OnPair(uint ChainIdxQ, uint ChainIdxR, bool IsT);
	const byte *GetMuSeq(uint ChainIdx, uint &L);
	bool SeedMatch(const byte *Q, int i, uint LQ,
				   const byte *R, int j, uint LR) const;
	bool GetKmer(const byte *Q, int i, uint QL, uint *Letters) const;
	};

bool Scop40_IsTP_SF(const string &Label1, const string &Label2);
