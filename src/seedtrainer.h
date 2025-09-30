#if 0 // @@DELETE
#pragma once

#include "trainer.h"
#include "dss.h"

class SeedTrainer : public Trainer
	{
public:
	const DSSParams *m_Params = 0;
	uint m_k = 0;
	float m_AaWeight = 1;
	uint m_MaxDiagDist = 40; // BLASTP value = 40
	vector<vector<uint> > m_AaKmersVec;
	vector<vector<uint> > m_MuKmersVec;
	DSS m_D;
	vector<float> m_PosScores;
	vector<float> m_NegScores;

	vector<float> m_TwoHitIdenticalMuKmer_MinScores_T;
	vector<float> m_TwoHitIdenticalMuKmer_MaxScores_F;

	uint m_NT = 0;
	uint m_NF = 0;

public:
	void Init(const DSSParams &Params, const string &PairAlnFN,
		const string &ChainsFN);
	void GetAlignedKmerStarts(uint PairIdx,
		vector<uint> &PosQs, vector<uint> &PosRs) const;
	float GetAaScore(byte Letter1, byte Letter2) const;
	float GetMuScore(byte Letter1, byte Letter2) const;
	float GetAaKmerScore(uint Kmer1, uint Kmer2) const;
	float GetMuKmerScore(uint Kmer1, uint Kmer2) const;
	void FindIdenticalMuKmerPairsSameDiag(uint ChainIdxQ, uint ChainIdxR,
	  float &MinScore, float &MaxScore) const;
	bool IsPairSameDiag(const vector<uint> &PosQs1, const vector<uint> &PosRs1,
		const vector<uint> &PosQs2, const vector<uint> &PosRs2) const;
	void OnPairT(uint PairIdx);
	void OnPairF(uint ChainIdxQ, uint ChainIdxR);
	};

bool Scop40_IsTP_SF(const string &Label1, const string &Label2);
#endif