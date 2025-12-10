#pragma once

#include "pdbchain.h"
#include "chaindata.h"
#include "maxbuff.h"

class FragAligner
	{
public:
	static int m_DistScale;

public:
	uint m_Length = 0;
	uint m_BandWidth = 0;
	uint m_IdxCount = 0;
	uint m_DistN = 0;
	const float *m_DistPairScoresPtr = 0;
	vector<uint> m_Idx_is;
	vector<uint> m_Idx_js;

public:
	FragAligner() {}
	~FragAligner() {}

	void Init(uint Length, uint BandWidth, uint DistN,
		const float *DistPairScoresPtr);
	float Align(const uint *DistsPtrA, const uint *DistsPtrB) const;
	uint *GetDistsPtr(const PDBChain &Chain, uint Pos) const;
	void LogScoreMx() const;

public:
	static float GetDistPairScore(uint Dist1, uint Dist2);
	static void GetIdx_ijs(uint Length, uint BandWidth,
		vector<uint> &Idx_is, vector<uint> &Idx_js);
	static float *MakeDistPairScoresPtr(uint DistN);
	};
