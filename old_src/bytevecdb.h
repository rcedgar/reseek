#pragma once

#include "features.h"
#include "dssparams.h"

class PDBChain;

class ByteVecDB
	{
public:
	FEATURE m_F = FEATURE(-1);
	uint m_AlphaSize = 0;
	const vector<PDBChain *> *m_ptrChains = 0;
	vector<vector<uint8_t> > m_ByteVecs;
	uint *m_CountsPtr = 0;

public:
	void Init(const vector<PDBChain *> &Chains, FEATURE F);
	uint GetTotalLetterCount() const;
	const vector<uint8_t> &GetByteVec(uint ChainIdx) const;
	uint GetAlphaSize() const { return m_AlphaSize; }
	uint GetSeqCount() const { return SIZE(m_ByteVecs); }
	void SetCounts();
	uint GetCount(uint8_t Letter) const;
	};