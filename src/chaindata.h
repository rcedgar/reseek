#pragma once

#include "pdbchain.h"
#include "nu.h"
#include "maxbuff.h"

const uint MAXL = 1024;

class ChainData
	{
public:
	const PDBChain *m_Chain = 0;
	uint m_L = 0;
	MaxBuff<float , MAXL> m_DistMxPtr;

public:
	static uint m_MaxAllocatedL;

public:
	void SetChain(const PDBChain &Chain)
		{
		m_Chain = &Chain;
		m_L = Chain.GetSeqLength();
		m_DistMxPtr.m_Filled = false;
		}

	const float* GetDistMxPtr();
	};