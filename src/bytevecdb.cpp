#include "myutils.h"
#include "bytevecdb.h"
#include "pdbchain.h"
#include "dss.h"
#include <omp.h>

void ByteVecDB::Init(
	const vector<PDBChain *> &Chains,
	FEATURE F)
	{
	const uint ThreadCount = GetRequestedThreadCount();
	m_F = F;
	m_ByteVecs.clear();
	m_ptrChains = &Chains;
	m_AlphaSize = DSSParams::GetAlphaSize(m_F);
	const uint ChainCount = SIZE(Chains);
	m_ByteVecs.resize(ChainCount);
	vector<DSS *> Ds(ThreadCount);
	for (uint i = 0; i < ThreadCount; ++i)
		Ds[i] = new DSS;
	Progress("Byte vecs...");
#pragma omp parallel for num_threads(ThreadCount)
	for (int ChainIdx = 0; ChainIdx < int(ChainCount); ++ChainIdx)
		{
		uint ThreadIndex = omp_get_thread_num();
		asserta(ThreadIndex < ThreadCount);
		DSS &D = *Ds[ThreadIndex];
		const PDBChain &Chain = *Chains[ChainIdx];
		D.Init(Chain);
		D.GetByteVec(F, m_ByteVecs[ChainIdx]);
		}
	Progress("done.\n");
	}

uint ByteVecDB::GetTotalLetterCount() const
	{
	uint n = 0;
	for (uint i = 0; i < SIZE(m_ByteVecs); ++i)
		n += SIZE(m_ByteVecs[i]);
	return n;
	}

const vector<uint8_t> &ByteVecDB::GetByteVec(uint ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_ByteVecs));
	return m_ByteVecs[ChainIdx];
	}

void ByteVecDB::SetCounts()
	{
	asserta(m_CountsPtr = 0);
	asserta(m_AlphaSize != 0);
	m_CountsPtr = myalloc(uint, m_AlphaSize);
	zero_array(m_CountsPtr, m_AlphaSize);
	for (uint i = 0; i < SIZE(m_ByteVecs); ++i)
		{
		const vector<uint8_t> &ByteVec = m_ByteVecs[i];
		const uint L = SIZE(ByteVec);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint8_t Letter = ByteVec[Pos];
			if (Letter == 0xff)
				continue;
			asserta(Letter < m_AlphaSize);
			m_CountsPtr[Letter] += 1;
			}
		}
	}

uint ByteVecDB::GetCount(uint8_t Letter) const
	{
	asserta(Letter < m_AlphaSize);
	asserta(m_CountsPtr != 0);
	return m_CountsPtr[Letter];
	}
