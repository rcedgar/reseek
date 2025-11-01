#include "myutils.h"
#include "scop40bench.h"

void SCOP40Bench::MakeSubset(SCOP40Bench &Subset, uint Pct) const
	{
	asserta(Pct > 0 && Pct <= 100);
	Subset.ClearHitsAndResults();

	vector<uint> ChainIdxs;
	const uint ChainCount = GetDBChainCount();
	const uint SubsetChainCount = (ChainCount*Pct)/100;
	asserta(SubsetChainCount > 0 && SubsetChainCount <= ChainCount);
	for (uint i = 0; i < ChainCount; ++i)
		ChainIdxs.push_back(i);
	Shuffle(ChainIdxs);

#define c(x)	Subset.x = x
	c(m_MinEvalue);
	c(m_Level);
	c(m_QuerySelf);
	c(m_SBS);
	c(m_LabelToChainIdx);
#undef c

// Vectors size=ChainCount indexed by ChainIdx
#define c(x)	{ \
	if (SIZE(x) != ChainCount) \
		Die("SCOP40Bench::MakeSubset() size"); \
	Subset.x.clear(); \
	for (uint i = 0; i < SubsetChainCount; ++i) \
		Subset.x.push_back(x[i]); \
	}

	c(m_DBChains);
	c(m_DBProfiles);
	c(m_DBMuLettersVec);
	c(m_DBMuKmersVec);
	c(m_DBSelfRevScores);
#undef c

	Subset.Setup();
	}
