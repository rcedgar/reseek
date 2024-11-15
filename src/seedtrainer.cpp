#include "myutils.h"
#include "seedtrainer.h"

void SeedTrainer::Init(const DSSParams &Params, const string &PairAlnFN,
	const string &ChainsFN)
	{
	m_Params = &Params;
	m_k = GetPatternOnes(Params.m_PatternStr);
	Trainer::Init(PairAlnFN, ChainsFN);
	uint ChainCount = GetChainCount();
	m_AaKmersVec.resize(ChainCount);
	m_MuKmersVec.resize(ChainCount);
	m_D.SetParams(*m_Params);
	const string &Pattern = m_Params->m_PatternStr;
	uint k = GetPatternOnes(Pattern);
	ProgressLog("k=%u, pattern='%s'\n", k, Pattern.c_str());
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Chain k-mers");
		const PDBChain &Chain = GetChain(ChainIdx);
		m_D.Init(Chain);

		vector<byte> AaLetters;
		vector<byte> MuLetters;
		m_D.GetAaLetters(AaLetters);
		m_D.GetMuLetters(MuLetters);

		m_D.GetAaKmers(AaLetters, m_AaKmersVec[ChainIdx]);
		m_D.GetMuKmers(AaLetters, m_MuKmersVec[ChainIdx]);
		}
	}
