#pragma once

class ChainReader2;
class PDBChain;
class DSSParams;
class DSSAligner;

class ProfileLoader
	{
public:
	const DSSParams *m_Params = 0;
	ChainReader2 *m_CR = 0;
	vector<PDBChain *> *m_Chains = 0;
	vector<vector<vector<byte> > *> *m_Profiles = 0;
	vector<vector<byte> *> *m_MuLetters = 0;
	vector<vector<uint> *> *m_KmersVec = 0;
	vector<vector<uint> *> *m_KmerBitsVec = 0;
	vector<float> *m_SelfRevScores = 0;
	mutex m_Lock;
	uint m_Count = 0;
	time_t m_LastProgress = 0;
	uint m_MinChainLength = 1;

public:
	void Load(
	  const DSSParams &Params,
	  ChainReader2 &CR,
	  vector<PDBChain *> *Chains,
	  vector<vector<vector<byte> > *> *Profiles,
	  vector<vector<byte> *> *MuLetters,
	  vector<vector<uint> *> *KmersVec,
	  vector<vector<uint> *> *KmerBitsVec,
	  vector<float> *SelfRevScores,
	  uint ThreadCount);
	void ThreadBody(uint ThreadIndex);

private:
	static void StaticThreadBody(uint ThreadIndex, ProfileLoader *PL);
	};
