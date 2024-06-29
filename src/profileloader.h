#pragma once

class ChainReader2;
class PDBChain;
class DSSParams;

class ProfileLoader
	{
public:
	const DSSParams *m_Params = 0;
	ChainReader2 *m_CR = 0;
	vector<PDBChain *> *m_Chains;
	vector<vector<vector<byte> > *> *m_Profiles;
	mutex m_Lock;

public:
	void Load(
	  ChainReader2 &CR,
	  uint ReserveSize,
	  vector<PDBChain *> &Chains,
	  vector<vector<vector<byte> > *> &m_Profiles,
	  const DSSParams &Params,
	  uint ThreadCount);
	void ThreadBody(uint ThreadIndex);

private:
	static void StaticThreadBody(uint ThreadIndex, ProfileLoader *PL);
	};
