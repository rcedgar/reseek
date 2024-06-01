#pragma once

#include "seqdb.h"
#include "pdbchain.h"
#include "logodds.h"

class Trainer;

typedef void TRAINER_ONPAIR(
  const Trainer &T, uint ChainIdxQ, uint ChainIdxR,
  const vector<uint> &PosQs, const vector<uint> &PosRs);

typedef void TRAINER_ONCOL(
  const Trainer &T, uint PosQ, uint PosR);

typedef void TRAINER_ALPHACOL(
  const Trainer &T, uint PosQ, uint PosR,
  uint &LetterQ, uint &LetterR);

class Trainer
	{
public:
	SeqDB m_PairAlnDB;
	vector<PDBChain *> m_Chains;
	map<string, uint> m_DomToChainIndex;
	vector<uint> m_ChainIdxsQ;
	vector<uint> m_ChainIdxsR;
	vector<vector<uint> > m_PosVecsQ;
	vector<vector<uint> > m_PosVecsR;

public:
	void Init(
	  const string &PairAlnFN,
	  const string &ChainsFN);
	uint GetPairCount() const;
	void Scan(TRAINER_ONPAIR OnPair, TRAINER_ONCOL OnCol) const;
	void TrainLogOdds(uint AlphaSize,
	  TRAINER_ONPAIR OnPair,
	  TRAINER_ALPHACOL AlphaCol,
	  LogOdds &LO) const;
	const PDBChain &GetChain(uint ChainIdx) const;
	};
