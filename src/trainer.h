#pragma once

#include "seqdb.h"
#include "pdbchain.h"
#include "logodds.h"
#include "dssparams.h"

class Trainer;

typedef bool TRAINER_IS_TP(const string &Label1, const string &Label2);

typedef void TRAINER_ONPAIR_T(Trainer &T, uint PairIdx);
typedef void TRAINER_ONPAIR_F(Trainer &T, uint ChainIdxQ, uint ChainIdxR);

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
	vector<string> m_RowsQ;
	vector<string> m_RowsR;
	vector<vector<uint> > m_PosVecsQ;
	vector<vector<uint> > m_PosVecsR;

public:
	void Init(const string &PairAlnFN, const string &ChainsFN);
	uint GetPairCount() const;
	const uint GetChainCount() const { return SIZE(m_Chains); }
	const string &GetLabel(uint ChainIdx) const
		{
		asserta(ChainIdx < SIZE(m_Chains));
		return m_Chains[ChainIdx]->m_Label;
		}
	uint GetSeqLength(uint ChainIdx) const;
	void Scan(TRAINER_ONPAIR OnPair, TRAINER_ONCOL OnCol) const;
	void TrainLogOdds(uint AlphaSize, TRAINER_ONPAIR OnPair,
	  TRAINER_ALPHACOL AlphaCol, LogOdds &LO) const;
	const PDBChain &GetChain(uint ChainIdx) const;
	void EnumChainPairsT(TRAINER_ONPAIR_T OnPair);
	void EnumChainPairsF(TRAINER_IS_TP IsTP, TRAINER_ONPAIR_F OnPair);
	};
