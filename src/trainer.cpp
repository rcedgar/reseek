#include "myutils.h"
#include "trainer.h"
#include "scop40bench.h"
#include "alpha.h"

void Trainer::Init(
  const string &PairAlnFN,
  const string &ChainsFN)
	{
	asserta(m_Chains.empty());
	asserta(m_PairAlnDB.GetSeqCount() == 0);
	asserta(m_ChainIdxsQ.empty());
	asserta(m_ChainIdxsR.empty());
	asserta(m_PosVecsQ.empty());
	asserta(m_PosVecsR.empty());

	m_PairAlnDB.FromFasta(PairAlnFN, true);

	ReadChains(ChainsFN, m_Chains);
	const uint ChainCount = SIZE(m_Chains);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const string &Label = m_Chains[ChainIndex]->m_Label;
		string Dom;
		SCOP40Bench::GetDomFromLabel(Label, Dom);
		m_DomToChainIndex[Dom] = ChainIndex;
		}

	const uint SeqCount = m_PairAlnDB.GetSeqCount();
	asserta(SeqCount%2 == 0);
	const uint PairCount = SeqCount/2;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Columns");
		const string &QLabel = m_PairAlnDB.GetLabel(2*PairIndex);
		const string &RLabel = m_PairAlnDB.GetLabel(2*PairIndex+1);
		string QDom;
		string RDom;
		SCOP40Bench::GetDomFromLabel(QLabel, QDom);
		SCOP40Bench::GetDomFromLabel(RLabel, RDom);
		uint QChainIndex = m_DomToChainIndex[QDom];
		uint RChainIndex = m_DomToChainIndex[RDom];
		m_ChainIdxsQ.push_back(QChainIndex);
		m_ChainIdxsR.push_back(RChainIndex);

		const PDBChain &QChain = *m_Chains[QChainIndex];
		const PDBChain &RChain = *m_Chains[RChainIndex];
		uint QL = QChain.GetSeqLength();
		uint RL = RChain.GetSeqLength();
		const string &QRow = m_PairAlnDB.GetSeq(2*PairIndex);
		const string &RRow = m_PairAlnDB.GetSeq(2*PairIndex+1);
		const uint ColCount = SIZE(QRow);
		asserta(SIZE(QRow) == SIZE(RRow));

		vector<uint> PosQs;
		vector<uint> PosRs;
		uint PosQ = 0;
		uint PosR = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char q = QRow[Col];
			char r = RRow[Col];
			if (!isgap(q) && !isgap(r))
				{
				PosQs.push_back(PosQ);
				PosRs.push_back(PosR);
				}
			if (!isgap(q))
				++PosQ;
			if (!isgap(r))
				++PosR;
			}
		m_PosVecsQ.push_back(PosQs);
		m_PosVecsR.push_back(PosRs);
		}
	}

const PDBChain &Trainer::GetChain(uint ChainIdx) const
	{
	asserta(ChainIdx < SIZE(m_Chains));
	return *m_Chains[ChainIdx];
	}

uint Trainer::GetPairCount() const
	{
	uint PairCount = SIZE(m_ChainIdxsQ);
	asserta(SIZE(m_ChainIdxsR) == PairCount);
	asserta(SIZE(m_PosVecsQ) == PairCount);
	asserta(SIZE(m_PosVecsR) == PairCount);
	return PairCount;
	}

void Trainer::Scan(TRAINER_ONPAIR OnPair,
  TRAINER_ONCOL OnCol) const
	{
	const uint PairCount = GetPairCount();
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Scanning");
		uint ChainIdxQ = m_ChainIdxsQ[PairIndex];
		uint ChainIdxR = m_ChainIdxsR[PairIndex];
		const vector<uint> &PosQs = m_PosVecsQ[PairIndex];
		const vector<uint> &PosRs = m_PosVecsR[PairIndex];
		const uint ColCount = SIZE(PosQs);
		asserta(SIZE(PosRs) == ColCount);
		if (OnPair != 0)
			OnPair(*this, ChainIdxQ, ChainIdxR, PosQs, PosRs);
		if (OnCol != 0)
			{
			for (uint Col = 0; Col < ColCount; ++Col)
				OnCol(*this, PosQs[Col], PosRs[Col]);
			}
		}
	}

void Trainer::TrainLogOdds(
  uint AlphaSize,
  TRAINER_ONPAIR OnPair,
  TRAINER_ALPHACOL AlphaCol,
  LogOdds &LO) const
	{
	LO.Init(AlphaSize);
	const uint PairCount = GetPairCount();
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Training");
		uint ChainIdxQ = m_ChainIdxsQ[PairIndex];
		uint ChainIdxR = m_ChainIdxsR[PairIndex];
		const vector<uint> &PosQs = m_PosVecsQ[PairIndex];
		const vector<uint> &PosRs = m_PosVecsR[PairIndex];
		const uint ColCount = SIZE(PosQs);
		asserta(SIZE(PosRs) == ColCount);
		OnPair(*this, ChainIdxQ, ChainIdxR, PosQs, PosRs);
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			uint LetterQ, LetterR;
			AlphaCol(*this, PosQs[Col], PosRs[Col],
			  LetterQ, LetterR);
			if (LetterQ < AlphaSize && LetterR < AlphaSize)
				LO.AddTruePair(LetterQ, LetterR);
			}
		}
	}
