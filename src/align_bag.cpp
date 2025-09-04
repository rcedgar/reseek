#include "myutils.h"
#include "dssaligner.h"

extern float GetSelfRevScore(DSSAligner &DA, DSS &D, const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers);


static ChainBag *MakeBag(
	const DSSParams &Params,
	MuKmerFilter &MKF,
	DSSAligner &DASelfRev,
	DSS &D,
	const PDBChain &QChain)
	{
	D.Init(QChain);

	vector<vector<byte> > *ptrQProfile = new vector<vector<byte> >;
	vector<byte> *ptrQMuLetters = new vector<byte>;
	vector<uint> *ptrQMuKmers = new vector<uint>;

	D.GetProfile(*ptrQProfile);
	D.GetMuLetters(*ptrQMuLetters);
	D.GetMuKmers(*ptrQMuLetters, *ptrQMuKmers, Params.m_MKFPatternStr);

	float QSelfRevScore = 
		GetSelfRevScore(DASelfRev, D, QChain,
				*ptrQProfile, ptrQMuLetters, ptrQMuKmers);

	uint16_t *HT = MKF.CreateEmptyHashTable();
	MKF.SetHashTable(*ptrQMuKmers, HT);

	ChainBag *ptrCBQ = new ChainBag;
	ptrCBQ->m_ptrChain = &QChain;
	ptrCBQ->m_ptrProfile = ptrQProfile;
	ptrCBQ->m_ptrMuLetters = ptrQMuLetters;
	ptrCBQ->m_ptrMuKmers = ptrQMuKmers;
	ptrCBQ->m_SelfRevScore = QSelfRevScore;
	ptrCBQ->m_ptrProfPara = DASelfRev.m_ProfPara;
	ptrCBQ->m_ptrProfParaRev = DASelfRev.m_ProfParaRev;
	ptrCBQ->m_ptrKmerHashTableQ = HT;

	return ptrCBQ;
	}

void cmd_align_bag()
	{
	if (!optset_input2)
		Die("Must specify -input2");

	const string &QFN = g_Arg1;
	const string &TFN = opt(input2);

	vector<PDBChain *> ChainsQ;
	vector<PDBChain *> ChainsT;
	ReadChains(QFN, ChainsQ);
	ReadChains(TFN, ChainsT);

	optset_sensitive = true;
	opt(sensitive) = true;
	DSSParams Params;
	Params.SetDSSParams(DM_AlwaysSensitive);
	Params.m_UsePara = false;
	Params.m_Omega = 0;

	const uint ChainCountQ = SIZE(ChainsQ);
	const uint ChainCountT = SIZE(ChainsT);
	asserta(ChainCountQ == 1);
	asserta(ChainCountT == 1);

	const PDBChain &ChainQ = *ChainsQ[0];
	const PDBChain &ChainT = *ChainsT[0];

	DSS D;
	MuKmerFilter MKF;
	DSSAligner DA;
	DSSAligner DASelfRev;
	D.SetParams(Params);
	MKF.SetParams(Params);
	DA.SetParams(Params);
	DASelfRev.SetParams(Params);

	ChainBag &BagA = *MakeBag(Params, MKF, DASelfRev, D, ChainQ);
	ChainBag &BagB = *MakeBag(Params, MKF, DASelfRev, D, ChainT);

	DA.AlignBags(BagA, BagB);
	DA.ToAln(g_fLog, true);
	}
