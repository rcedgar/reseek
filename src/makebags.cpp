#include "myutils.h"
#include "dssaligner.h"

float GetSelfRevScore(DSSAligner &DA, DSS &D, const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers);

static const vector<PDBChain *> *s_ptr_Chains;
static vector<ChainBag *> *s_ptr_Bags;
static atomic<uint> s_QueryIdx;

void MakeBag(const PDBChain &Chain, ChainBag &CB,
	DSS &D, DSSAligner &DASelfRev, MuKmerFilter &MKF)
	{
	D.Init(Chain);

	vector<vector<byte> > *ptrProfile = new vector<vector<byte> >;
	vector<byte> *ptrMuLetters = new vector<byte>;
	vector<uint> *ptrMuKmers = new vector<uint>;

	D.GetProfile(*ptrProfile);
	D.GetMuLetters(*ptrMuLetters);
	D.GetMuKmers(*ptrMuLetters, *ptrMuKmers, DSSParams::m_MKFPatternStr);

	float SelfRevScore =  GetSelfRevScore(DASelfRev, D, Chain,
		*ptrProfile, ptrMuLetters, ptrMuKmers);

	uint16_t *HT = MKF.CreateEmptyHashTable();
	MKF.SetHashTable(*ptrMuKmers, HT);

	CB.m_ptrChain = &Chain;
	CB.m_ptrProfile = ptrProfile;
	CB.m_ptrMuLetters = ptrMuLetters;
	CB.m_ptrMuKmers = ptrMuKmers;
	CB.m_SelfRevScore = SelfRevScore;
	CB.m_ptrProfPara8 = DASelfRev.m_ProfPara8;
	CB.m_ptrProfPara16 = DASelfRev.m_ProfPara16;
	CB.m_ptrProfParaRev8 = DASelfRev.m_ProfParaRev8;
	CB.m_ptrProfParaRev16 = DASelfRev.m_ProfParaRev16;
	CB.m_ptrKmerHashTableQ = HT;

// Transfer ownsership of parasail profiles to ChainBag
	DASelfRev.m_ProfPara8 = 0;
	DASelfRev.m_ProfPara16 = 0;
	DASelfRev.m_ProfParaRev8 = 0;
	DASelfRev.m_ProfParaRev16 = 0;
	}

static void ThreadBody(uint ThreadIndex)
	{
	DSS D;
	DSSAligner DASelfRev;
	MuKmerFilter MKF;

	vector<ChainBag *> &Bags = *s_ptr_Bags;
	const vector<PDBChain *> &Chains = *s_ptr_Chains;
	const uint ChainCount = SIZE(Chains);

	if (ThreadIndex == 0 && ChainCount > 1)
		ProgressStep(0, ChainCount, "Index chains");
	for (;;)
		{
		uint ChainIdx = s_QueryIdx++;
		if (ChainIdx >= ChainCount)
			{
			if (ThreadIndex == 0 && ChainCount > 1)
				ProgressStep(ChainCount-1, ChainCount, "Index chains");
			return;
			}
		if (ThreadIndex == 0)
			ProgressStep(ChainIdx, ChainCount, "Index chains");

		const PDBChain &Chain = *Chains[ChainIdx];
		ChainBag &CB = *new ChainBag;
		MakeBag(Chain, CB, D, DASelfRev, MKF);
		asserta(Bags[ChainIdx] == 0);
		Bags[ChainIdx] = &CB;
		}
	}

void MakeBags(const vector<PDBChain *> Chains, vector<ChainBag *> &Bags)
	{
	s_ptr_Chains = &Chains;
	s_ptr_Bags = &Bags;
	s_QueryIdx = 0;

	Bags.clear();
	Bags.resize(SIZE(Chains));

	uint ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	}