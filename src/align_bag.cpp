#include "myutils.h"
#include "dssaligner.h"
#include "statsig.h"

static ChainBag *MakeBag(
	const DSSParams &Params,
	MuKmerFilter &MKF,
	DSSAligner &DA_selfrev,
	DSS &D,
	const PDBChain &QChain)
	{
	D.Init(QChain, Params);

	vector<vector<byte> > *ptrQProfile = new vector<vector<byte> >;
	vector<byte> *ptrQMuLetters = new vector<byte>;
	vector<uint> *ptrQMuKmers = new vector<uint>;

	D.GetProfile(*ptrQProfile);
	D.GetMuLetters(*ptrQMuLetters);
	D.GetMuKmers(*ptrQMuLetters, *ptrQMuKmers, Params.m_MKFPatternStr);

	float QSelfRevScore = 
		GetSelfRevScore(DA_selfrev, Params, QChain,
				*ptrQProfile, ptrQMuLetters, ptrQMuKmers);

	uint16_t *HT = MKF.CreateEmptyHashTable();
	MKF.SetHashTable(*ptrQMuKmers, HT);

	ChainBag *ptrCBQ = new ChainBag;
	ptrCBQ->m_ptrChain = &QChain;
	ptrCBQ->m_ptrProfile = ptrQProfile;
	ptrCBQ->m_ptrMuLetters = ptrQMuLetters;
	ptrCBQ->m_ptrMuKmers = ptrQMuKmers;
	ptrCBQ->m_SelfRevScore = QSelfRevScore;
	ptrCBQ->m_ptrProfPara = DA_selfrev.m_ProfPara;
	ptrCBQ->m_ptrProfParaRev = DA_selfrev.m_ProfParaRev;
	ptrCBQ->m_ptrKmerHashTableQ = HT;

	return ptrCBQ;
	}

// Align pair, exactly one chain in each input file.
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
	DSSAligner DA_selfrev;
	//D.SetParams(Params);
	MKF.SetParams(Params);
	DA.SetParams(Params);
	DA_selfrev.SetParams(Params);

	ChainBag &BagA = *MakeBag(Params, MKF, DA_selfrev, D, ChainQ);
	ChainBag &BagB = *MakeBag(Params, MKF, DA_selfrev, D, ChainT);

	DA.AlignBagsMKF(BagA, BagB);
	if (DA.m_Path == "")
		ProgressLog("No alignment found\n");
	else
		DA.ToAln(g_fLog, true);
	}

// Align all-vs-all, compare with full S-W
void cmd_align_bags()
	{
	const string &QFN = g_Arg1;
	FILE *f = CreateStdioFile(opt(output));

	vector<PDBChain *> Chains;
	ReadChains(QFN, Chains);

	optset_sensitive = true;
	opt(sensitive) = true;
	DSSParams Params;
	Params.SetDSSParams(DM_AlwaysSensitive);
	Params.m_UsePara = false;
	Params.m_Omega = 0;
	StatSig::InitSensitive(SCOP40_DBSIZE);

	const uint ChainCount = SIZE(Chains);

	DSS D;
	MuKmerFilter MKF;
	DSSAligner DA_sw;
	DSSAligner DA_bag;
	DSSAligner DA_selfrev;
	//D.SetParams(Params);
	MKF.SetParams(Params);
	DA_selfrev.SetParams(Params);
	DA_sw.SetParams(Params);
	DA_bag.SetParams(Params);

	uint PairCount = ChainCount + ChainCount*(ChainCount-1)/2;
	uint PairIndex = 0;
	for (uint ChainIndexA = 0; ChainIndexA < ChainCount; ++ChainIndexA)
		{
		const PDBChain &ChainA = *Chains[ChainIndexA];
		ChainBag &BagA = *MakeBag(Params, MKF, DA_selfrev, D, ChainA);

		vector<vector<byte> > ProfileA;
		D.Init(ChainA, Params);
		D.GetProfile(ProfileA);
		float SelfRevScoreA = GetSelfRevScore(DA_selfrev, Params, ChainA, ProfileA, 0, 0);
		DA_sw.SetQuery(ChainA, &ProfileA, 0, 0, SelfRevScoreA);

		for (uint ChainIndexB = ChainIndexA; ChainIndexB < ChainCount; ++ChainIndexB)
			{
			ProgressStep(PairIndex++, PairCount, "Aligning");
			const PDBChain &ChainB = *Chains[ChainIndexB];

			uint LA = ChainA.GetSeqLength();
			uint LB = ChainB.GetSeqLength();
			if (LA < 400 || LB < 400)
				continue;

			vector<vector<byte> > ProfileB;
			D.Init(ChainB, Params);
			D.GetProfile(ProfileB);
			float SelfRevScoreB = GetSelfRevScore(DA_selfrev, Params, ChainB, ProfileB, 0, 0);

			DA_sw.SetTarget(ChainB, &ProfileB, 0, 0, SelfRevScoreB);
			DA_sw.Align_NoAccel();
			float E_sw = DA_sw.GetEvalue(true);
			if (E_sw > 1)
				continue;

			ChainBag &BagB = *MakeBag(Params, MKF, DA_selfrev, D, ChainB);
			DA_bag.AlignBagsMKF(BagA, BagB);

			if (f == 0)
				continue;

			bool problem = false;
			bool b = (DA_bag.m_MKF.m_BestChainScore > 0);

			fprintf(f, "%s", ChainA.m_Label.c_str());
			fprintf(f, "\t%s", ChainB.m_Label.c_str());

			fprintf(f, "\t%.2e", DA_sw.GetEvalue(true));
			if (b)
				fprintf(f, "\t%.2e", DA_bag.GetEvalue(true));
			else
				{
				if (E_sw < 0.01)
					problem = true;
				fprintf(f, "\tPROBE");
				}

			float PctId_sw = DA_sw.GetPctId();
			float PctId_bag = DA_bag.GetPctId();
			fprintf(f, "\t%.1f", PctId_sw);
			if (b)
				{
				if (PctId_sw - PctId_bag > 5)
					problem = true;
				fprintf(f, "\t%.1f", DA_bag.GetPctId());
				}
			else
				fprintf(f, "\tnobag");
			if (problem)
				fprintf(f, "\tPROBLEM");
			fprintf(f, "\n");
			fflush(f);
			}
		}
	CloseStdioFile(f);
	}
