#include "myutils.h"
#include "dssaligner.h"
#include "chainreader2.h"

/***
[c300a5f] Add bags to postmufilter but not used, 
	SEPQ0.1=0.2109 SEPQ1=0.3142 SEPQ10=0.3878 S1FP=0.3347 N1FP=152204 area=7.14
***/

float GetSelfRevScore(DSSAligner &DA, DSS &D, const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers);

// Query & DB need C-alpha
void cmd_postmufilter()
	{
	asserta(optset_db);
	asserta(optset_filin);
	asserta(optset_dbsize);
	const string &QueryCAFN = g_Arg1;
	FILE *fTsv = CreateStdioFile(opt_output);
	float MaxEvalue = 10;
	if (optset_evalue)
		MaxEvalue = (float) opt_evalue;

	DSSParams Params;
	Params.SetFromCmdLine((uint) opt_dbsize);

	//DBSearcher DBS;
	//DBS.m_Params = &Params;
	//DBS.LoadDB(QueryCAFN);
	//const uint QueryCount = DBS.GetDBChainCount();

	vector<PDBChain *> QChains;
	ReadChains(QueryCAFN, QChains);
	const uint QueryCount = SIZE(QChains);

	DSS D;
	//MuKmerFilter MKF;
	DSSAligner DASelfRev;
	D.SetParams(Params);
	DASelfRev.SetParams(Params);
	//MKF.SetParams(Params);

	vector<DSSAligner *> DAs;
	vector<ChainBag *> ChainBagsQ;
	DAs.reserve(QueryCount);
	ChainBagsQ.reserve(QueryCount);
	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		{
		ProgressStep(QueryIdx, QueryCount, "Index query");
		const PDBChain &QChain = *QChains[QueryIdx];
		D.Init(QChain);

		DSSAligner *ptrDA = new DSSAligner;
		ptrDA->SetParams(Params);
		vector<vector<byte> > *ptrQProfile = new vector<vector<byte> >;
		vector<byte> *ptrQMuLetters = new vector<byte>;
		vector<uint> *ptrQMuKmers = new vector<uint>;

		D.GetProfile(*ptrQProfile);
		D.GetMuLetters(*ptrQMuLetters);
		D.GetMuKmers(*ptrQMuLetters, *ptrQMuKmers);
		float QSelfRevScore = GetSelfRevScore(DASelfRev, D, QChain, *ptrQProfile, ptrQMuLetters, ptrQMuKmers);

		//MKF.SetQ(QChain.m_Label, ptrQMuLetters, ptrQMuKmers);

		ptrDA->SetQuery(QChain, ptrQProfile, ptrQMuLetters, ptrQMuKmers, QSelfRevScore);
		DAs.push_back(ptrDA);

		ChainBag *ptrCBQ = new ChainBag;
		ptrCBQ->m_ptrChain = &QChain;
		ptrCBQ->m_ptrProfile = ptrQProfile;
		ptrCBQ->m_ptrMuLetters = ptrQMuLetters;
		ptrCBQ->m_ptrMuKmers = ptrQMuKmers;
		ptrCBQ->m_SelfRevScore = QSelfRevScore;
		ptrCBQ->m_ptrProfPara = ptrDA->m_ProfPara;
		ptrCBQ->m_ptrProfParaRev = ptrDA->m_ProfParaRev;
		ptrCBQ->m_ptrKmerHashTableQ = ptrDA->m_MKF.GetHashTableQ();
		ptrCBQ->Validate("CBQ");
		ChainBagsQ.push_back(ptrCBQ);
		}
	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		ChainBagsQ[QueryIdx]->Validate("CBQ1");

	BCAData DB;
	DB.Open(opt_db);

	LineReader2 LR;
	LR.Open(opt_filin);
	string Line;
	vector<string> Fields;
	uint PairCount = 0;
	uint ScannedCount = 0;
	time_t TimeLastProgress = 0;
	bool Ok = LR.ReadLine(Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	uint FieldCount = SIZE(Fields);
	asserta(FieldCount == 2);
	const uint LineCount  = StrToUint(Fields[1]);
	vector<vector<byte> > DBProfile;
	vector<byte> DBMuLetters;
	vector<uint> DBMuKmers;
	float SelfRevScore = 0;
	ChainBag CBT;
	DSSAligner TheDA;
	TheDA.SetParams(Params);
	for (uint LineIdx = 0; LineIdx < LineCount; ++LineIdx)
		{
		ProgressStep(LineIdx, LineCount, "Scanning");
		bool Ok = LR.ReadLine(Line);
		asserta(Ok);
		double Pct = LR.GetPctDone();
		Split(Line, Fields, '\t');
		FieldCount = SIZE(Fields);
		asserta(FieldCount > 2);
		const uint TargetIdx = StrToUint(Fields[0]);

		PDBChain DBChain;
		DB.ReadChain(TargetIdx, DBChain);

		D.Init(DBChain);
		D.GetProfile(DBProfile);
		D.GetMuLetters(DBMuLetters);
		D.GetMuKmers(DBMuLetters, DBMuKmers);

		float DBSelfRevScore = GetSelfRevScore(DASelfRev, D, DBChain, DBProfile,
										   &DBMuLetters, &DBMuKmers);
		DASelfRev.UnsetQuery();

		CBT.m_ptrChain = &DBChain;
		CBT.m_ptrProfile = &DBProfile;
		CBT.m_ptrMuLetters = &DBMuLetters;
		CBT.m_ptrMuKmers = &DBMuKmers;
		CBT.m_SelfRevScore = DBSelfRevScore;
		CBT.m_ptrProfPara = 0;
		CBT.m_ptrProfParaRev = 0;
		CBT.Validate("CBT");

		const uint FilHitCount = StrToUint(Fields[1]);
		asserta(FilHitCount + 2 == FieldCount);
		for (uint FilHitIdx = 0; FilHitIdx < FilHitCount; ++FilHitIdx)
			{
			uint QueryIdx = StrToUint(Fields[FilHitIdx+2]);
			asserta(QueryIdx < QueryCount);

			//DSSAligner &DA = *DAs[QueryIdx];
			//DA.SetTarget(DBChain, &DBProfile, &DBMuLetters, &DBMuKmers, DBSelfRevScore);
			//DA.AlignQueryTarget();
			//if (DA.m_EvalueA <= MaxEvalue)
			//	DA.ToTsv(fTsv, true);
			++ScannedCount;

			//float E1 = DA.m_EvalueA;
			asserta(QueryIdx < SIZE(ChainBagsQ));
			const ChainBag &CBQ = *ChainBagsQ[QueryIdx];
			CBQ.Validate("Q");
			CBT.Validate("T");
			TheDA.AlignBags(CBQ, CBT);
			float E2 = TheDA.m_EvalueA;
			//asserta(feq(E1, E2));
			if (TheDA.m_EvalueA <= MaxEvalue)
				TheDA.ToTsv(fTsv, true);
			}
		}
	Ok = LR.ReadLine(Line);
	asserta(!Ok);
	LR.Close();
	CloseStdioFile(fTsv);
	}
