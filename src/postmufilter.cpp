#include "myutils.h"
#include "dbsearcher.h"
#include "chainreader2.h"

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

	DSSParams Params;
	Params.SetFromCmdLine((uint) opt_dbsize);

	DBSearcher DBS;
	DBS.m_Params = &Params;
	DBS.LoadDB(QueryCAFN);
	const uint QueryCount = DBS.GetDBChainCount();

	DSS D;
	D.SetParams(Params);

	DSSAligner DA;
	DA.SetParams(Params);

	vector<DSSAligner *> DAs;
	DAs.reserve(QueryCount);
	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		{
		ProgressStep(QueryIdx, QueryCount, "Index query");
		DSSAligner *DA = new DSSAligner;
		DA->SetParams(Params);
		PDBChain &QChain = *DBS.m_DBChains[QueryIdx];
		const vector<vector<byte> > *ptrQProfile = DBS.m_DBProfiles[QueryIdx];
		const vector<byte> *ptrQMuLetters = DBS.m_DBMuLettersVec[QueryIdx];
		const vector<uint> *ptrQMuKmers = DBS.m_DBMuKmersVec[QueryIdx];
		float SelfRevScore = DBS.m_DBSelfRevScores[QueryIdx];
		DA->SetQuery(QChain, ptrQProfile, ptrQMuLetters, ptrQMuKmers, SelfRevScore);
		DAs.push_back(DA);
		}

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

		float DBSelfRevScore = GetSelfRevScore(DA, D, DBChain, DBProfile,
											   &DBMuLetters, &DBMuKmers);

		//DA.SetTarget(DBChain, &DBProfile, &DBMuLetters, &DBMuKmers, DBSelfRevScore);

		const uint FilHitCount = StrToUint(Fields[1]);
		asserta(FilHitCount + 2 == FieldCount);
		for (uint FilHitIdx = 0; FilHitIdx < FilHitCount; ++FilHitIdx)
			{
			uint QueryIdx = StrToUint(Fields[FilHitIdx+2]);
			asserta(QueryIdx < QueryCount);

			DSSAligner &DA = *DAs[QueryIdx];
			DA.SetTarget(DBChain, &DBProfile, &DBMuLetters, &DBMuKmers, DBSelfRevScore);
			DA.AlignQueryTarget();
			++ScannedCount;
			}
		}
	Ok = LR.ReadLine(Line);
	asserta(!Ok);
	LR.Close();
	}
