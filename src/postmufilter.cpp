#include "myutils.h"
#include "dssaligner.h"
#include "chainreader2.h"

#define JOIN	0

#if JOIN
struct JoinData
	{
	uint idx1 = UINT_MAX;
	uint idx2 = UINT_MAX;
	double E;
	};
static map<pair<string, string>, JoinData> s_JoinMap;
static atomic<uint> s_nAB;
static atomic<uint> s_nBA;

//               0                 1         2      3      4
//d1w6ga1/b.30.2.1  d2oqea1/b.30.2.1  5.07e-20   1099   9211
//d2oqea1/b.30.2.1  d1w6ga1/b.30.2.1  5.07e-20   9211   1099
//d1r1ha_/d.92.1.4  d3dwba_/d.92.1.0  5.95e-20   2392  10429
//d3dwba_/d.92.1.0  d1r1ha_/d.92.1.4  5.95e-20  10429   2392
// d1gwea_/e.5.1.1   d1m7sa_/e.5.1.1  1.57e-19   4645  11069
// d1m7sa_/e.5.1.1   d1gwea_/e.5.1.1  1.57e-19  11069   4645

static void ReadJoin()
	{
	Progress("Reading join...");
	FILE *f = OpenStdioFile("join.tsv");
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 5);
		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];

		JoinData JD;
		JD.E = StrToFloat(Fields[2]);
		JD.idx1 = StrToUint(Fields[3]);
		JD.idx2 = StrToUint(Fields[4]);

		pair<string, string> Labels(Label1, Label2);
		s_JoinMap[Labels] = JD;
		}
	Progress(" done.\n");
	}
#endif

float GetSelfRevScore(DSSAligner &DA, DSS &D, const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers);

// Query & DB need C-alpha
void cmd_postmufilter()
	{
#if JOIN
	ReadJoin();
	FILE *fJoin2 = CreateStdioFile("join2.tsv");
#endif
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
	D.SetParams(Params);

	DSSAligner DA;
	DA.SetParams(Params);

	vector<DSSAligner *> DAs;
	DAs.reserve(QueryCount);
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
		float QSelfRevScore = GetSelfRevScore(DA, D, QChain, *ptrQProfile, ptrQMuLetters, ptrQMuKmers);

		ptrDA->SetQuery(QChain, ptrQProfile, ptrQMuLetters, ptrQMuKmers, QSelfRevScore);
		DAs.push_back(ptrDA);
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
#if JOIN
		ProgressStep(LineIdx, LineCount, "Scanning AB %u, BA %u", s_nAB.load(), s_nBA.load());
#else
		ProgressStep(LineIdx, LineCount, "Scanning");
#endif
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

		DA.SetTarget(DBChain, &DBProfile, &DBMuLetters, &DBMuKmers, DBSelfRevScore);

		const uint FilHitCount = StrToUint(Fields[1]);
		asserta(FilHitCount + 2 == FieldCount);
		for (uint FilHitIdx = 0; FilHitIdx < FilHitCount; ++FilHitIdx)
			{
			uint QueryIdx = StrToUint(Fields[FilHitIdx+2]);
			asserta(QueryIdx < QueryCount);

			DSSAligner &DA = *DAs[QueryIdx];
			DA.SetTarget(DBChain, &DBProfile, &DBMuLetters, &DBMuKmers, DBSelfRevScore);
			DA.AlignQueryTarget();
#if JOIN
			{
			DA.AlignQueryTarget_Trace();
			const string &LabelA = DA.m_ChainA->m_Label;
			const string &LabelB = DA.m_ChainB->m_Label;
			pair<string, string> AB(LabelA, LabelB);
			pair<string, string> BA(LabelB, LabelA);
			map<pair<string, string>, JoinData>::const_iterator iterAB = s_JoinMap.find(AB);
			map<pair<string, string>, JoinData>::const_iterator iterBA = s_JoinMap.find(BA);
			if (iterAB != s_JoinMap.end())
				{
				++s_nAB;
				const JoinData &JD = iterAB->second;
				fprintf(fJoin2, "AB");
				fprintf(fJoin2, "\t%s", LabelA.c_str());
				fprintf(fJoin2, "\t%s", LabelB.c_str());
				fprintf(fJoin2, "\t%.3g", DA.m_EvalueA);
				fprintf(fJoin2, "\t%.3g", JD.E);
				fprintf(fJoin2, "\t%u", JD.idx1);
				fprintf(fJoin2, "\t%u", JD.idx2);
				fprintf(fJoin2, "\n");
				}
			else if (iterBA != s_JoinMap.end())
				{
				++s_nBA;
				const JoinData &JD = iterBA->second;
				fprintf(fJoin2, "BA");
				fprintf(fJoin2, "\t%s", LabelA.c_str());
				fprintf(fJoin2, "\t%s", LabelB.c_str());
				fprintf(fJoin2, "\t%.3g", DA.m_EvalueA);
				fprintf(fJoin2, "\t%.3g", JD.E);
				fprintf(fJoin2, "\t%u", JD.idx1);
				fprintf(fJoin2, "\t%u", JD.idx2);
				fprintf(fJoin2, "\n");
				}
			DA.Align_NoAccel();
			Log("NoAccel E=%.3g\n", DA.m_EvalueA);
			}
#endif
			if (DA.m_EvalueA <= MaxEvalue)
				DA.ToTsv(fTsv, true);
			++ScannedCount;
			}
		}
	Ok = LR.ReadLine(Line);
	asserta(!Ok);
	LR.Close();
	CloseStdioFile(fTsv);
	}
