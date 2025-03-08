#include "myutils.h"
#include "dssaligner.h"
#include "chainreader2.h"

float GetSelfRevScore(DSSAligner &DA, DSS &D, const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers);

static mutex s_IndexQueryLock;
static mutex s_ScanLock;
static uint s_QueryIdx;
static uint s_QueryCount;
static vector<PDBChain *> *s_ptrQChains;
static vector<ChainBag *> *s_ptrChainBagsQ;
static const DSSParams *s_ptrParams;
static BCAData *s_ptrDB;
static uint s_LineIdx;
static uint s_LineCount;
static uint s_ScannedCount;
static LineReader2 *s_ptrLR;
static float s_MaxEvalue = 10;
static FILE *s_fTsv;
static FILE *s_fTsv2;

static void ThreadBody_IndexQuery(uint ThreadIndex)
	{
	const DSSParams &Params = *s_ptrParams;
	DSS D;
	DSSAligner DASelfRev;
	D.SetParams(Params);
	DASelfRev.SetParams(Params);
	MuKmerFilter MKF;
	MKF.SetParams(Params);
	vector<ChainBag *> &ChainBagsQ = *s_ptrChainBagsQ;
	vector<PDBChain *> &QChains = *s_ptrQChains;

	if (ThreadIndex == 0)
		ProgressStep(0, s_QueryCount, "Index query");
	for (;;)
		{
		s_IndexQueryLock.lock();
		uint QueryIdx = s_QueryIdx++;
		s_IndexQueryLock.unlock();

		if (QueryIdx >= s_QueryCount)
			{
			if (ThreadIndex == 0)
				ProgressStep(s_QueryCount-1, s_QueryCount, "Index query");
			return;
			}

		if (ThreadIndex == 0)
			ProgressStep(QueryIdx, s_QueryCount, "Index query");

		const PDBChain &QChain = *QChains[QueryIdx];
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
		asserta(ChainBagsQ[QueryIdx] == 0);
		s_IndexQueryLock.lock();
		ChainBagsQ[QueryIdx] = ptrCBQ;
		s_IndexQueryLock.unlock();

		//DASelfRev.m_MKF.ForceZero();
		DASelfRev.m_ProfPara = 0;
		DASelfRev.m_ProfParaRev = 0;
		}
	}

static void ThreadBody_Scan(uint ThreadIndex)
	{
	const DSSParams &Params = *s_ptrParams;
	DSS D;
	DSSAligner DASelfRev;
	D.SetParams(Params);
	DASelfRev.SetParams(Params);
	MuKmerFilter MKF;
	MKF.SetParams(Params);
	vector<ChainBag *> &ChainBagsQ = *s_ptrChainBagsQ;
	vector<PDBChain *> &QChains = *s_ptrQChains;
	const BCAData &DB = *s_ptrDB;
	vector<vector<byte> > DBProfile;
	vector<byte> DBMuLetters;
	vector<uint> DBMuKmers;
	float SelfRevScore = 0;
	ChainBag CBT;
	DSSAligner TheDA;
	TheDA.SetParams(Params);
	LineReader2 &LR = *s_ptrLR;

	string Line;
	vector<string> Fields;

	if (ThreadIndex == 0)
		ProgressStep(0, s_LineCount, "Scanning");
	for (;;)
		{
		s_ScanLock.lock();
		uint LineIdx = s_LineIdx++;
		bool Ok = LR.ReadLine(Line);
		s_ScanLock.unlock();
		if (LineIdx >= s_LineCount)
			{
			if (ThreadIndex == 0)
				ProgressStep(s_LineCount-1, s_LineCount, "Scanning");
			return;
			}
		asserta(Ok);
		if (ThreadIndex == 0 && LineIdx+1 < s_LineCount)
			ProgressStep(LineIdx, s_LineCount, "Scanning");
		double Pct = LR.GetPctDone();
		Split(Line, Fields, '\t');
		uint FieldCount = SIZE(Fields);
		asserta(FieldCount > 2);
		const uint TargetIdx = StrToUint(Fields[0]);

		PDBChain DBChain;
		DB.ReadChain(TargetIdx, DBChain);

		D.Init(DBChain);
		D.GetProfile(DBProfile);
		D.GetMuLetters(DBMuLetters);
		D.GetMuKmers(DBMuLetters, DBMuKmers, Params.m_MKFPatternStr);

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

		const uint FilHitCount = StrToUint(Fields[1]);
		asserta(FilHitCount + 2 == FieldCount);
		for (uint FilHitIdx = 0; FilHitIdx < FilHitCount; ++FilHitIdx)
			{
			uint QueryIdx = StrToUint(Fields[FilHitIdx+2]);
			asserta(QueryIdx < SIZE(ChainBagsQ));
			const ChainBag &CBQ = *ChainBagsQ[QueryIdx];
			TheDA.AlignBags(CBQ, CBT);
			if (TheDA.m_EvalueA <= s_MaxEvalue)
				TheDA.ToTsv(s_fTsv, true);
			s_ScanLock.lock();
			if (s_fTsv2)
				{
				fprintf(s_fTsv2, "%s\t%s\t%c\n",
						DBChain.m_Label.c_str(),
						CBQ.m_ptrChain->m_Label.c_str(),
						tof(TheDA.m_EvalueA <= s_MaxEvalue));
				}
			++s_ScannedCount;
			s_ScanLock.unlock();
			}
		}
	}

// Query & DB need C-alpha
void PostMuFilter(const DSSParams &Params,
				  const string &MuFilterTsvFN,
				  const string &QueryCAFN,
				  const string &DBBCAFN,
				  float MaxEvalue,
				  const string &HitsFN)
	{
	time_t t0 = time(0);
	LineReader2 LR;
	LR.Open(MuFilterTsvFN);
	s_ptrLR = &LR;
	string Line;
	vector<string> Fields;
	bool Ok = LR.ReadLine(Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	uint FieldCount = SIZE(Fields);
	asserta(FieldCount == 2);
	s_LineCount  = StrToUint(Fields[1]);
	if (s_LineCount == 0)
		{
		Warning("No hits found by mufilter pass");
		return;
		}

	s_fTsv = CreateStdioFile(HitsFN);
	s_fTsv2 = CreateStdioFile(opt(output2));
	s_MaxEvalue = MaxEvalue;

	vector<PDBChain *> QChains;
	ReadChains(QueryCAFN, QChains);
	s_QueryCount = SIZE(QChains);

	DSS D;
	DSSAligner DASelfRev;
	D.SetParams(Params);
	DASelfRev.SetParams(Params);

	vector<ChainBag *> ChainBagsQ;
	ChainBagsQ.resize(s_QueryCount, 0);

	s_ptrQChains = &QChains;
	s_ptrChainBagsQ = &ChainBagsQ;
	s_ptrParams = &Params;

	uint ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody_IndexQuery, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	asserta(SIZE(*s_ptrChainBagsQ) == s_QueryCount);

	BCAData DB;
	DB.Open(DBBCAFN);
	s_ptrDB = &DB;

	ts.clear();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody_Scan, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	Ok = LR.ReadLine(Line);
	asserta(!Ok);
	LR.Close();
	CloseStdioFile(s_fTsv);
	CloseStdioFile(s_fTsv2);
	time_t t1 = time(0);
	Progress("Alignments %u secs\n", uint(t1 - t0));
	}

void cmd_postmuprefilter()
	{
	const string &QueryCAFN = g_Arg1;
	const string &HitsFN = opt(output);

	asserta(optset_db);
	const string &DBCAFN = opt(db);

	asserta(optset_filin);
	const string &MuFilterTsvFN = opt(filin);

	s_fTsv = CreateStdioFile(opt(output));
	float MaxEvalue = 10;
	if (optset_evalue)
		MaxEvalue = (float) opt(evalue);

	DSSParams Params;
	asserta(optset_dbsize);
	Params.SetDSSParams(DM_AlwaysSensitive, SCOP40_DBSIZE);

	PostMuFilter(Params,
				 opt(filin),
				 QueryCAFN,
				 DBCAFN,
				 MaxEvalue,
				 HitsFN);
	}
