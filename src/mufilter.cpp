#include "myutils.h"
#include "seqdb.h"
#include "mukmerfilter.h"
#include "dss.h"
#include "alpha.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "fastaseqsource.h"

static atomic<uint> s_PairCount;
static atomic<uint> s_BCSCount;
static uint s_TargetCount;
static mutex s_Lock;
static time_t s_last_progress;

static void GetMuKmers(const string &PatternStr,
					   const vector<byte> &Letters,
					   vector<uint> &Kmers)
	{
	Kmers.clear();
	const uint PatternLength = SIZE(PatternStr);
	const uint L = SIZE(Letters);
	Kmers.reserve(L);
	for (uint Pos = 0; Pos + PatternLength <= L; ++Pos)
		{
		uint Kmer = 0;
		for (uint j = 0; j < PatternLength; ++j)
			{
			if (PatternStr[j] == '1')
				{
				assert(Pos + j < SIZE(Letters));
				byte Letter = Letters[Pos + j];
				assert(Letter < 36);
				Kmer = Kmer*36 + Letter;
				}
			}
		assert(Kmer != UINT_MAX);
		Kmers.push_back(Kmer);
		}
	}

static void ThreadBody(uint ThreadIndex, 
					   const DSSParams *ptrParams,
					   FASTASeqSource *ptrFSS, 
					   vector<MuKmerFilter *> *ptrMKFs)
	{
	const DSSParams &Params = *ptrParams;
	FASTASeqSource &FSS = *ptrFSS;
	vector<MuKmerFilter *> &MKFs = *ptrMKFs;
	const string &PatternStr = Params.m_PatternStr;
	const uint QueryCount = SIZE(MKFs);

	ObjMgr OM;

	SeqInfo *SI = OM.GetSeqInfo();
	vector<byte> MuLettersT;
	vector<uint> MuKmersT;
	for (;;)
		{
		bool Ok = FSS.GetNext(SI);
		if (!Ok)
			return;
		uint TL = SI->m_L;
		const byte *T = SI->m_Seq;
		MuLettersT.clear();
		MuLettersT.reserve(TL);
		for (uint i = 0; i < TL; ++i)
			{
			byte c = T[i];
			byte Letter = g_CharToLetterMu[c];
			asserta(Letter < 36);
			MuLettersT.push_back(Letter);
			}
		GetMuKmers(PatternStr, MuLettersT, MuKmersT);
		bool DoProgress = false;
		s_Lock.lock();
		++s_TargetCount;
		if (s_TargetCount%100 == 0)
			{
			time_t now = time(0);
			if (now != s_last_progress)
				{
				DoProgress = true;
				s_last_progress = now;
				}
			}
		s_Lock.unlock();
		if (DoProgress)
			Progress("%u chains scanned   \r", s_TargetCount);

		for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
			{
			MuKmerFilter &MKF = *MKFs[QueryIdx];
			MKF.Align(MuLettersT, MuKmersT);
			++s_PairCount;
			if (MKF.m_BestChainScore > 0)
				{
				++s_BCSCount;
				Log("%d\n", MKF.m_BestChainScore);
				}
			}
		}
	}

void cmd_mufilter()
	{
	asserta(optset_db);
	const string &QueryFN = g_Arg1;
	const string &DBFN = string(opt_db);
	FASTASeqSource FSS;
	FSS.Open(DBFN);

	DSSParams Params;
	Params.SetFromCmdLine(10000);
	const string &PatternStr = Params.m_PatternStr;

	SeqDB QueryDB;
	QueryDB.FromFasta(QueryFN);
	const uint QueryCount = QueryDB.GetSeqCount();
	vector<MuKmerFilter *> MKFs;
	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		{
		ProgressStep(QueryIdx, QueryCount, "Indexing query");
		MuKmerFilter *MKF = new MuKmerFilter;
		MKF->SetParams(Params);
		const uint QL = QueryDB.GetSeqLength(QueryIdx);
		const string &Q = QueryDB.GetSeq(QueryIdx);
		vector<byte> *ptrMuLetters = new vector<byte>;
		vector<uint> *ptrMuKmers = new vector<uint>;
		ptrMuLetters->reserve(QL);
		for (uint i = 0; i < QL; ++i)
			{
			byte c = Q[i];
			byte Letter = g_CharToLetterMu[c];
			asserta(Letter < 36);
			ptrMuLetters->push_back(Letter);
			}
		GetMuKmers(PatternStr, *ptrMuLetters, *ptrMuKmers);
		MKF->ResetQ();
		MKF->SetQ(ptrMuLetters, ptrMuKmers);
		MKFs.push_back(MKF);
		}

	const uint ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody, ThreadIndex, &Params, &FSS, &MKFs);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	ProgressLog("%u / %u with +ve HSP score (%.1f%%)\n",
				s_BCSCount.load(), s_PairCount.load(),
				GetPct(s_BCSCount, s_PairCount));
	}
