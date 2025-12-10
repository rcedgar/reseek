#include "myutils.h"
#include "fragaligner.h"
#include "sort.h"

static uint SSCharToInt(char c)
	{
	switch (c)
		{
	case 'h': return 0;
	case 's': return 1;
	case 't': return 2;
	case '~': return 3;
		}
	asserta(false);
	return UINT_MAX;
	}

static void GetSSCounts(string &SS, uint Counts[4])
	{
	memset(Counts, 0, sizeof(Counts));
	const uint L = SIZE(SS);
	for (uint i = 0; i < L; ++i)
		Counts[SSCharToInt(SS[i])] += 1;
	}

static bool SelectCandidateS3E(const string &SS3, uint M, uint &Lo)
	{
	uint L = SIZE(SS3);
	if (M < 2*L)
		return false;
	Lo = randu32()%(L - M);
	return true;
	}

static uint s_M;
static vector<PDBChain *> *s_ptrFrags;
static vector<PDBChain *> *s_ptrChains;
static uint s_CandidateCount;
static mutex s_Lock;
static vector<float> *s_ptrTotalScores;
static float s_MinScore = 10;
static uint s_NextChainIdx;

static void ThreadBody(uint ThreadIndex)
	{
	vector<float> TotalScores(s_CandidateCount);
	vector<ChainData *> FragCDs(s_CandidateCount);
	for (uint CandidateIdx = 0; CandidateIdx < s_CandidateCount; ++CandidateIdx)
		{
		ChainData &CDF = *new ChainData;
		CDF.SetChain(*(*s_ptrFrags)[CandidateIdx]);
		FragCDs[CandidateIdx] = &CDF;
		}

	FragAligner FA;
	const uint ChainCount = SIZE(*s_ptrChains);
	for (;;)
		{
		s_Lock.lock();
		uint MyChainIdx = s_NextChainIdx++;
		s_Lock.unlock();
		if (MyChainIdx >= ChainCount)
			break;
		ProgressStep(MyChainIdx, ChainCount, "Aligning");
		PDBChain &Q = *(*s_ptrChains)[MyChainIdx];

		ChainData CDQ;
		CDQ.SetChain(Q);
		const uint LQ = Q.GetSeqLength();

		for (uint CandidateIdx = 0; CandidateIdx < s_CandidateCount; ++CandidateIdx)
			{
			ChainData &CDF = *FragCDs[CandidateIdx];
			for (uint PosQ = 0; PosQ + s_M < LQ; ++PosQ)
				{
				FA.Align(CDQ, CDF, PosQ);
				float Score = FA.m_DALIScore;
				asserta(!isnan(Score));
				asserta(!isinf(Score));
				if (Score >= s_MinScore)
					TotalScores[CandidateIdx] += Score;
				}
			float x = TotalScores[CandidateIdx];
			asserta(!isnan(x));
			asserta(!isinf(x));
			}
		}

	s_Lock.lock();
	for (uint CandidateIdx = 0; CandidateIdx < s_CandidateCount; ++CandidateIdx)
		(*s_ptrTotalScores)[CandidateIdx] += TotalScores[CandidateIdx];
	s_Lock.unlock();
	}

void cmd_build_s3e_library()
	{
	const string &ChainsFN = g_Arg1;
	asserta(optset_alpha_size);
	asserta(optset_subsample);
	asserta(optset_m);
	const uint AS = opt(alpha_size);
	s_CandidateCount = opt(subsample);
	s_M = opt(m);
	s_MinScore = 10;
	if (optset_minscore)
		s_MinScore = (float) opt(minscore);

	vector<PDBChain *> Chains;
	Progress("Reading chains...");
	ReadChains(ChainsFN, Chains);
	const uint ChainCount = SIZE(Chains);
	Progress(" %u done.\n", ChainCount);
	s_ptrChains = &Chains;

	vector<string> SSs(ChainCount);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		Chains[ChainIdx]->GetSS(SSs[ChainIdx]);

	vector<PDBChain *> Frags(s_CandidateCount);
	vector<string> FragSSs(s_CandidateCount);
	s_ptrFrags = &Frags;
	for (uint CandidateIdx = 0; CandidateIdx < s_CandidateCount; ++CandidateIdx)
		{
		ProgressStep(CandidateIdx, s_CandidateCount, "Candidates");
		uint Lo;
		uint TryCount = 0;
		uint ChainIdx = UINT_MAX;
		for (;;)
			{
			ChainIdx = randu32()%ChainCount;
			const PDBChain &Chain = *Chains[ChainIdx];
			uint L = Chain.GetSeqLength();
			if (L < 2*s_M)
				{
				if (++TryCount > 100)
					Die("tries");
				continue;
				}
			Lo = randu32()%(L - s_M);
			break;
			}
		PDBChain *Frag = new PDBChain;
		Chains[ChainIdx]->GetSubChain(*Frag, Lo, s_M);
		Frags[CandidateIdx] = Frag;
		string &FragSS = FragSSs[CandidateIdx];
		const string &SS = SSs[ChainIdx];
		for (uint PosF = 0; PosF < s_M; ++PosF)
			FragSS.push_back(SS[Lo+PosF]);
		}

	vector<float> TotalScores(s_CandidateCount);
	s_ptrTotalScores = &TotalScores;

	vector<thread *> ts;
	uint ThreadCount = GetRequestedThreadCount();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	vector<uint> Order(s_CandidateCount);
	QuickSortOrderDesc(TotalScores.data(), s_CandidateCount, Order.data());

	for (uint k = 0; k < s_CandidateCount; ++k)
		{
		uint CandidateIdx = Order[k];
		const PDBChain &Frag = *Frags[CandidateIdx];
		const string &FragSS = FragSSs[CandidateIdx];
		Log("%10.3g  %s  %s\n",
			TotalScores[CandidateIdx],
			FragSS.c_str(),
			Frag.m_Label.c_str());
		}
	}
