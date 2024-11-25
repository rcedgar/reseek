#include "myutils.h"
#include "mermx.h"
#include "quarts.h"

/***
PairScoreDist()
Mu 5mer pair scores
Self: N=60466176	Min=20  LoQ=43  Med=47  HiQ=51  Max=75 Avg=18    StdDev=24.1
Rand: N=60466176	Min=-94 LoQ=-36 Med=-25 HiQ=-15 Max=56 Avg=-15.3 StdDev=16.6
***/

static void TestNbr(MerMx &MM, const string &sKmer, uint FSn, uint FSTicks)
	{
	uint Kmer = MM.StrToKmer(sKmer);
	uint DictSize = myipow(20, 6);

	uint *Kmers = myalloc(uint, DictSize);
	uint *Kmers_Brute = myalloc(uint, DictSize);

	uint n_Brute = MM.GetHighScoring6mers_Brute(Kmer, 78, Kmers_Brute, true);

	//short *Work = myalloc(short, 2*MM.m_AS3);

	auto c0 = std::chrono::high_resolution_clock::now();
	uint n = MM.GetHighScoring6mers(Kmer, 78, Kmers);
	auto c1 = std::chrono::high_resolution_clock::now();
	auto elapsed = c1 - c0;
	uint32_t cticks = (uint32_t) elapsed.count();
	ProgressLog("%s cticks=%u,%u, n=%u,%u\n",
				sKmer.c_str(), cticks, FSTicks, n, FSn);
	asserta(n == n_Brute);
	}

static void PairScoreDist()
	{
	MerMx MM;
	extern const short * const *Mu_S_ij_short;
	MM.Init(Mu_S_ij_short, 5, 36, 2);

	vector<float> Scores;
	uint DictSize = MM.m_AS_pow[5];
	Scores.reserve(DictSize);
	for (uint kmer = 0; kmer < DictSize; ++kmer)
		{
		ProgressStep(kmer, DictSize, "self scores");
		int self_score = MM.GetScoreKmerPair(kmer, kmer);
		Scores.push_back(float(self_score));
		}
	QuartsFloat QF;
	GetQuartsFloat(Scores, QF);
	ProgressLog("Self: ");
	QF.ProgressLogMe();

	Scores.clear();
	Scores.reserve(DictSize);
	for (uint i = 0; i < DictSize; ++i)
		{
		ProgressStep(i, DictSize, "pairscores");
		uint kmer1 = randu32()%DictSize;
		uint kmer2 = randu32()%DictSize;
		int self_score = MM.GetScoreKmerPair(kmer1, kmer2);
		Scores.push_back(float(self_score));
		}
	GetQuartsFloat(Scores, QF);
	ProgressLog("Rand: ");
	QF.ProgressLogMe();
	}

static void Test(const MerMx &MM, const char *KmerStr, int ScoreDelta)
	{
	assert(MM.m_AS_pow[5] == 60466176);
	asserta(strlen(KmerStr) == 5);
	uint *Kmers = myalloc(uint, 60466176);
	uint Kmer = MM.StrToKmer(KmerStr);
	assert(Kmer < MM.m_AS_pow[5]);
	int SelfScore = MM.GetScoreKmerPair(Kmer, Kmer);
	int MinScore = SelfScore - ScoreDelta;
	ProgressLog("Self score %s = %d, MinScore = %d\n", KmerStr, SelfScore, MinScore);
	uint n_Brute = MM.GetHighScoring5mers_Brute(Kmer, MinScore, Kmers, false);
	ProgressLog("n_Brute=%u\n", n_Brute);
	uint n_Fast = MM.GetHighScoring5mers(Kmer, MinScore, Kmers);
	ProgressLog("n_Fast=%u\n", n_Brute);
	myfree(Kmers);
	}

void cmd_testmu5mers()
	{
	MerMx MM;
	extern const short * const *Mu_S_ij_short;
	MM.Init(Mu_S_ij_short, 5, 36, 2);

	Test(MM, "AAAAA", 4);
	Test(MM, "ABCDE", 4);
	Test(MM, "abcde", 4);
	}
