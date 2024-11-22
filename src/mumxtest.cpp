#include "myutils.h"
#include "alpha.h"
#include "quarts.h"
#include "mermx.h"
#include "mumx.h"
#include <chrono>

/***
k=3 self N=46656	Min=12	LoQ=25	Med=28	HiQ=32	Max=45 Avg=28.4 StdDev=4.54
k=3 pair N=1000000	Min=-60	LoQ=-23 Med=-15 HiQ=-7	Max=38 Avg=-15.1 StdDev=12
k=4 self N=1679616	Min=16	LoQ=34	Med=38	HiQ=41	Max=60 Avg=37.8 StdDev=5.24
k=4 pair N=1000000	Min=-73	LoQ=-30	Med=-20 HiQ=-11 Max=40 Avg=-20.1 StdDev=13.8
k=5 self N=60466176	Min=20	LoQ=43	Med=47	HiQ=51	Max=75 Avg=18 StdDev=30.1
k=5 pair N=1000000	Min=-84 LoQ=-36	Med=-25 HiQ=-15	Max=47 Avg=-25.2 StdDev=15.4
k=6 self N=1000000	Min=30	LoQ=52	Med=57	HiQ=61	Max=88 Avg=57 StdDev=6.41
k=6 pair N=1000000	Min=-103 LoQ=-42 Med=-30 HiQ=-19 Max=47 Avg=-30.2 StdDev=16.9
***/

static uint s_dict_size;

static char *kmer2str(uint k, uint kmer)
	{
	static char s[7];
	for (uint i = 0; i < k; ++i)
		{
		s[i] = g_LetterToCharAmino[kmer%36];
		kmer /= 36;
		}
	return s;
	}

static uint str2kmer(const string &s)
	{
	uint k = SIZE(s);
	uint kmer = 0;
	for (uint i = 0; i < k; ++i)
		{
		byte c = s[i];
		uint code = g_CharToLetterMu[c];
		kmer = kmer*36 + code;
		}
	return kmer;
	}

static int get_kmer_self_score(uint k, uint kmer)
	{
	int sum = 0;
	for (uint i = 0; i < k; ++i)
		{
		uint code = kmer%36;
		sum += Mu_S_ij_short[code][code];
		kmer /= 36;
		}
	return sum;
	}

static int get_kmer_pair_score(uint k, uint kmer1, uint kmer2)
	{
	int sum = 0;
	for (uint i = 0; i < k; ++i)
		{
		uint code1 = kmer1%36;
		uint code2 = kmer2%36;
		sum += Mu_S_ij_short[code1][code2];
		kmer1 /= 36;
		kmer2 /= 36;
		}
	return sum;
	}

static void test_self_score(const string &s, int should_be)
	{
	uint k = SIZE(s);
	asserta(k == 6);
	uint kmer = str2kmer(s);
	int sc = get_kmer_self_score(k, kmer);
	ProgressLog("Self: %s = %3d (%3d)\n", s.c_str(), sc, should_be);
	}

static void test_pair_score(const string &s, int should_be,
							bool Trace = false)
	{
	uint k = SIZE(s);
	asserta(k == 6);
	uint kmer = str2kmer(s);
	int sc = get_kmer_self_score(k, kmer);
	uint nge78 = 0;
	for (uint kmer2 = 0; kmer2 < s_dict_size; ++kmer2)
		{
		int sc = get_kmer_pair_score(k, kmer, kmer2);
		if (sc >= 78)
			{
			if (Trace)
				Log(" %s = %d\n", kmer2str(k, kmer2), sc);
			++nge78;
			}
		}
	ProgressLog("Nbrs: %s = %5d (%5d)\n", s.c_str(), nge78, should_be);
	}

static void Info()
	{
	int scmin = Mu_S_ij_short[0][0];
	int scmax = Mu_S_ij_short[0][0];
	int scmaxd = Mu_S_ij_short[0][0];
	double scsum = 0;
	for (uint i = 0; i < 20; ++i)
		{
		for (uint j = 0; j < 20; ++j)
			{
			int sc = Mu_S_ij_short[i][j];
			scmin = min(scmin, sc);
			scmax = max(scmax, sc);
			scsum += sc;
			if (i == j) scmaxd = max(scmaxd, sc);
			}
		}
	ProgressLog("min %d, max %d, maxd %d, avg %.1f\n",
				scmin, scmax, scmaxd, scsum/400);
	}

static uint GetRandomKmer(uint k)
	{
	uint Kmer = 0;
	for (uint i = 0; i < k; ++i)
		Kmer = Kmer*36 + randu32()%36;
	return Kmer;
	}

static void KmerSelfScoreDist(uint k)
	{
	vector<float> Scores;
	uint DictSize64 = 1;
	for (uint i = 0; i < k; ++i)
		DictSize64 *= 36;
	asserta(DictSize64 < UINT32_MAX);
	uint DictSize = uint(DictSize64);
	Scores.reserve(DictSize);
	for (uint Kmer = 0; Kmer < DictSize; ++Kmer)
		{
		float Score = (float) get_kmer_pair_score(k, Kmer, Kmer);
		Scores.push_back(Score);
		}
	QuartsFloat QF;
	GetQuartsFloat(Scores, QF);
	ProgressLog("k=%u self ", k);
	QF.ProgressLogMe();
	}

static void KmerSelfScoreDistSample(uint k, uint N)
	{
	vector<float> Scores;
	Scores.reserve(N);
	uint DictSize64 = 1;
	for (uint i = 0; i < k; ++i)
		DictSize64 *= 36;
	asserta(DictSize64 < UINT32_MAX);
	uint DictSize = uint(DictSize64);
	for (uint i = 0; i < N; ++i)
		{
		uint Kmer = randu32()%DictSize;
		float Score = (float) get_kmer_pair_score(k, Kmer, Kmer);
		Scores.push_back(Score);
		}
	QuartsFloat QF;
	GetQuartsFloat(Scores, QF);
	ProgressLog("k=%u self(sample) ", k);
	QF.ProgressLogMe();
	}

static void KmerScoreDist(uint k, uint N)
	{
	vector<float> Scores;
	Scores.reserve(N);
	for (uint i = 0; i < N; ++i)
		{
		uint Kmer1 = GetRandomKmer(k);
		uint Kmer2 = GetRandomKmer(k);
		float Score = (float) get_kmer_pair_score(k, Kmer1, Kmer2);
		Scores.push_back(Score);
		}
	QuartsFloat QF;
	GetQuartsFloat(Scores, QF);
	ProgressLog("k=%u pair ", k);
	QF.ProgressLogMe();
	}

void cmd_mumxtest()
	{
	Info();
	const uint SAMPLE_SIZE = 1000*1000;
	for (uint k = 3; k <= 6; ++k)
		{
		if (k >= 6)
			KmerSelfScoreDistSample(k, SAMPLE_SIZE);
		else
			KmerSelfScoreDist(k);
		KmerScoreDist(k, SAMPLE_SIZE);
		}
	}
