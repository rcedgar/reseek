#include "myutils.h"
#include "alpha.h"
#include "quarts.h"
#include "mermx.h"
#include "mumx.h"
#include <chrono>

static uint s_dict_size;

static char *kmer2str(uint kmer)
	{
	static char s[7];
	for (uint i = 0; i < 6; ++i)
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

static int get_kmer_self_score(uint kmer)
	{
	int sum = 0;
	for (uint i = 0; i < 6; ++i)
		{
		uint code = kmer%36;
		sum += Mu_S_ij_short[code][code];
		kmer /= 36;
		}
	return sum;
	}

static int get_kmer_pair_score(uint kmer1, uint kmer2)
	{
	int sum = 0;
	for (uint i = 0; i < 6; ++i)
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
	asserta(s.size() == 6);
	uint kmer = str2kmer(s);
	int sc = get_kmer_self_score(kmer);
	ProgressLog("Self: %s = %3d (%3d)\n", s.c_str(), sc, should_be);
	}

static void test_pair_score(const string &s, int should_be,
							bool Trace = false)
	{
	asserta(s.size() == 6);
	uint kmer = str2kmer(s);
	int sc = get_kmer_self_score(kmer);
	uint nge78 = 0;
	for (uint kmer2 = 0; kmer2 < s_dict_size; ++kmer2)
		{
		int sc = get_kmer_pair_score(kmer, kmer2);
		if (sc >= 78)
			{
			if (Trace)
				Log(" %s = %d\n", kmer2str(kmer2), sc);
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

static void Test()
	{
	s_dict_size = myipow(36, 5);
	asserta(s_dict_size == 60466176);
	vector<float> kmer_self_scores;
	for (uint kmer = 0; kmer < s_dict_size; ++kmer)
		{
		ProgressStep(kmer, s_dict_size, "self scores");
		int self_score = get_kmer_self_score(kmer);
		kmer_self_scores.push_back(float(self_score));
		}
	QuartsFloat QF;
	GetQuartsFloat(kmer_self_scores, QF);
	QF.ProgressLogMe();

	test_self_score("VVVVV", 72);
	test_self_score("NVDDN", 112);
	test_self_score("WDTRD", 136);

	test_pair_score("CVPVV", 0);
	test_pair_score("VVVVC", 0);
	test_pair_score("SLVVV", 1);
	test_pair_score("VSVVC", 31);
	test_pair_score("SVVVQ", 72);
	test_pair_score("AVNPK", 4660);
	test_pair_score("VNPHD", 4299);
	test_pair_score("CQVND", 1118);
	}

static void TestNbr(MerMx &MM, const string &sKmer, uint FSn, uint FSTicks)
	{
	uint Kmer = MM.StrToKmer(sKmer);
	uint DictSize = myipow(36, 5);
	asserta(DictSize == 60466176);

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

void cmd_mumxtest()
	{
	Info();
	return;
	MerMx MM;
	MM.Init((short **) Mu_S_ij_short, 5, 36, 2);
	TestNbr(MM, "PALVV", 18, 9999);
	}
