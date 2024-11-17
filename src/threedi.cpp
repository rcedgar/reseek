#include "myutils.h"
#include "alpha.h"
#include "quarts.h"

static int8_t threedi_substmx[20][20] = {
{ 24, -12, 4,8, 12, -8, -8, -29, -12, -12, -41, -20, -4, 4, -16, -29, -20, -25, 0, -8 }, // i=0
{ -12, 24, -8, -33, -20, -16, -16, -49, -53, 4, -57, 0,0, 4, -4, 0, -33, 4, -29, -37 }, // i=1
{  4, -8, 16, -12, 0,4, 4, -12, -20, -16, -20, -8, 4, -4, -4, -16, -8, -12, -8, -8 }, // i=2
{  8, -33, -12, 36, -8, -29, -16, -49, -41, -29, -69, -33, -25, -12, -33, -41, -41, -53, -25, -12 }, // i=3
{ 12, -20, 0, -8, 28, -12, -12, -20, 4, -12, -37, -20, -8, 8, -20, -33, -12, -29, 16, -16 }, // i=4
{ -8, -16, 4, -29, -12, 24,12, 0, -29, -29, -4, -8, -8, -16, 12, -12, 16, -25, -16, -8 }, // i=5
{ -8, -16, 4, -16, -12, 12,24, -16, -29, -25, -25, 0, -4, -12, 4, -12, -4, -20, -20,12 }, // i=6
{ -29, -49, -12, -49, -20, 0, -16, 32, -20, -45, 28, -29, -25, -25, -12, -37, 24, -49, -20, -33 }, // i=7
{ -12, -53, -20, -41, 4, -29, -29, -20, 36, -45, -33, -49, -25, -20, -37, -57, -20, -61, 20, -33 }, // i=8
{ -12, 4, -16, -29, -12, -29, -25, -45, -45, 24, -65, -12, -8, 8, -16, -16, -37, 0, -33, -37 }, // i=9
{ -41, -57, -20, -69, -37, -4, -25, 28, -33, -65, 40, -37, -37, -41, -20, -41, 12, -65, -25, -37 }, // i=10
{ -20, 0, -8, -33, -20, -8, 0, -29, -49, -12, -37, 28, 0, -8, 8, 12, -16, 0, -33, -20 }, // i=11
{ -4, 0,4, -25, -8, -8, -4, -25, -25, -8, -37, 0, 16, 0,0, -8, -16, 0, -16, -20 }, // i=12
{  4, 4, -4, -12, 8, -16, -12, -25, -20, 8, -41, -8, 0, 20, -8, -16, -20, -4, -8, -20 }, // i=13
{ -16, -4, -4, -33, -20, 12, 4, -12, -37, -16, -20, 8,0, -8, 24, 8,0, -4, -25, -12 }, // i=14
{ -29, 0, -16, -41, -33, -12, -12, -37, -57, -16, -41, 12, -8, -16, 8, 24, -25, 0, -45, -37 }, // i=15
{ -20, -33, -8, -41, -12, 16, -4, 24, -20, -37, 12, -16, -16, -20, 0, -25, 32, -37, -20, -20 }, // i=16
{ -25, 4, -12, -53, -29, -25, -20, -49, -61, 0, -65, 0,0, -4, -4, 0, -37, 12, -41, -45 }, // i=17
{  0, -29, -8, -25, 16, -16, -20, -20, 20, -33, -25, -33, -16, -8, -25, -45, -20, -41, 32, -25 }, // i=18
{ -8, -37, -8, -12, -16, -8, 12, -33, -33, -37, -37, -20, -20, -20, -12, -37, -20, -45, -25,36 }, // i=19
};

static uint s_dict_size;

static char *kmer2str(uint kmer)
	{
	static char s[7];
	for (uint i = 0; i < 6; ++i)
		{
		s[i] = g_LetterToCharAmino[kmer%20];
		kmer /= 20;
		}
	}

static uint str2kmer(const string &s)
	{
	uint k = SIZE(s);
	uint kmer = 0;
	for (uint i = 0; i < k; ++i)
		{
		byte c = s[i];
		uint code = g_CharToLetterAmino[c];
		kmer = kmer*20 + code;
		}
	return kmer;
	}

static int get_kmer_self_score(uint kmer)
	{
	int sum = 0;
	for (uint i = 0; i < 6; ++i)
		{
		uint code = kmer%20;
		sum += threedi_substmx[code][code];
		kmer /= 20;
		}
	return sum;
	}

static int get_kmer_pair_score(uint kmer1, uint kmer2)
	{
	int sum = 0;
	for (uint i = 0; i < 6; ++i)
		{
		uint code1 = kmer1%20;
		uint code2 = kmer2%20;
		sum += threedi_substmx[code1][code2];
		kmer1 /= 20;
		kmer2 /= 20;
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

static void test_pair_score(const string &s, int should_be)
	{
	asserta(s.size() == 6);
	uint kmer = str2kmer(s);
	int sc = get_kmer_self_score(kmer);
	uint nge78 = 0;
	for (uint kmer2 = 0; kmer2 < s_dict_size; ++kmer2)
		{
		int sc = get_kmer_pair_score(kmer, kmer2);
		if (sc >= 78)
			++nge78;
		}
	ProgressLog("Nbrs: %s = %3d (%3d)\n", s.c_str(), nge78, should_be);
	}

// Foldseek high-scoring k-mer pair threshold kmerThr=78
// min -69, max 40, maxd 40, avg -15.2
// VVVVVV =  72 ( 72)
// NVDDNV = 112 (112)
// WDTRDD = 136 (136)
// Randomly selected 6-mer pairs:
//	min -358, max 188, avg -91.3, nge78=0.000342, npk=21912

void cmd_threedi()
	{
	int scmin = threedi_substmx[0][0];
	int scmax = threedi_substmx[0][0];
	int scmaxd = threedi_substmx[0][0];
	double scsum = 0;
	for (uint i = 0; i < 20; ++i)
		{
		for (uint j = 0; j < 20; ++j)
			{
			int sc = threedi_substmx[i][j];
			scmin = min(scmin, sc);
			scmax = max(scmax, sc);
			scsum += sc;
			if (i == j) scmaxd = max(scmaxd, sc);
			}
		}
	ProgressLog("min %d, max %d, maxd %d, avg %.1f\n",
				scmin, scmax, scmaxd, scsum/400);

	s_dict_size = myipow(20, 6);
	asserta(s_dict_size == 64000000);
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

	test_self_score("VVVVVV", 72);
	test_self_score("NVDDNV", 112);
	test_self_score("WDTRDD", 136);

	test_pair_score("VSVVCQ", 31);
	test_pair_score("SVVVQA", 72);
	test_pair_score("AVNPKD", 4660);
	test_pair_score("VNPHDT", 4299);
	test_pair_score("CQVNDH", 1118);
	}
