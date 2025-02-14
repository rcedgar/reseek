#include "myutils.h"
#include "alpha.h"
#include "quarts.h"
#include "mermx.h"
#include <chrono>

const MerMx &Get3DiMerMx();

int8_t threedi_substmx2[20][20] = {
{  6, -3, 1, 2, 3, -2, -2, -7, -3, -3, -10, -5, -1, 1, -4, -7, -5, -6, 0, -2}, // 0
{  -3, 6, -2, -8, -5, -4, -4, -12, -13, 1, -14, 0, 0, 1, -1, 0, -8, 1, -7, -9}, // 1
{  1, -2, 4, -3, 0, 1, 1, -3, -5, -4, -5, -2, 1, -1, -1, -4, -2, -3, -2, -2}, // 2
{  2, -8, -3, 9, -2, -7, -4, -12, -10, -7, -17, -8, -6, -3, -8, -10, -10, -13, -6, -3}, // 3
{  3, -5, 0, -2, 7, -3, -3, -5, 1, -3, -9, -5, -2, 2, -5, -8, -3, -7, 4, -4}, // 4
{  -2, -4, 1, -7, -3, 6, 3, 0, -7, -7, -1, -2, -2, -4, 3, -3, 4, -6, -4, -2}, // 5
{  -2, -4, 1, -4, -3, 3, 6, -4, -7, -6, -6, 0, -1, -3, 1, -3, -1, -5, -5, 3}, // 6
{  -7, -12, -3, -12, -5, 0, -4, 8, -5, -11, 7, -7, -6, -6, -3, -9, 6, -12, -5, -8}, // 7
{  -3, -13, -5, -10, 1, -7, -7, -5, 9, -11, -8, -12, -6, -5, -9, -14, -5, -15, 5, -8}, // 8
{  -3, 1, -4, -7, -3, -7, -6, -11, -11, 6, -16, -3, -2, 2, -4, -4, -9, 0, -8, -9}, // 9
{  -10, -14, -5, -17, -9, -1, -6, 7, -8, -16, 10, -9, -9, -10, -5, -10, 3, -16, -6, -9}, // 10
{  -5, 0, -2, -8, -5, -2, 0, -7, -12, -3, -9, 7, 0, -2, 2, 3, -4, 0, -8, -5}, // 11
{  -1, 0, 1, -6, -2, -2, -1, -6, -6, -2, -9, 0, 4, 0, 0, -2, -4, 0, -4, -5}, // 12
{  1, 1, -1, -3, 2, -4, -3, -6, -5, 2, -10, -2, 0, 5, -2, -4, -5, -1, -2, -5}, // 13
{  -4, -1, -1, -8, -5, 3, 1, -3, -9, -4, -5, 2, 0, -2, 6, 2, 0, -1, -6, -3}, // 14
{  -7, 0, -4, -10, -8, -3, -3, -9, -14, -4, -10, 3, -2, -4, 2, 6, -6, 0, -11, -9}, // 15
{  -5, -8, -2, -10, -3, 4, -1, 6, -5, -9, 3, -4, -4, -5, 0, -6, 8, -9, -5, -5}, // 16
{  -6, 1, -3, -13, -7, -6, -5, -12, -15, 0, -16, 0, 0, -1, -1, 0, -9, 3, -10, -11}, // 17
{  0, -7, -2, -6, 4, -4, -5, -5, 5, -8, -6, -8, -4, -2, -6, -11, -5, -10, 8, -6}, // 18
{  -2 -9, -2, -3, -4, -2, 3, -8, -8, -9, -9, -5, -5, -5, -3, -9, -5, -11, -6, 9}, // 19
};

int8_t threedi_substmx[20][20] = {
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
	return s;
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

/***
min -69, max 40, maxd 40, avg -15.2
Randomly selected 6-mer pairs:
	min -358, max 188, avg -91.3, nge78=0.000342, npk=21912
***/

/***
Self: VVVVVV =  72 ( 72)
Self: NVDDNV = 112 (112)
Self: WDTRDD = 136 (136)
***/

/***
Foldseek high-scoring k-mer pair threshold kmerThr=78

Nbrs: CVPVVV =     5 (    0)
Nbrs: VVVVCV =     1 (    0)
Nbrs: SLVVVV =    38 (    1)
Nbrs: VSVVCQ =   186 (   31)
Nbrs: SVVVQA =   206 (   72)
Nbrs: AVNPKD =  4660 ( 4660)
Nbrs: VNPHDT =  4323 ( 4299)
Nbrs: CQVNDH =  1891 ( 1118)
***/

/***
ticks for k-mer neighborhood in foldseek
ticks: std::chrono::high_resolution_clock::now() /auto elapsed = t1 - t0;

30 37 41 45 46 46 47 49 54 54 55 59 66 111 111 124 124 127 131 132 134 138 147 157 
160 171 176 193 199 201 216 218 224 225 226 237 239 244 245 246 250 253 255 261 266
277 280 284 285 287 293 297 306 307 314 337 338 356 367 385 396 417 418 418 440 463
489 505 549 616 618 631 647 689 690 711 726 782 810 863 866 878 911 922 942 976 978
1005 1014 1112 1294 1362 1418 1436 1569 1757 1770 1831 1892 2110 2140 2340 2391 2557
2690 2802 2815 2895 2984 3323 3705 3850 3864 3983 4070 4377 4852 5550 5645 5778 6342
8857 12780 14104 14499 21780 48680 82394
***/

static void Test()
	{
	const MerMx &MM = Get3DiMerMx();

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

	const int16_t *SelfScores = MM.BuildSelfScores_6mers();

	s_dict_size = myipow(20, 6);
	asserta(s_dict_size == 64000000);
	vector<float> kmer_self_scores;
	for (uint kmer = 0; kmer < s_dict_size; ++kmer)
		{
		ProgressStep(kmer, s_dict_size, "self scores");
		int self_score = get_kmer_self_score(kmer);
		asserta(SelfScores[kmer] == self_score);
		kmer_self_scores.push_back(float(self_score));
		}
	QuartsFloat QF;
	GetQuartsFloat(kmer_self_scores, QF);
	QF.ProgressLogMe();

	test_self_score("VVVVVV", 72);
	test_self_score("NVDDNV", 112);
	test_self_score("WDTRDD", 136);

	test_pair_score("CVPVVV", 0);
	test_pair_score("VVVVCV", 0);
	test_pair_score("SLVVVV", 1);
	test_pair_score("VSVVCQ", 31, true);
	test_pair_score("SVVVQA", 72);
	test_pair_score("AVNPKD", 4660);
	test_pair_score("VNPHDT", 4299);
	test_pair_score("CQVNDH", 1118);
	}

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

static MerMx *s_ptrMM3Di = 0;

const MerMx &Get3DiMerMx()
	{
	if (s_ptrMM3Di != 0)
		return *s_ptrMM3Di;
	s_ptrMM3Di = new MerMx;

	short **MxPtrs = myalloc(short *, 20);
	for (uint i = 0; i < 20; ++i)
		{
		short *Row = myalloc(short, 20);
		for (uint j = 0; j < 20; ++j)
			Row[j] = threedi_substmx[i][j];
		MxPtrs[i] = Row;
		}
	(*s_ptrMM3Di).Init(MxPtrs, 6, 20, 3);
	return *s_ptrMM3Di;
	}

void cmd_threedi()
	{
	Test();
	return;

	short **MxPtrs = myalloc(short *, 20);
	for (uint i = 0; i < 20; ++i)
		{
		short *Row = myalloc(short, 20);
		for (uint j = 0; j < 20; ++j)
			Row[j] = threedi_substmx[i][j];
		MxPtrs[i] = Row;
		}
	MerMx MM;
	MM.Init(MxPtrs, 6, 20, 3);
	//MM.LogMe();

	//uint Kmer = MM.StrToKmer("SLVVVV");
	//uint DictSize = myipow(20, 6);

	//uint *Kmers = myalloc(uint, DictSize);
	//uint *Kmers_Brute = myalloc(uint, DictSize);

	//uint n_Brute = MM.GetHighScoring6mers_Brute(Kmer, 78, Kmers_Brute, true);
	//ProgressLog("n_Brute=%u\n", n_Brute);

	//short *Work = myalloc(short, 2*MM.m_AS3);

	//auto c0 = std::chrono::high_resolution_clock::now();
	//uint n = MM.GetHighScoring6mers(Kmer, 78, Work, Kmers);
	//auto c1 = std::chrono::high_resolution_clock::now();
	//auto elapsed = c1 - c0;
	//uint32_t cticks = (uint32_t) elapsed.count();
	//ProgressLog("cticks=%u, n=%u\n", cticks, n);

	TestNbr(MM, "PALVVV", 18, 9999);
	//TestNbr(MM, "AVNPKD", 4660, 4020);
	//TestNbr(MM, "VNPHDT", 4299, 4017);
	//TestNbr(MM, "CQVNDH", 1118, 1429);

	//TestNbr(MM, "VVLVVV", 0, 110);
	//TestNbr(MM, "DLVVLV", 20, 321);
	//TestNbr(MM, "ADVLVL", 131, 484);
	//TestNbr(MM, "DDVSSN", 781, 1685);
	}
