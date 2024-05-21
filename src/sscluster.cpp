#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"
// #include "outputfiles.h"
#include "trainer.h"
#include "sort.h"

void GetScopDomFromLabel(const string &Label, string &Dom);

static vector<vector<double> > Means;
static DSS g_DSSQ;
static DSS g_DSSR;
static uint K;
static uint M;
static vector<int> ivalues;
static vector<int> jvalues;

static void Getv(const PDBChain &Chain, uint Pos,
  vector<double> &v)
	{
	v.clear();
	const uint L = Chain.GetSeqLength();
	if (Pos < 3 || Pos + 3 >= int(L))
		return;

	for (uint m = 0; m < M; ++m)
		{
		int i = ivalues[m];
		int j = jvalues[m];
		double d = Chain.GetDist(Pos+i, Pos+j);
		v.push_back(d);
		}
	asserta(SIZE(v) == M);
	}

static double GetDist(
  const vector<double> &v1,
  const vector<double> &v2)
	{
	const uint n = SIZE(v1);
	asserta(SIZE(v2) == n);
	double Sum2 = 0;
	for (uint i = 0; i < n; ++i)
		{
		double diff = v1[i] - v2[i];
		Sum2 += diff*diff;
		}
	return sqrt(Sum2);
	}

static uint GetLetter(const vector<double> &v)
	{
	asserta(SIZE(v) == M);
	double MinDist = DBL_MAX;
	uint BestCluster = UINT_MAX;
	for (uint k = 0; k < K; ++k)
		{
		double d = GetDist(v, Means[k]);
		if (k == 0 || d < MinDist)
			{
			BestCluster = k;
			MinDist = d;
			}
		}
	return BestCluster;
	}

static const PDBChain *ChainQ;
static const PDBChain *ChainR;

static void TrainerOnPair(
  const Trainer &T, uint ChainIdxQ, uint ChainIdxR)
	{
	ChainQ = &T.GetChain(ChainIdxQ);
	ChainR = &T.GetChain(ChainIdxR);
	}

static void TrainerAlphaCol(
  const Trainer &T, uint PosQ, uint PosR,
  uint &LetterQ, uint &LetterR)
	{
	vector<double> vQ;
	vector<double> vR;
	Getv(*ChainQ, PosQ, vQ);
	Getv(*ChainR, PosR, vR);
	if (vQ.empty() || vR.empty())
		{
		LetterQ = UINT_MAX;
		LetterR = UINT_MAX;
		return;
		}

	LetterQ = GetLetter(vQ);
	LetterR = GetLetter(vR);
	}

static void GetMeans(uint N, uint K, uint M,
  const vector<vector<double> > &vs,
  const vector<uint> &ClusterIdxs,
  vector<vector<double> > &Means)
	{
	Means.clear();
	Means.resize(K);
	vector<vector<double> > Sums(K);
	vector<uint> Counts(K);
	for (uint k = 0; k < K; ++k)
		{
		Sums[k].resize(M, 0);
		Means[k].resize(M, DBL_MAX);
		}

	for (uint i = 0; i < N; ++i)
		{
		uint ClusterIdx = ClusterIdxs[i];
		asserta(ClusterIdx < K);
		Counts[ClusterIdx] += 1;
		const vector<double> &v = vs[i];
		asserta(SIZE(v) == M);
		for (uint j = 0; j < M; ++j)
			{
			double x = v[j];
			Sums[ClusterIdx][j] += x;
			}
		}

	for (uint k = 0; k < K; ++k)
		{
		uint Count = Counts[k];
		asserta(Count > 0);
		for (uint j = 0; j < M; ++j)
			Means[k][j] = Sums[k][j]/Count;
		}
	}

static uint Assign(uint N, uint K, uint M,
  const vector<vector<double> > &vs,
  const vector<vector<double> > &Means,
  vector<uint> &ClusterIdxs,
  vector<uint> &Sizes)
	{
	Sizes.clear();
	Sizes.resize(K, 0);
	uint ChangeCount = 0;
	for (uint i = 0; i < N; ++i)
		{
		const vector<double> &v = vs[i];
		uint OldCluster = ClusterIdxs[i];
		double MinDist = DBL_MAX;
		uint BestCluster = UINT_MAX;
		for (uint k = 0; k < K; ++k)
			{
			double d = GetDist(v, Means[k]);
			if (k == 0 || d < MinDist)
				{
				BestCluster = k;
				MinDist = d;
				}
			}
		Sizes[BestCluster] += 1;
		if (BestCluster != OldCluster)
			{
			ClusterIdxs[i] = BestCluster;
			++ChangeCount;
			}
		}
	return ChangeCount;
	}

void cmd_sscluster()
	{
	Trainer Tr;
	Tr.Init(g_Arg1, opt_train_cal);

	vector<PDBChain *> Chains = Tr.m_Chains;
	const uint ChainCount = SIZE(Chains);
	vector<vector<double> > vs;
	const uint N = optset_n ? opt_n : 100000;
	asserta(optset_k);
	K = opt_k;
	vector<char> sss;
	const double MAXD = 16;

	uint m = 0;
	for (int i = -2; i <= 2; ++i)
		for (int j = i+1; j <= 2; ++j)
			{
			int minij = min(i, j);
			int maxij = max(i, j);
			if (maxij - minij == 1)
				Log("EXCLUDE Pos%d,%d\n", i, j);
			else
				{
				Log("Pos%d,%d\n", i, j);
				ivalues.push_back(minij);
				jvalues.push_back(maxij);
				}
			++m;
			}
	if (string(opt_myss3) == "Y")
		{
		ivalues.push_back(-3);
		jvalues.push_back(3);

		ivalues.push_back(0);
		jvalues.push_back(3);

		ivalues.push_back(-3);
		jvalues.push_back(0);
		}

	M = SIZE(ivalues);

	FILE *ftsv = CreateStdioFile(opt_output);
	map<string, uint> DomToChainIndex;
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Label = Chain.m_Label;

		string Dom;
		GetScopDomFromLabel(Label, Dom);
		DomToChainIndex[Dom] = ChainIndex;

		string SS;
		Chain.GetSS(SS);
		const uint L = Chain.GetSeqLength();
		for (int Pos = 0; Pos < int(L); ++Pos)
			{
			vector<double> v;
			Getv(Chain, Pos, v);
			if (v.empty())
				continue;
			vs.push_back(v);
			char ss = SS[Pos];
			sss.push_back(ss);
			if (ftsv != 0)
				{
				fprintf(ftsv, "%c", ss);
				for (uint i = 0; i < M; ++i)
					fprintf(ftsv, "\t%.4g", v[i]);
				fprintf(ftsv, "\n");
				}
			if (SIZE(vs) == N)
				break;
			}
		}
	CloseStdioFile(ftsv);
	ftsv = 0;

	vector<uint> ClusterIdxs;
	for (uint i = 0; i < N; ++i)
		{
		uint ClusterIdx = randu32()%K;
		ClusterIdxs.push_back(ClusterIdx);
		}

	const uint ITERS = 100;
	vector<uint> Sizes;
	for (uint Iter = 0; Iter < ITERS; ++Iter)
		{
		GetMeans(N, K, M, vs, ClusterIdxs, Means);
		Log("\n");
		Log("Iter %u\n", Iter);
		for (uint k = 0; k < K; ++k)
			{
			Log("Mean[%3u] ", k);
			for (uint j = 0; j < M; ++j)
				Log(" %10.4g", Means[k][j]);
			Log("\n");
			}
		uint ChangeCount = 
		  Assign(N, K, M, vs, Means, ClusterIdxs, Sizes);
		ProgressLog("Iter %u, %u changes\n", Iter, ChangeCount);
		Log("Sizes");
		for (uint k = 0; k < K; ++k)
			Log(" %u", Sizes[k]);
		Log("\n");
		if (ChangeCount == 0)
			{
			ProgressLog("=== CONVERGED ===\n");
			break;
			}
		}

	vector<uint> Order(K);
	QuickSortOrderDesc(Sizes.data(), K, Order.data());

	vector<vector<uint> > Correl(K);
	vector<vector<uint> > InvCorrel(4);
	for (uint k = 0; k < K; ++k)
		Correl[k].resize(4, 0);

	for (uint m = 0; m < 4; ++m)
		InvCorrel[m].resize(K, 0);

	for (uint i = 0; i < N; ++i)
		{
		char ss = sss[i];
		uint k = ClusterIdxs[i];
		asserta(k < K);
		switch (ss)
			{
		case 'h': Correl[k][0] += 1; InvCorrel[0][k] += 1; break;
		case 's': Correl[k][1] += 1; InvCorrel[1][k] += 1; break;
		case '~': Correl[k][2] += 1; InvCorrel[2][k] += 1; break;
		case 't': Correl[k][3] += 1; InvCorrel[3][k] += 1; break;
			}
		}

	Log("=========================================================\n");
	double TopPct = GetPct(Sizes[Order[0]], N);
	Log("Seed %u K=%u\n", opt_randseed, K);
	Log("Sizes");
	for (uint k = 0; k < K; ++k)
		Log(" %.1f", GetPct(Sizes[Order[k]], N));
	Log("\n");

	Log("\n");
	Log("//                     ");
	for (uint j = 0; j < M; ++j)
		{
		string tmp;
		Ps(tmp, "%d,%d", ivalues[j], jvalues[j]);
		Log("  %10.10s", tmp.c_str());
		}
	Log("\n");
	for (uint kk = 0; kk < K; ++kk)
		{
		uint k = Order[kk];
		Log("SSKMEAN(%3u", kk);
		Log(", %10u", Sizes[k]);
		for (uint j = 0; j < M; ++j)
			Log(", %10.4g", Means[k][j]);
		Log(");\n");
		}

	Log("\n");
	for (uint kk = 0; kk < K; ++kk)
		{
		uint k = Order[kk];
		Log("%2u: ", kk);
		for (uint m = 0; m < 4; ++m)
			Log(" %c(%7u)", "hs~t"[m], Correl[k][m]);
		Log("\n");
		}

	Log("\n");
	for (uint m = 0; m < 4; ++m)
		{
		Log("%c: ", "hs~t"[m]);
		uint Sum = 0;
		for (uint kk = 0; kk < K; ++kk)
			{
			uint k = Order[kk];
			uint n = InvCorrel[m][k];
			Log(" %2u(%7u)", k, n);
			Sum += n;
			}
		Log("  = %u\n", Sum);
		}

	LogOdds LO;
	Tr.TrainLogOdds(K, TrainerOnPair, TrainerAlphaCol, LO);
	vector<vector<double> > ScoreMx;
	double ExpectedScore = LO.GetLogOddsMx(ScoreMx);
	ProgressLog("K=%u myss3=%s M=%u N=%u seed=%u top=%.1f%% ES=%.3g\n",
	  K, opt_myss3, M, N, opt_randseed, TopPct, ExpectedScore);
	LO.MxToSrc(g_fLog, "MySS", ScoreMx);
	}
