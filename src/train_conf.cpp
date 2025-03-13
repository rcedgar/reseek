#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"
#include "kmeans.h"
#include "sort.h"
#include "trainer.h"

static vector<int> ivalues;
static vector<int> jvalues;
static uint s_r;

// All-vs-all distances in range -r .. +r
//  excluding adjacent pairs
static void Setijs(int r)
	{
	for (int i = -r; i <= r; ++i)
		{
		for (int j = i+1; j <= r; ++j)
			{
			int minij = min(i, j);
			int maxij = max(i, j);
			if (maxij - minij == 1)
				continue;
			ivalues.push_back(minij);
			jvalues.push_back(maxij);
			}
		}
	}

static void Getv(const PDBChain &Chain, uint Pos, uint r, vector<float> &v)
	{
	v.clear();
	const uint L = Chain.GetSeqLength();
	if (Pos < r || Pos + r >= L)
		return;
	const uint M = SIZE(ivalues);
	asserta(SIZE(jvalues) == M);
	int IntPos = Pos;
	for (uint m = 0; m < M; ++m)
		{
		int i = ivalues[m];
		int j = jvalues[m];
		float d = Chain.GetDist(IntPos+i, IntPos+j);
		v.push_back(d);
		}
	}

static void Getvs(const vector<PDBChain *> &Chains, uint N, uint r,
				  vector<vector<float> > &vs)
	{
	vs.clear();
	const uint ChainCount = SIZE(Chains);
	vs.reserve(N);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Label = Chain.m_Label;

		const uint L = Chain.GetSeqLength();
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			vector<float> v;
			Getv(Chain, Pos, r, v);
			if(v.empty())
				continue;
			if (SIZE(vs) == N)
				return;
			vs.push_back(v);
			}
		}
	}

static KMeans *s_ptrKM;
static const PDBChain *ChainQ;
static const PDBChain *ChainR;

static uint GetLetter(const vector<float> &v)
	{
	uint Letter = s_ptrKM->GetBestFitClusterIdx(v);
	return Letter;
	}

static void TrainerOnPair(
  const Trainer &T, uint ChainIdxQ, uint ChainIdxR,
  const vector<uint> &PosQs, const vector<uint> &PosRs)
	{
	ChainQ = &T.GetChain(ChainIdxQ);
	ChainR = &T.GetChain(ChainIdxR);
	}

static void TrainerAlphaCol(
  const Trainer &T, uint PosQ, uint PosR,
  uint &LetterQ, uint &LetterR)
	{
	vector<float> vQ;
	vector<float> vR;
	Getv(*ChainQ, PosQ, s_r, vQ);
	Getv(*ChainR, PosR, s_r, vR);
	if (vQ.empty() || vR.empty())
		{
		LetterQ = UINT_MAX;
		LetterR = UINT_MAX;
		return;
		}

	LetterQ = GetLetter(vQ);
	LetterR = GetLetter(vR);
	}

static void LogSrc(const KMeans &KM)
	{
	const uint K = KM.m_K;
	const uint M = KM.m_M;
	const uint N = KM.m_N;
	vector<uint> Order(K);
	const vector<uint> &Sizes = KM.m_Sizes;
	const vector<vector<float> > &Means = KM.m_Centroids;
	QuickSortOrderDesc(Sizes.data(), K, Order.data());

	Log("=========================================================\n");
	double TopPct = GetPct(Sizes[Order[0]], N);
	Log("Seed %u K=%u\n", opt(randseed), K);
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
	}

void cmd_train_conf()
	{
// Max number of positions to consider, for efficiency
	const uint N = optset_n ? opt(n) : 100000;
	asserta(optset_k);
	uint K = opt(k);
	s_r = (optset_r ? opt(r) : 3);

	if (!optset_randseed)
		{
		optset_randseed = true;
		opt_randseed = 1;
		}

	string FeatureName;
	if (optset_feature)
		FeatureName = string(opt(feature));
	else
		Ps(FeatureName, "Conf_K%u_r%u_seed%u", K, s_r, opt_randseed);

	Trainer Tr;
	Tr.Init(g_Arg1, opt(db));
	vector<PDBChain *> Chains = Tr.m_Chains;
	const uint ChainCount = SIZE(Chains);

	Setijs(s_r);

	vector<vector<float> > vs;
	Getvs(Chains, N, s_r, vs);

	KMeans KM;
	s_ptrKM = &KM;
	const uint MaxIters = 1000;
	KM.Run(vs, K, MaxIters);

	LogOdds LO;
	Tr.TrainLogOdds(K, TrainerOnPair, TrainerAlphaCol, LO);

	vector<vector<float> > ScoreMx;
	float ES = LO.GetLogOddsMx(ScoreMx);
	ProgressLog("K=%u, N=%u, seed=%u, r=%d ES=%.4g\n", K, N, opt_randseed, s_r, ES);
	LO.MxToSrc(g_fLog, "Conf", ScoreMx);
	LogSrc(KM);

	FILE *fOut = CreateStdioFile(opt(output));
	LO.WriteFeature(fOut, FeatureName);
	CloseStdioFile(fOut);

	FILE *fOut2 = CreateStdioFile(opt(output2));
	KM.WriteCentroids(fOut2);
	CloseStdioFile(fOut2);
	}
