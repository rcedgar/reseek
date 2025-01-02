#include "myutils.h"
#include "dssaligner.h"
#include "pdbchain.h"
#include "sort.h"

float GetSelfRevScore(DSSAligner &DA, DSS &D, const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers);

/***
[0d97074]++
reseek \
	-mukmerfilter c:/int/reseek_bench/data/scop40.cal \
	-output mukmerfilter.tsv \
	-log mukmerfilter.log

sort -gk3 mukmerfilter.tsv > mukmerfilter_sorted.tsv

$c/int/reseek_bench/py/analyze_hits.py \
  --input mukmerfilter_sorted.tsv \
  --type evalue \
  --fields 1,2,3 \
  --lookup $c/int/reseek_bench/data/dom_scopid.tsv

SEPQ0.1=0.2161 SEPQ1=0.3341 SEPQ10=0.4516 S1FP=0.3691 N1FP=167859 area=9.69
***/
void cmd_mukmerfilter()
	{
	opt_fast = true;
	optset_fast = true;
	DSSParams Params;
	Params.SetFromCmdLine(10000);
	asserta(optset_output);
	FILE *fOut = CreateStdioFile(opt_output);
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	const uint N = SIZE(Chains);
	vector<vector<byte> > MuLettersVec(N);
	vector<vector<uint> > MuKmersVec(N);
	vector<vector<vector<byte> > > Profiles(N);
	vector<float> SelfScores;

	DSSParams SelfParams = Params;
	SelfParams.m_UsePara = false;
	SelfParams.m_Omega = 0;
	SelfParams.m_OwnScoreMxs = false;

	DSSAligner DASelf;
	DASelf.SetParams(SelfParams);

	DSSAligner DA;
	DA.SetParams(Params);

	DSS D;
	D.SetParams(Params);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Mu k-mers");
		const PDBChain &Chain = *Chains[i];
		D.Init(Chain);
		D.GetMuLetters(MuLettersVec[i]);
		D.GetMuKmers(MuLettersVec[i], MuKmersVec[i]);
		D.GetProfile(Profiles[i]);

		float SelfScore = GetSelfRevScore(DASelf, D, Chain,
										  Profiles[i], &MuLettersVec[i], &MuKmersVec[i]);
		SelfScores.push_back(SelfScore);
		}

	MuKmerFilter MKF;
	MKF.SetParams(Params);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Filtering");
		const PDBChain &Chaini = *Chains[i];
		const vector<vector<byte> > &Profilei = Profiles[i];
		const vector<byte> &MuLettersi = MuLettersVec[i];
		const vector<uint> &MuKmersi = MuKmersVec[i];
		float SelfScorei = SelfScores[i];
		const string &Labeli = Chaini.m_Label;

		MKF.ResetQ();
		MKF.SetQ(&MuLettersi, &MuKmersi);

		DA.SetQuery(Chaini, &Profilei, &MuLettersi, &MuKmersi, SelfScorei);

		vector<int> Scorejs;
		vector<uint> js;
		for (uint j = 0; j < N; ++j)
			{
			if (j == i)
				continue;
			int Scorej = MKF.GetMaxHSPScore(MuLettersVec[j], MuKmersVec[j]);
			if (Scorej > 30)
				{
				js.push_back(j);
				Scorejs.push_back(Scorej);
				}
			}

		uint n = SIZE(js);
		asserta(SIZE(Scorejs) == n);
		if (n == 0)
			continue;
		uint *Order = myalloc(uint, n);
		QuickSortOrderDesc(Scorejs.data(), n, Order);

		uint m = n;
		if (m > 1000)
			m = 1000;

		int LastScore = INT_MAX;
		for (uint k = 0; k < m; ++k)
			{
			uint Idx = Order[k];
			uint j = js[Idx];
			asserta(j != i);
			const int Score = Scorejs[Idx];
			asserta(Score <= LastScore);
			LastScore = Score;

			const PDBChain &Chainj = *Chains[j];
			const string &Labelj = Chainj.m_Label;
			const vector<vector<byte> > &Profilej = Profiles[j];
			const vector<byte> &MuLettersj = MuLettersVec[j];
			const vector<uint> &MuKmersj = MuKmersVec[j];
			float SelfScorej = SelfScores[j];
			DA.SetTarget(Chainj, &Profilej, &MuLettersj, &MuKmersj, SelfScorej);
			DA.Align_NoAccel();
			float E = DA.m_EvalueA;
			if (E < 10)
				{
				fprintf(fOut, "%s", Labeli.c_str());
				fprintf(fOut, "\t%s", Labelj.c_str());
				fprintf(fOut, "\t%.3g", E);
				fprintf(fOut, "\n");
				}
			}

		myfree(Order);
		}
	CloseStdioFile(fOut);
	}
