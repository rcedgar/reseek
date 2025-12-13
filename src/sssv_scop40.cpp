#include "myutils.h"
#include "fastbench.h"
#include "dssparams.h"
#include "triangle.h"

static void MakeSig(const vector<byte> &Seq, uint AS,
	vector<byte> &Sig)
	{
	Sig.clear();
	Sig.resize(AS, 0);
	const uint L = SIZE(Seq);
	for (uint i = 0; i < L; ++i)
		{
		byte Letter = Seq[i];
		asserta(Letter < AS);
		if (Sig[Letter] < 0xff)
			Sig[Letter] += 1;
		}
	}

float Align(const vector<vector<byte> > &Sigs, uint AS,
	const vector<vector<float> > &Mx, uint i, uint j)
	{
	asserta(i < SIZE(Sigs));
	asserta(j < SIZE(Sigs));
	const vector<byte> &Sig_i = Sigs[i];
	const vector<byte> &Sig_j = Sigs[j];
	asserta(SIZE(Sig_i) == AS);
	asserta(SIZE(Sig_j) == AS);
	float Score = 0;
	for (uint Letter = 0; Letter < AS; ++Letter)
		{
		uint n_i = Sig_i[Letter];
		uint n_j = Sig_j[Letter];
		if (n_i > 0 || n_j > 0)
			Score += Mx[Letter][Letter]*float(min(n_i, n_j))/max(n_i, n_j);
		}
	return Score;
	}

void cmd_sssv_scop40()
	{
	DSSParams::Init(DM_AlwaysVerysensitive);

	uint AS = DSSParams::m_AlphaSize_SSSA;
	asserta(AS > 0 && AS != UINT_MAX);

	const vector<vector<float> > &ScoreMx = DSSParams::m_ScoreMx_SSSA;
	asserta(SIZE(ScoreMx) == AS);

	const vector<vector<byte> > &ByteSeqs = DSSParams::m_IntSeqs_SSSA;
	const uint SeqCount = SIZE(ByteSeqs);

	const vector<string> &Labels = DSSParams::m_Labels_SSSA;
	asserta(SIZE(Labels) == SeqCount);

	vector<vector<byte> > Sigs(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		const string &Label = Labels[i];
		MakeSig(ByteSeqs[i], AS, Sigs[i]);
		}

	FastBench FB;
	FB.m_Labels = Labels;
	FB.SetLookupFromLabels();
	FB.Alloc();
	
	uint Counter = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		for (uint j = i; j < SeqCount; ++j)
			{
			ProgressStep(Counter++, FB.m_PairCount, "Aligning");
			float Score = Align(Sigs, AS, ScoreMx, i, j);
			FB.AppendHit(i, j, Score);
			}
		}

	FB.SetScoreOrder();

	uint nt = 0;
	uint nf = 0;
	for (uint k = 0; k < FB.m_PairCount; ++k)
		{
		uint HitIdx = FB.m_ScoreOrder[FB.m_PairCount-k-1];
		uint LabelIdx_i, LabelIdx_j;
		triangle_k_to_ij(HitIdx, FB.m_SeqCount, LabelIdx_i, LabelIdx_j);
		if (LabelIdx_i == LabelIdx_j)
			continue;
		uint SFIdx_i = FB.m_LabelIdxToSFIdx[LabelIdx_i];
		uint SFIdx_j = FB.m_LabelIdxToSFIdx[LabelIdx_j];
		if (SFIdx_i == SFIdx_j)
			nt += 2;
		else
			nf += 2;
		if (nt >= 100000)
			{
			ProgressLog("nt=%u nf=%u = %.1f%%\n", nt, nf, GetPct(nf, 69248284));
			break;
			}
		}

	FB.Bench();
	FB.WriteHits(opt(output));
	}
