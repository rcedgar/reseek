#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "featuretrainer.h"

void DSSParams::LoadFeatures(const vector<string> &FNs,
							 const vector<float> &Weights)
	{
	m_Features.clear();
	m_Weights.clear();

	const uint N = SIZE(FNs);
	asserta(SIZE(Weights) == N);
	FeatureTrainer FT;
	vector<float> Freqs;
	vector<vector<float> > FreqMx;
	vector<vector<float> > ScoreMx;
	for (uint i = 0; i < N ; ++i)
		{
		const string &FN = FNs[i];
		float w = Weights[i];
		FEATURE F = StrToFeature(FN.c_str());
		FT.FromTsv(FN);

		FT.GetFreqs(Freqs);
		FT.GetFreqMx(FreqMx);
		FT.GetLogOddsMx(ScoreMx);

		DSS::SetFeature(F, Freqs, FreqMx, ScoreMx);

		m_Features.push_back(F);
		m_Weights.push_back(w);
		}
	}
