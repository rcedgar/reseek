#include "myutils.h"
#include "dss.h"
#include "alpha.h"

char GetFeatureChar(byte Letter, uint AlphaSize);

void cmd_mu_mapping()
	{
	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast);

	DSS D;
	D.SetParams(Params);
	//void GetMuLetters(uint MuLetter, vector<uint> &Letters) const;
	//uint GetMuLetter(const vector<uint> &Letters) const;
	uint AS = Params.m_MuAlphaSize;
	const uint N = Params.m_MuFeatureCount;
	Log("Mu");
	for (uint i = 0; i < N; ++i)
		Log("\t%s", FeatureToStr(Params.m_MuFeatures[i]));
	Log("\n");

	vector<uint> ASs;
	for (uint i = 0; i < N; ++i)
		ASs.push_back(DSS::GetAlphaSize(Params.m_MuFeatures[i]));

	for (uint Letter = 0; Letter < AS; ++Letter)
		{
		vector<uint> Letters;
		D.GetMuLetters(Letter, Letters);
		asserta(SIZE(Letters) == N);
		char c = GetFeatureChar(Letter, AS);
		Log("%c", c);
		for (uint i = 0; i < N; ++i)
			{
			uint AS = ASs[i];
			uint Letter = Letters[i];
			asserta(Letter < AS);
			char c = GetFeatureChar(Letter, AS);
			Log("\t%c", c);
			}
		Log("\n");
		}
	}
