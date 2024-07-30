#include "myutils.h"
#include "dss.h"
#include "alpha.h"

char GetFeatureChar(byte Letter, uint AlphaSize);

void cmd_mu_mapping()
	{
	DSSParams Params;
	Params.SetFromCmdLine(10000);

	DSS D;
	D.SetParams(Params);
	//void GetComboLetters(uint ComboLetter, vector<uint> &Letters) const;
	//uint GetComboLetter(const vector<uint> &Letters) const;
	uint AS = Params.m_ComboAlphaSize;
	const uint N = SIZE(Params.m_ComboFeatures);
	Log("Mu");
	for (uint i = 0; i < N; ++i)
		Log("\t%s", FeatureToStr(Params.m_ComboFeatures[i]));
	Log("\n");

	vector<uint> ASs;
	for (uint i = 0; i < N; ++i)
		ASs.push_back(DSS::GetAlphaSize(Params.m_ComboFeatures[i]));

	for (uint Letter = 0; Letter < AS; ++Letter)
		{
		vector<uint> Letters;
		D.GetComboLetters(Letter, Letters);
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
