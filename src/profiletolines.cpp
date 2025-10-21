#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "alpha.h"

char GetFeatureChar(byte Letter, uint AlphaSize);

void ProfileToLines(const vector<vector<byte> > &Profile,
	vector<string> &Lines)
	{
	Lines.clear();
	const uint FeatureCount = DSSParams::GetFeatureCount();
	asserta(SIZE(Profile) == FeatureCount);
	asserta(FeatureCount > 0);
	const uint L = SIZE(Profile[0]);
	Lines.resize(FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = DSSParams::m_Features[FeatureIdx];
		uint AS = DSSParams::GetAlphaSize(F);
		asserta(SIZE(Profile[FeatureIdx]) == L);
		string &Line = Lines[FeatureIdx];
		Line.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			byte Letter = Profile[FeatureIdx][Pos];
			Line.push_back(GetFeatureChar(Letter, AS));
			}
		asserta(SIZE(Line) == L);
		}
	}
