#include "myutils.h"
#include "dss.h"

void cmd_dump_features()
	{
	const string &Prefix = g_Arg1;

	DSSParams::Init(DM_AlwaysSensitive);
	const uint FeatureCount = DSSParams::GetFeatureCount();
	for (uint FIdx = 0; FIdx < FeatureCount; ++FIdx)
		{
		FEATURE F = DSSParams::m_Features[FIdx];
		const string &FN = Prefix + FeatureToStr(F);
		Progress("%s\n", FN.c_str());
		DSSParams::DumpFeature(F, FN);
		}
	}
