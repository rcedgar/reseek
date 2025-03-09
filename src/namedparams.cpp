#include "myutils.h"
#include "dssparams.h"

void DSSParams::FromParamStr(const string &Str)
	{
	Clear();
	vector<string> Fields;
	Split(Str, Fields, '_');

	m_GapOpen = -1.5f;
	m_GapExt = -0.42f;
	//m_FwdMatchScore = 0;
	m_MinFwdScore = 0;
	m_Omega = 0;
	m_MKFPatternStr = "*";
	m_MuPrefilterPatternStr = "*";

	const uint N = SIZE(Fields);
	for (uint i = 0; i < N; ++i)
		{
		const string &Field = Fields[i];
		vector<string> Fields2;
		Split(Field, Fields2, ':');
		asserta(SIZE(Fields2) == 2);
		const string &Name = Fields2[0];
		double w = StrToFloat(Fields2[1]);
		FEATURE F = StrToFeature(Name.c_str());
		AddFeature(F, w);
		}
	}

void DSSParams::SetDefaults()
	{
	Clear();

	AddFeature(FEATURE_AA,			0.398145);
	AddFeature(FEATURE_NENDist,		0.129367);
	AddFeature(FEATURE_Conf,		0.202354);
	AddFeature(FEATURE_NENConf,		0.149383);
	AddFeature(FEATURE_RENDist,		0.0937677);
	AddFeature(FEATURE_DstNxtHlx,	0.00475462);
	AddFeature(FEATURE_StrandDens,	0.0183853);
	AddFeature(FEATURE_NormDens,	0.00384384);

	m_GapOpen = -0.685533f;
	m_GapExt = -0.051881f;
	m_MinFwdScore = 7.0f;
	m_MuPrefilterPatternStr = "1110011";
	m_MKFPatternStr = "111";
	}
