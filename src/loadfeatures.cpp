#include "myutils.h"
#include "dssparams.h"
#include "featuretrainer.h"

void FeatureTrainer::FromTsv(const string &FN)
	{
	FILE *f = OpenStdioFile(FN);

	string FeatureName;
	ReadStringValue(f, "feature", FeatureName);
	m_FeatureName = mystrsave(FeatureName.c_str());
	m_F = StrToFeature(FeatureName.c_str());

	string Type;
	ReadStringValue(f, "type", Type);
	if (Type == "int")
		m_IsInt = true;
	else if (Type == "float")
		m_IsInt = false;
	else
		Die("Bad feature type '%s'", Type.c_str());

	if (!m_IsInt)
		{
		m_MinValue = ReadFloatValue(f, "min");
		m_MedValue = ReadFloatValue(f, "med");
		m_MaxValue = ReadFloatValue(f, "max");
		m_UndefFreq = ReadFloatValue(f, "undef");
		}

	LogOdds::FromTsv(f);

	m_BinTs.clear();
	if (!m_IsInt)
		{
		for (uint i = 0; i + 2 < m_AlphaSize; ++i)
			m_BinTs.push_back(ReadFloatValue(f, "bint", i));
		}

	CloseStdioFile(f);
	}

void cmd_test()
	{
	FeatureTrainer FT;
	FT.FromTsv(g_Arg1);
	FT.ToTsv(opt(output));
	}
