#pragma once

static const uint SCOP40_DBSIZE = 11211;

enum SEARCH_MODE
	{
	SM_undefined = 0,
	SM_fast = 1,
	SM_sensitive = 2,
	SM_verysensitive = 3
	};

enum CALIBRATION_REF
	{
	REF_undefined = 0,
	REF_SCOP40 = 1,
	REF_SCOP40c = 2
	};

class StatSig
	{
public:
	static uint m_DBSize;
	static SEARCH_MODE m_Mode;
	static CALIBRATION_REF m_Ref;

public:
	static void Init(uint DBSize)
		{
		m_Mode = GetSearchModeFromCmdLine();
		if (optset_dbsize)
			m_DBSize = opt(dbsize);
		else
			m_DBSize = DBSize;
		}

	static void InitSensitive(uint DBSize)
		{
		m_Mode = SM_sensitive;
		if (optset_dbsize)
			m_DBSize = opt(dbsize);
		else
			m_DBSize = DBSize;
		}

	static SEARCH_MODE GetSearchModeFromCmdLine()
		{
		if (optset_fast)
			return SM_fast;
		else if (optset_sensitive)
			return SM_sensitive;
		else if (optset_verysensitive)
			return SM_verysensitive;
		Die("Unknown search mode");
		return SM_undefined;
		}

	static double GetQual(double TS)
		{
		const double a = 5.0f;
		const double b = -40.0f;
		double logE = a + b*TS;
		double DBSize = 11211;
		double Qual = 0;
		if (logE < -20)
			Qual = 1;
		else
			{
			double x = pow(10, logE/10);
			Qual = 1/(1 + x/2);
			}
		return Qual;
		}

	static double GetEvalue(double TS);
	static double GetPvalue(double TS);
	};
