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
	static void SetMode(SEARCH_MODE Mode)
		{
		m_Mode = Mode;
		}

	static void SetDBSize(uint DBSize)
		{
		if (optset_dbsize)
			m_DBSize = opt_dbsize;
		else
			m_DBSize = DBSize;
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
	};
