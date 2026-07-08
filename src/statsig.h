#pragma once

static const uint SCOP40c_DBSIZE = 8340;

class StatSig
	{
public:
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
