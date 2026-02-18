#include "myutils.h"
#include "statsig.h"

/***
C:\src\null_model3\py\fitted_p_value_params.py

                                        fitted log-linear elbow
                                        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
set_pvalue_params(name='ts', ref='sfc', x1=0.11, m0=-80, c0=-0.58, m=-52, c=-3.7)

# p-value calculation
def get_pvalue(name, ref, score):
	x1 = name_ref2x1[(name, ref)]
	if score < x1:
		m = name_ref2m0[(name, ref)]
		c = name_ref2c0[(name, ref)]
	else:
		m = name_ref2m[(name, ref)]
		c = name_ref2c[(name, ref)]
	p = math.exp(m*score + c)
	assert p >= 0
	if p > 1:
		p = 1
	return p
***/

double StatSig::GetPvalue(double TS)
	{
	static const double x1 = 0.11;
	static const double m0 =-80;
	static const double c0 =-0.58;
	static const double m =-52;
	static const double c =-3.7;

	double log10C = 0;
	if (TS < x1)
		log10C = m0*TS + c0;
	else
		log10C = m*TS + c;
	double P = pow(10, log10C);
	if (P > 1)
		P = 1;
	return P;
	}

double StatSig::GetEvalue(double TS)
	{
	double P = GetPvalue(TS);
	return P*SCOP40c_DBSIZE;
	}
