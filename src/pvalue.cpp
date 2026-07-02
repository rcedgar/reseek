#include "myutils.h"
#include "flat_params.h"
#include "hitdata.h"

/***
src/reseek_tune2/fold_sf_fam_fit_loglin/
========================================
grep fit_lo *.txt

fam_fit.txt		fit_lo=200	fit_hi=500 p=10^(m*score + c) m=-0.0113 c=-0.00966
sf_fit1.txt		fit_lo=10	fit_hi=25 p=10^(m*score + c) m=-0.105 c=-0.796
sf_fit2.txt		fit_lo=25	fit_hi=50 p=10^(m*score + c) m=-0.0769 c=-1.5
fold_fit.txt	fit_lo=40	fit_hi=100 p=10^(m*score + c) m=-0.0448 c=-0.718
***/

double hitdata::calc_pvalue(double TS, PVALUE_MODE pvm)
	{
	double m = DBL_MAX;
	double c = DBL_MAX;
	switch (pvm)
		{
	case PVM_fam:
		m = -0.0113;
		c = -0.00966;
		break;

	case PVM_sf:
		if (TS < 25)
			{
			m = -0.105;
			c = -0.796;
			}
		else
			{
			m = -0.0769;
			c = -1.5;
			}
		break;

	case PVM_fold:
		m = -0.0448;
		c = -0.718;
		break;

	default:
		asserta(false);
		}

	double p = pow(10, m*TS + c);
	if (p > 1)
		p = 1;
	return p;
	}
