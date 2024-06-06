#pragma once

static inline double fast_erf(double x)
	{
    const double a1 =  0.254829592;
    const double a2 = -0.284496736;
    const double a3 =  1.421413741;
    const double a4 = -1.453152027;
    const double a5 =  1.061405429;
    const double p  =  0.3275911;

    double t = 1.0 / (1.0 + p * abs(x));
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x*x);

    asserta(y > 0.0 && y < 1.0);
    return (x < 0.0 ? -y : y);
	}
