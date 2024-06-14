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

    if (y < 0.0 || y > 1.0)
		Die("fast_erf(%.8g) y=%.8g", x, y);
    return (x < 0.0 ? -y : y);
	}

// https://en.wikipedia.org/wiki/Q-function
// Q(x) is the probability that a standard normal variable
//   takes a value larger than x.
static inline double Q_func_standard(double x)
	{
	static const double SQRT2 = sqrt(2);
	//double Q = 0.5 - 0.5*erf(x/SQRT2);
	//double Q = 0.5 - 0.5*SUN_erf(x/SQRT2);
	double Q = 0.5 - 0.5*fast_erf(x/SQRT2);
	return Q;
	}

static inline double Q_func(double x, double mu, double sigma)
	{
	double x_standard = (mu - x)/sigma;
	return Q_func_standard(x_standard);
	}
