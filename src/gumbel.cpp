#include "myutils.h"
#include "calibratesearcher.h"
#include <math.h>

double gumbel(double mu, double beta, double x)
	{
	double z = (x - mu)/beta;
	double e_z = exp(-z);
	double y = (1/beta)*exp(-(z + e_z));
	if(isnan(y))
		Die("gumbel(mu=%.3g, beta=%.3g, x=%.3g) = nan",
		  mu, beta, x);
	return y;
	}

double gumbel_cdf(double mu, double beta, double x)
	{
	double a = (x - mu)/beta;
	double y = exp(-exp(-a));
	return y;
	}

#if 0
// https://en.wikipedia.org/wiki/Gumbel_distribution
// reproduce Wikipedia figure f(x, mu=1.0, beta=2.0)
void cmd_test_gumbel()
	{
	FILE *f = CreateStdioFile(g_Arg1);
	double mu = 1;
	double beta = 2;
	double x = -5;
	while (x <= 20.0001)
		{
		double y = gumbel(mu, beta, x);
		fprintf(f, "%.2f\t%.4g\n", x, y);
		x += 0.1;
		}

	CloseStdioFile(f);
	}
#endif // 0

static double GetRMSE(double x0, double dx, const vector<double> &ys,
  double Mu, double Beta)
	{
	asserta(Beta > 0);
	const uint N = SIZE(ys);
	double Sume2 = 0;
	double x = x0;
	for (uint i = 0; i < N; ++i)
		{
		double y = ys[i];
		double yfit = gumbel(Mu, Beta, x);
		if (isnan(yfit))
			Die("gumbel(Mu=%.3g, Beta=%.3g, x=%.3g) = nan", Mu, Beta, x);
		double e = y*fabs((yfit - y));
		Sume2 += e*2;
		asserta(!isnan(Sume2));
		x += dx;
		}
	double RMSE = sqrt(Sume2/N);
	asserta(!isnan(RMSE));
	return RMSE;
	}

static double getmode(double x0, double dx,
  const vector<double> &ys)
	{
	const uint N = SIZE(ys);
	asserta(N > 0);
	double maxy = ys[0];
	double maxx = x0;
	double x = x0;
	for (uint i = 1; i < N; ++i)
		{
		if (ys[i] > maxy)
			{
			maxy = ys[i];
			maxx = x;
			}
		x += dx;
		}
	return maxx;
	}

static double getmean(double x0, double dx,
  const vector<double> &ys)
	{
	const uint N = SIZE(ys);
	asserta(N > 0);
	double sumy = 0;
	double sumxy = 0;
	double x = x0;
	for (uint i = 0; i < N; ++i)
		{
		double y = ys[i];
		sumxy += x*y;
		sumy += y;
		x += x0;
		}
	double mean = sumxy/sumy;
	return mean;
	}

/***
For perfect gumbel:
	mode = mu
	median = mu - beta*ln(ln(2))
	mean = mu + lambda*beta where lambda=0.5772156649

	beta = (mean - mu)/lambda
***/
void fit_gumbel(double x0, double dx,
  const vector<double> &ys, double &Mu, double &Beta)
	{
	const double lambda = 0.5772156649;
	const uint MAXITERS = 100;
	const double mean = getmean(x0, dx, ys);
	const double mode = getmode(x0, dx, ys);
	Mu = mode;
	double dMu = abs(Mu)/10;
	double d = mean - Mu;
	Beta = abs(d/lambda);
	if (Beta < 0.1)
		Beta = 0.1;
	double dBeta = Beta/4;
	uint StalledIters = 0;
	uint Iter = 0;
	for (;;)
		{
		if (++Iter > MAXITERS)
			break;
		double RMSE = GetRMSE(x0, dx, ys, Mu, Beta);
		double MuPlus = Mu + dMu;

		double MuMinus = Mu - dMu;
		if (MuMinus < 0.1)
			MuMinus = 0.1;

		double RMSE_MuPlus = GetRMSE(x0, dx, ys, MuPlus, Beta);
		double RMSE_MuMinus = GetRMSE(x0, dx, ys, MuMinus, Beta);

		if (RMSE <= RMSE_MuPlus && RMSE <= RMSE_MuMinus)
			{
			++StalledIters;
			dMu /= 2;
			}
		else if (RMSE_MuPlus <= RMSE_MuMinus)
			{
			StalledIters = 0;
			RMSE = RMSE_MuPlus;
			Mu = MuPlus;
			}
		else
			{
			asserta(RMSE_MuMinus <= RMSE_MuPlus);
			StalledIters = 0;
			RMSE = RMSE_MuMinus;
			Mu = MuMinus;
			}

		double BetaPlus = Beta + dBeta;
		double BetaMinus = Beta - dBeta;
		if (BetaMinus < 0.01)
			BetaMinus = 0.01;
		double RMSE_BetaPlus = GetRMSE(x0, dx, ys, Mu, BetaPlus);
		double RMSE_BetaMinus = GetRMSE(x0, dx, ys, Mu, BetaMinus);

		if (RMSE <= RMSE_BetaPlus && RMSE <= RMSE_BetaMinus)
			{
			++StalledIters;
			dBeta /= 2;
			}
		else if (RMSE_BetaPlus <= RMSE_BetaMinus)
			{
			StalledIters = 0;
			RMSE = RMSE_BetaPlus;
			Beta = BetaPlus;
			}
		else
			{
			asserta(RMSE_BetaMinus <= RMSE_BetaPlus);
			StalledIters = 0;
			RMSE = RMSE_BetaMinus;
			Beta = BetaMinus;
			}

		if (StalledIters > 2)
			break;
		}
	}

void cmd_test_gumbel()
	{
	double mu = 1.3;
	double beta = 0.8;
	double x0 = -5;
	double dx = 0.1;
	vector<double> xs;
	vector<double> ys;
	double x = x0;
	while (x < 20)
		{
		double y = gumbel(mu, beta, x);
		xs.push_back(x);
		ys.push_back(y);
		x += dx;
		}

	double FitMu, FitBeta;
	fit_gumbel(x0, dx, ys, FitMu, FitBeta);
	ProgressLog("FitMu %.3g, FitBeta %.3g\n", FitMu, FitBeta);
	}

void cmd_fit_gumbel()
	{
	vector<string> Lines;
	vector<string> Fields;
	ReadLinesFromFile(g_Arg1, Lines);
	vector<double> xs;
	vector<double> ys;
	const uint N = SIZE(Lines);

	const string &Line = Lines[0];
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 2);
	double x0 = StrToFloat(Fields[0]);
	double dx = StrToFloat(Fields[1]);

	for (uint i = 1; i < N; ++i)
		{
		const string &Line = Lines[i];
		double y = StrToFloat(Line);
		ys.push_back(y);
		}

	double FitMu, FitBeta;
	fit_gumbel(x0, dx, ys, FitMu, FitBeta);
	ProgressLog("FitMu %.3g, FitBeta %.3g\n", FitMu, FitBeta);
	if (!optset_output) 
		return;

	FILE *fOut = CreateStdioFile(opt_output);
	fprintf(fOut, "x\ty\tfity\n");
	double x = x0;
	for (uint i = 0; i < N; ++i)
		{
		double y = ys[i];
		double fity = gumbel(FitMu, FitBeta, x);
		fprintf(fOut, "%.3g\t%.3g\t%.3g\n", x, y, fity);
		x += dx;
		}
	CloseStdioFile(fOut);
	}
