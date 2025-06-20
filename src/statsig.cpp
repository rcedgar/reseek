#include "myutils.h"
#include "statsig.h"

#define	fast	0
#define	sensitive 1
#define	verysensitive 2

#define scop40	0
#define bfvd	1
#define pdb		2
#define afdb50	3
#define scop40c	4

uint StatSig::m_DBSize = UINT_MAX;
SEARCH_MODE StatSig::m_Mode = SM_undefined;
CALIBRATION_REF StatSig::m_Ref = REF_SCOP40;

static const char *searchdb2str(int searchdb)
	{
	switch (searchdb)
		{
	case scop40: return "SCOP40";
	case bfvd: return "BFVD";
	case pdb: return "PDB";
	case afdb50: return "AFDB50";
		}
	Die("searchdb2str(%d)", searchdb);
	return 0;
	}

static vector<int> dbsizes;
static vector<double> hitrates_fast;
static vector<double> hitrates_sensitive;
static vector<double> C_prefilter_ms_scop40;
static vector<double> C_prefilter_ms_scop40c;
static vector<double> C_prefilter_cs_scop40;
static vector<double> C_prefilter_cs_scop40c;
static double C_score_F_m_scop40;
static double C_score_F_m_scop40c;
static double C_score_F_c_scop40;
static double C_score_F_c_scop40c;

static void add_dbsize(int searchdb, int n)
	{
	assert(searchdb == dbsizes.size());
	dbsizes.push_back(n);
	}

static void add_hitrate(int mode, int refdb, double h)
	{
	if (mode == fast)
		hitrates_fast.push_back(h);
	else if (mode == sensitive)
		hitrates_sensitive.push_back(h);
	else
		Die("add_hitrate(mode=%d, refdb=%d, h=%.3g)",
		  mode, refdb, h);
	}

static void add_C_score_F_mc(int refdb, double m, double c)
	{
	if (refdb == scop40)
		{
		C_score_F_m_scop40 = m;
		C_score_F_c_scop40 = c;
		}
	else if (refdb == scop40c)
		{
		C_score_F_m_scop40c = m;
		C_score_F_c_scop40c = c;
		}
	else
		Die("add_C_score_F_mc(%d, %.3g, %.3g)", refdb, m, c);
	}

static void add_prefilter_mc(int searchdb, int refdb, double m, double c)
	{
	if (refdb == scop40)
		{
		assert(SIZE(C_prefilter_ms_scop40) == searchdb);
		assert(SIZE(C_prefilter_cs_scop40) == searchdb);
		C_prefilter_ms_scop40.push_back(m);
		C_prefilter_cs_scop40.push_back(c);
		}
	else if (refdb == scop40c)
		{
		assert(SIZE(C_prefilter_ms_scop40c) == searchdb);
		assert(SIZE(C_prefilter_cs_scop40c) == searchdb);
		C_prefilter_ms_scop40c.push_back(m);
		C_prefilter_cs_scop40c.push_back(c);
		}
	else
		Die("add_prefilter_mc(searchdb=%d, refdb=%d, m=%.3g, c=%.3g)",
		  searchdb, refdb, m, c);
	}

#define dbsize(searchdb, size) add_dbsize(searchdb, size);
#define hitrate(mode, refdb, h) add_hitrate(mode, refdb, h);
#define C_score_F_mc(refdb, m, c) add_C_score_F_mc(refdb, m, c);
#define prefilter_mc(searchdb, refdb, m, c) add_prefilter_mc(searchdb, refdb, m, c);

static bool Init()
	{
#include "reseek_calibrate.h"

#undef hitrate
#undef C_score_F_mc
#undef prefilter_mc
	return true;
	}
static bool InitDone = Init();

double interpolate(int query_x, const std::vector<int>& V, const std::vector<double>& D)
	{
    assert(V.size() == D.size() && V.size() > 1);
    // Check if V is sorted in ascending order
    for (size_t i = 0; i < V.size() - 1; ++i)
		if (V[i+1] <= V[i])
			Die("Vector V must be sorted in ascending order.");

    // --- Handle Out-of-Range Cases ---
    // Case 1: query_x is less than the first x value in V (extrapolation)
    if (query_x <= V[0])
		{
        // Use the gradient of the first interval (V[0], D[0]) and (V[1], D[1])
        double gradient = (D[1] - D[0]) / (V[1] - V[0]);
        return D[0] + gradient * (query_x - V[0]);
	    }

    // Case 2: query_x is greater than the last x value in V (extrapolation)
    if (query_x >= V.back()) // V.back() is equivalent to V[V.size() - 1]
		{
        // Use the gradient of the last interval (V[n-2], D[n-2]) and (V[n-1], D[n-1])
        double gradient = (D.back() - D[D.size() - 2]) / (V.back() - V[V.size() - 2]);
        return D.back() + gradient * (query_x - V.back());
		}

    // --- Handle In-Range Interpolation ---
    // Find the iterator to the first element in V that is strictly greater than query_x.
    // This gives us the upper bound of our interpolation interval.
    auto it_upper = std::upper_bound(V.begin(), V.end(), query_x);

    // The index of the upper bound point (x2, y2)
    // If query_x is exactly V[i], upper_bound will point to V[i+1], which is correct.
    size_t index2 = std::distance(V.begin(), it_upper);

    // The index of the lower bound point (x1, y1)
    // If query_x is exactly V[0], index2 will be 1, so index1 will be 0, which is correct.
    size_t index1 = index2 - 1;

    // Get the coordinates of the two points that bracket query_x
    int x1 = V[index1];
    double y1 = D[index1];
    int x2 = V[index2];
    double y2 = D[index2];

    // If query_x is exactly one of the data points, return its corresponding y value
    if (query_x == x1)
        return y1;
    if (query_x == x2)
        return y2;

    // Perform linear interpolation: y = y1 + (y2 - y1) * ((x - x1) / (x2 - x1))
    return y1 + (y2 - y1) * ((static_cast<double>(query_x) - x1) / (x2 - x1));
	}

static double get_hitrate(int mode, uint dbsize)
	{
	if (mode == fast)
		return interpolate(dbsize, dbsizes, hitrates_fast);
	else if (mode == sensitive)
		return interpolate(dbsize, dbsizes, hitrates_sensitive);
	else if (mode == verysensitive)
		return 1;
	else
		Die("get_hitrate(mode=%d, dbsize=%u)", mode, dbsize);
	return 0;
	}

static double get_Pvalue(double ts, int refdb)
	{
	double m = DBL_MAX;
	double c = DBL_MAX;
	if (refdb == scop40)
		{
		m = C_score_F_m_scop40;
		c = C_score_F_c_scop40;
		}
	else if (refdb == scop40c)
		{
		m = C_score_F_m_scop40c;
		c = C_score_F_c_scop40c;
		}
	else
		Die("get_C_score_F(ts=%.3g, refdb=%d)\n", ts, refdb);

	double log10C = -(m*ts + c);
	double P = min(pow(10, log10C), 1.0);
	return P;
	}

static double get_C_prefilter(double ts, int refdb, uint dbsize)
	{
	double m = DBL_MAX;
	double c = DBL_MAX;
	if (refdb == scop40)
		{
		m = interpolate(dbsize, dbsizes, C_prefilter_ms_scop40);
		c = interpolate(dbsize, dbsizes, C_prefilter_cs_scop40);
		}
	else if (refdb == scop40c)
		{
		m = interpolate(dbsize, dbsizes, C_prefilter_ms_scop40c);
		c = interpolate(dbsize, dbsizes, C_prefilter_cs_scop40c);
		}
	else
		Die("get_C_score_F(ts=%.3g, refdb=%d, dbsize=%u)", ts, refdb, dbsize);

	double log10C = m*ts + c;
	double C = min(pow(10, log10C), 1.0);
	return C;
	}

static double get_PF(int mode, uint dbsize)
	{
	if (mode == fast)
		return 0.5;
	else if (mode == sensitive)
		return 1;
	else if (mode == verysensitive)
		return 1;
	else
		Die("getPF(mode=%d, dbsize=%u)", mode, dbsize);
	return DBL_MAX;
	}

/***
C:\src\null_model\py\reseek_calibrate_v2.py

def estimate_evalue(ts, mode, searchdb, refdb):
	D = dbname2size.dbname2size[searchdb]
	PF = estimate_PF(mode, searchdb)
	h = get_hitrate(mode, searchdb)
	C_score_F_m, C_score_F_c = get_C_score_F_mc(refdb)
	C_score_F = get_loglin_C_score_F(ts, C_score_F_m, C_score_F_c)
	if mode == "fast":
		CDF_prefilter_m, CDF_prefilter_c = get_CDF_prefilter_mc(searchdb, refdb)
		CDF_prefilter = get_loglin_CDF_prefilter(ts, CDF_prefilter_m, CDF_prefilter_c)
		evalue = D*h*PF*CDF_prefilter
	elif mode == "sensitive":
		evalue = D*h*PF*C_score_F
	else:
		assert False, "mode=" + mode
	return evalue
***/

double get_Bayesian_Evalue(double ts, uint dbsize, int mode, int refdb)
	{
	double E = DBL_MAX;
	double P = get_Pvalue(ts, refdb);
	double h = get_hitrate(mode, dbsize);
	double PF = get_PF(mode, dbsize);
	if (mode == sensitive || mode == verysensitive)
		E = dbsize*h*PF*P;
	else if (mode == fast)
		{
		double C_prefilter = get_C_prefilter(ts, refdb, dbsize);
		E = dbsize*h*PF*C_prefilter;
		}
	else
		Die("get_Bayesian_Evalue(ts=%.3g, dbsize=%u, mode=%d)\n",
		  ts, dbsize, mode);
	return E;
	}

double StatSig::GetEvalue(double TS)
	{
	asserta(m_DBSize != DBL_MAX);
	asserta(m_Mode != SM_undefined);
	asserta(m_Ref != REF_undefined);

	int mode = -1;
	switch (StatSig::m_Mode)
		{
	case SM_fast: mode = 0; break;
	case SM_sensitive : mode = 1; break;
	case SM_verysensitive : mode = 2; break;
	default: asserta(false);
		}

	int refdb = -1;
	switch (StatSig::m_Ref)
		{
	case REF_SCOP40: refdb = 0; break;
	case REF_SCOP40c: refdb = 4; break;
	default: asserta(false);
		}

	double E = get_Bayesian_Evalue(TS, StatSig::m_DBSize, mode, refdb);
	return E;
	}

static void print_table(FILE *f, bool tsv, const char *title, bool show_scop40, bool show_scop40c,
						bool show_modes[3])
	{
	bool show_fast = show_modes[0];
	bool show_sensitive = show_modes[1];
	bool show_verysensitive = show_modes[2];
	asserta(show_scop40 || show_scop40c);
	asserta(show_fast || show_sensitive || show_verysensitive);

	fprintf(f, "\n\n");
	if (tsv)
		fprintf(f, "%s", "TS");
	else
		{
		fprintf(f, "%s\n", title);
		fprintf(f, "%4.4s", "TS");
		}

	if (tsv)
		{
		if (show_scop40)
			fprintf(f, "\t%s", "P-value");
		if (show_scop40c)
			fprintf(f, "\t%s", "P-value/c");
		}
	else
		{
		if (show_scop40)
			fprintf(f, "  %10.10s", "P-value");
		if (show_scop40c)
			fprintf(f, "  %10.10s", "P-value/c");
		fprintf(f, " |");
		}

	for (int searchdb = 0; searchdb < dbsizes.size(); ++searchdb)
		{
		const string name = searchdb2str(searchdb);
		const char *fmt = (tsv ? "\t%s" : "  %10.10s");
		if (show_scop40 && show_fast)
			fprintf(f, fmt, (name + "/F").c_str());
		if (show_scop40c && show_fast)
			fprintf(f, fmt, (name + "/Fc").c_str());
		if (show_scop40 && show_sensitive)
			fprintf(f, fmt, (name + "/S").c_str());
		if (show_scop40c && show_sensitive)
			fprintf(f, fmt, (name + "/Sc").c_str());
		if (show_scop40 && show_verysensitive)
			fprintf(f, fmt, (name + "/V").c_str());
		if (show_scop40c && show_verysensitive)
			fprintf(f, fmt, (name + "/Vc").c_str());
		if (!tsv)
			fprintf(f, " |");
		}
	fprintf(f, "\n");

	if (!tsv)
		{
		fprintf(f, "%4.4s", "----");
		if (show_scop40)
			fprintf(f, "  %10.10s", "----------");
		if (show_scop40c)
			fprintf(f, "  %10.10s", "----------");
		fprintf(f, " |");
		for (int searchdb = 0; searchdb < dbsizes.size(); ++searchdb)
			{
			const string name = searchdb2str(searchdb);
			if (show_scop40 && show_fast)
				fprintf(f, "  %10.10s", "----------");
			if (show_scop40c && show_fast)
				fprintf(f, "  %10.10s", "----------");
			if (show_scop40 && show_sensitive)
				fprintf(f, "  %10.10s", "----------");
			if (show_scop40c && show_sensitive)
				fprintf(f, "  %10.10s", "----------");
			if (show_scop40 && show_verysensitive)
				fprintf(f, "  %10.10s", "----------");
			if (show_scop40c && show_verysensitive)
				fprintf(f, "  %10.10s", "----------");
			fprintf(f, " |");
			}
		fprintf(f, "\n");
		}

	for (uint binidx = 1; binidx < 16; ++binidx)
		{
		double ts = binidx*0.05;
		if (tsv)
			fprintf(f, "%4.2f", ts);
		else
			fprintf(f, "%.2f", ts);
		double P_scop40 = get_Pvalue(ts, scop40);
		double P_scop40c = get_Pvalue(ts, scop40c);
		const char *efmt = (tsv ? "\t%.3e" : "  %10.3e");
		const char *sfmt = (tsv ? "\t%s" : "  %10.10s");
		if (show_scop40)
			fprintf(f, efmt, P_scop40);
		if (show_scop40c)
			fprintf(f, efmt, P_scop40c);
		if (!tsv)
			fprintf(f, " |");

		// style: 0=everything, 1=fast only, 2=sensitive only
		for (int searchdb = 0; searchdb < dbsizes.size(); ++searchdb)
			{
			int dbsize = dbsizes[searchdb];
			for (int mode = 0; mode < 3; ++mode)
				{
				if (!show_modes[mode])
					continue;
				if (show_scop40)
					{
					double E_scop40 = get_Bayesian_Evalue(ts, dbsize, mode, scop40);
					if (E_scop40 > 10)
						fprintf(f, sfmt, ">10");
					else
						fprintf(f, efmt, E_scop40);
					}
				if (show_scop40c)
					{
					double E_scop40c = get_Bayesian_Evalue(ts, dbsize, mode, scop40c);
					if (E_scop40c > 10)
						fprintf(f, sfmt, ">10");
					else
						fprintf(f, efmt, E_scop40c);
					}
				}
			if (!tsv)
				fprintf(f, " |");
			}
		fprintf(f, "\n");
		}
	}

static void print_hitrates(FILE *f)
	{
	fprintf(f, "Hit-rate (h)\n");
	fprintf(f, "%9.9s", "Mode");
	for (int searchdb = 0; searchdb < dbsizes.size(); ++searchdb)
		fprintf(f, "  %8.8s", searchdb2str(searchdb));
	fprintf(f, "\n");
	for (int mode = 0; mode < 2; ++mode)
		{
		fprintf(f, "%9.9s", mode == 0 ? "fast" : "sensitive");
		for (int searchdb = 0; searchdb < dbsizes.size(); ++searchdb)
			{
			int dbsize = dbsizes[searchdb];
				{
				double h = get_hitrate(mode, dbsize);
				fprintf(f, "  %8.3g", h);
				}
			}
		fprintf(f, "\n");
		}
	}

void cmd_bayes_report()
	{
	FILE *f = CreateStdioFile(g_Arg1);
	print_hitrates(f);

	bool all_modes[3] = { true, true, true };
	bool fast_only[3] = { true, false, false };
	bool sensitive_only[3] = { false, true, false };
	bool verysensitive_only[3] = { false, false, true };

	print_table(f, false, "All, F=fast, S=sensitive, V=verysensitive", true, true, all_modes);
	print_table(f, false, "fast, SCOP40 (/X) and SCOP40c (/Xc) reference, F=fast, S=sensitive, V=verysensitive", true, true, fast_only);
	print_table(f, false, "sensitive, SCOP40 and SCOP40c reference, F=fast, S=sensitive, V=verysensitive", true, true, sensitive_only);
	print_table(f, false, "verysensitive, SCOP40 (/X) and SCOP40c (/Xc) reference, F=fast, S=sensitive, V=verysensitive", true, true, verysensitive_only);
	print_table(f, false, "SCOP40 reference, F=fast, S=sensitive, V=verysensitive", true, false, all_modes);
	print_table(f, false, "SCOP40c reference, F=fast, S=sensitive, V=verysensitive", false, true, all_modes);

	print_table(f, true, "TSV", true, true, all_modes);
	CloseStdioFile(f);
	}
