#pragma once

#include "pdbchain.h"

class TMA
	{
public:
	const PDBChain *m_Q = 0;
	const PDBChain *m_R = 0;
	string m_QRow;
	string m_RRow;
	double m_TM1 = 0;
	double m_TM2 = 0;

public:
	void NWDP_TM1(double **score, bool **path, double **val,
		int len1, int len2, double gap_open, int j2i[]);

	void NWDP_SE(bool **path, double **val, double **x, double **y,
		int len1, int len2, double d02, double gap_open, int j2i[]);

	void NWDP_TM2(bool **path, double **val, double **x, double **y,
		int len1, int len2, double t[3], double u[3][3],
		double d02, double gap_open, int j2i[]);

	void NWDP_TM3(bool **path, double **val, const char *secx, const char *secy,
		const int len1, const int len2, const double gap_open, int j2i[]);

	void clean_up_after_approx_TM(int* invmap0, int* invmap,
		double** score, bool** path, double** val, double** xtm, double** ytm,
		double** xt, double** r1, double** r2, const int xlen, const int minlen);

	double approx_TM(const int xlen, const int ylen, const int a_opt,
		double** xa, double** ya, double t[3], double u[3][3],
		const int invmap0[], const int mol_type);

	void copy_t_u(double t[3], double u[3][3], double t0[3], double u0[3][3]);

	double DP_iter(double** r1, double** r2, double** xtm, double** ytm,
		double** xt, bool** path, double** val, double** x, double** y,
		int xlen, int ylen, double t[3], double u[3][3], int invmap0[],
		int g1, int g2, int iteration_max, double local_d0_search,
		double D0_MIN, double Lnorm, double d0, double score_d8);

	double standard_TMscore(double** r1, double** r2, double** xtm, double** ytm,
		double** xt, double** x, double** y, int xlen, int ylen, int invmap[],
		int& L_ali, double& RMSD, double D0_MIN, double Lnorm, double d0,
		double d0_search, double score_d8, double t[3], double u[3][3],
		const int mol_type);

	void find_max_frag(double** x, int len, int* start_max,
		int* end_max, double dcu0);

	double get_initial_fgt(double** r1, double** r2, double** xtm, double** ytm,
		double** x, double** y, int xlen, int ylen,
		int* y2x, double d0, double d0_search,
		double dcu0, double t[3], double u[3][3]);

	//void get_initial_ssplus(double** r1, double** r2, double** score, bool** path,
	//	double** val, const char* secx, const char* secy, double** x, double** y,
	//	int xlen, int ylen, int* y2x0, int* y2x, const double D0_MIN, double d0);

	void score_matrix_rmsd_sec(double** r1, double** r2, double** score,
		const char* secx, const char* secy, double** x, double** y,
		int xlen, int ylen, int* y2x, const double D0_MIN, double d0);

	double TMscore8_search(double** r1, double** r2, double** xtm, double** ytm,
		double** xt, int Lali, double t0[3], double u0[3][3], int simplify_step,
		int score_sum_method, double* Rcomm, double local_d0_search, double Lnorm,
		double score_d8, double d0);

	double TMscore8_search_standard(double** r1, double** r2,
		double** xtm, double** ytm, double** xt, int Lali,
		double t0[3], double u0[3][3], int simplify_step, int score_sum_method,
		double* Rcomm, double local_d0_search, double score_d8, double d0);

	double detailed_search(double** r1, double** r2, double** xtm, double** ytm,
		double** xt, double** x, double** y, int xlen, int ylen,
		int invmap0[], double t[3], double u[3][3], int simplify_step,
		int score_sum_method, double local_d0_search, double Lnorm,
		double score_d8, double d0);

	double detailed_search_standard(double** r1, double** r2,
		double** xtm, double** ytm, double** xt, double** x, double** y,
		int xlen, int ylen, const int invmap0[], double t[3], double u[3][3],
		int simplify_step, int score_sum_method, double local_d0_search,
		const bool& bNormalize, double Lnorm, double score_d8, double d0);

	double get_score_fast(double** r1, double** r2, double** xtm, double** ytm,
		double** x, double** y, int xlen, int ylen, int invmap[],
		double d0, double d0_search, double t[3], double u[3][3]);

	double get_initial(double** r1, double** r2, double** xtm, double** ytm,
		double** x, double** y, int xlen, int ylen, int* y2x,
		double d0, double d0_search, double t[3], double u[3][3]);

	void smooth(int* sec, int len);

	char sec_str(double dis13, double dis14, double dis15,
		double dis24, double dis25, double dis35);

//	void make_sec(double** x, int len, char* sec);

	bool get_initial5(double** r1, double** r2, double** xtm, double** ytm,
		bool** path, double** val,
		double** x, double** y, int xlen, int ylen, int* y2x,
		double d0, double d0_search, const double D0_MIN);

	//void get_initial_ss(bool** path, double** val,
	//	const char* secx, const char* secy, int xlen, int ylen, int* y2x);

	bool Kabsch(double **x, double **y, int n, int mode, double *rms,
		double t[3], double u[3][3]);

	void parameter_set4search(const int xlen, const int ylen,
		double &D0_MIN, double &Lnorm,
		double &score_d8, double &d0, double &d0_search, double &dcu0);
	void parameter_set4final(const double len, double &D0_MIN, double &Lnorm,
		double &d0, double &d0_search, const int mol_type);
	void parameter_set4scale(const int len, const double d_s, double &Lnorm,
		double &d0, double &d0_search);
	int score_fun8( double **xa, double **ya, int n_ali, double d, int i_ali[],
		double *score1, int score_sum_method, const double Lnorm, 
		const double score_d8, const double d0);
	int score_fun8_standard(double **xa, double **ya, int n_ali, double d,
		int i_ali[], double *score1, int score_sum_method,
		double score_d8, double d0);
	int TMalign_main(double** xa, double** ya,
		const char* seqx, const char* secy,
		double& TM1, double& TM2,
		double& d0A, double& d0B,
		string& seqM, string& seqxA, string& seqyA,
		const int xlen, const int ylen);

	void WriteAln(FILE *f,
		const int xlen, const int ylen,
		const double TM1, const double TM2,
		double d0A, double d0B,
		const string &seqM, const string &seqxA, const string &seqyA) const;

	void LogTU(const char *Msg, const double t[3], const double u[3][3]) const;

	uint ReadCal(const string &FileName, char *Seq, double **a);

	double AlignChains(const PDBChain &Q, const PDBChain &R);
	};
