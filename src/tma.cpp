#include "myutils.h"
#include "tma.h"

#pragma warning(disable: 4267)	// size_t -> int
#pragma warning(disable: 4244)	// double -> int

static void SetRCEtu(double t[3], double u[3][3])
	{
	t[0] = 0;
	t[1] = 0;
	t[2] = 0;

	u[0][0] = 1;
	u[0][1] = 0;
	u[0][2] = 0;

	u[1][0] = 0;
	u[1][1] = 1;
	u[1][2] = 0;

	u[2][0] = 0;
	u[2][0] = 0;
	u[2][2] = 1;
	}

void TMA::LogTU(const char *Msg, const double t[3], const double u[3][3]) const
	{
	Log("LogTU(%s):\n  t(%10.3f,  %10.3f,  %10.3f)\n", Msg, t[0], t[1], t[2]);
	Log(" u_x(%10.3f, %10.3f, %10.3f)\n", u[0][0], u[0][1], u[0][2]);
	Log(" u_y(%10.3f, %10.3f, %10.3f)\n", u[1][0], u[1][1], u[1][2]);
	Log(" u_z(%10.3f, %10.3f, %10.3f)\n", u[2][0], u[2][1], u[2][2]);
	}

static void split_white(const string& line, vector<string>& line_vec)
	{
	const char delimiter = ' ';
	bool within_word = false;
	for (int pos = 0; pos < line.size(); pos++)
		{
		if (line[pos] == delimiter)
			{
			within_word = false;
			continue;
			}
		if (!within_word)
			{
			within_word = true;
			line_vec.push_back("");
			}
		line_vec.back() += line[pos];
		}
	}

/* strip white space at the begining or end of string */
static string TrimWhiteSpace(const string& inputString)
	{
	string result = inputString;
	int idxBegin = inputString.find_first_not_of(" \n\r\t");
	int idxEnd = inputString.find_last_not_of(" \n\r\t");
	if (idxBegin >= 0 && idxEnd >= 0)
		result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
	return result;
	}

double TMA::dist(double x[3], double y[3])
	{
	double d1 = x[0] - y[0];
	double d2 = x[1] - y[1];
	double d3 = x[2] - y[2];

	return (d1 * d1 + d2 * d2 + d3 * d3);
	}

static double dot(double* a, double* b)
	{
	return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
	}

static void transform(double t[3], double u[3][3], double* x, double* x1)
	{
	x1[0] = t[0] + dot(&u[0][0], x);
	x1[1] = t[1] + dot(&u[1][0], x);
	x1[2] = t[2] + dot(&u[2][0], x);
	}

void TMA::do_rotation(double** x, double** x1, int len, double t[3], double u[3][3])
	{
	for (int i = 0; i < len; i++)
		{
		transform(t, u, &x[i][0], &x1[i][0]);
		}
	}

double TMA::TMscore8_search(
	double** r1, double** r2,
	double** xtm, double** ytm,
	double** xt, int Lali,
	double t0[3], double u0[3][3],
	int simplify_step, int score_sum_method,
	double* Rcomm, double local_d0_search,
	double Lnorm, double score_d8, double d0)
	{
	int i, m;
	double score_max, score, rmsd;
	const int kmax = Lali;
	//    int k_ali[kmax], ka, k;

	int ka, k;
	int* k_ali = myalloc(int, kmax);

	double t[3];
	double u[3][3];
	double d;

	//iterative parameters
	int n_it = 20;            //maximum number of iterations
	int n_init_max = 6; //maximum number of different fragment length 
	//    int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4   

	int* L_ini = myalloc(int, n_init_max);

	int L_ini_min = 4;
	if (Lali < L_ini_min) L_ini_min = Lali;

	int n_init = 0, i_init;
	for (i = 0; i < n_init_max - 1; i++)
		{
		n_init++;
		L_ini[i] = (int)(Lali / pow(2.0, (double)i));
		if (L_ini[i] <= L_ini_min)
			{
			L_ini[i] = L_ini_min;
			break;
			}
		}
	if (i == n_init_max - 1)
		{
		n_init++;
		L_ini[i] = L_ini_min;
		}

	score_max = -1;
//find the maximum score starting from local structures superposition
//    int i_ali[kmax], n_cut;
	int n_cut;
	int* i_ali = myalloc(int, kmax);
	int L_frag; //fragment length
	int iL_max; //maximum starting postion for the fragment

	for (i_init = 0; i_init < n_init; i_init++)
		{
		L_frag = L_ini[i_init];
		iL_max = Lali - L_frag;

		i = 0;
		while (1)
			{
			//extract the fragment starting from position i 
			ka = 0;
			for (k = 0; k < L_frag; k++)
				{
				int kk = k + i;
				r1[k][0] = xtm[kk][0];
				r1[k][1] = xtm[kk][1];
				r1[k][2] = xtm[kk][2];

				r2[k][0] = ytm[kk][0];
				r2[k][1] = ytm[kk][1];
				r2[k][2] = ytm[kk][2];

				k_ali[ka] = kk;
				ka++;
				}

			//extract rotation matrix based on the fragment
			Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);
			if (simplify_step != 1)
				*Rcomm = 0;
			do_rotation(xtm, xt, Lali, t, u);

			//get subsegment of this fragment
			d = local_d0_search - 1;
			n_cut = score_fun8(xt, ytm, Lali, d, i_ali, &score,
				score_sum_method, Lnorm, score_d8, d0);
			if (score > score_max)
				{
				score_max = score;

				//save the rotation matrix
				for (k = 0; k < 3; k++)
					{
					t0[k] = t[k];
					u0[k][0] = u[k][0];
					u0[k][1] = u[k][1];
					u0[k][2] = u[k][2];
					}
				}

			//try to extend the alignment iteratively            
			d = local_d0_search + 1;
			for (int it = 0; it < n_it; it++)
				{
				ka = 0;
				for (k = 0; k < n_cut; k++)
					{
					m = i_ali[k];
					r1[k][0] = xtm[m][0];
					r1[k][1] = xtm[m][1];
					r1[k][2] = xtm[m][2];

					r2[k][0] = ytm[m][0];
					r2[k][1] = ytm[m][1];
					r2[k][2] = ytm[m][2];

					k_ali[ka] = m;
					ka++;
					}
				//extract rotation matrix based on the fragment                
				Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);
				do_rotation(xtm, xt, Lali, t, u);
				n_cut = score_fun8(xt, ytm, Lali, d, i_ali, &score,
					score_sum_method, Lnorm, score_d8, d0);
				if (score > score_max)
					{
					score_max = score;
#if 0
					if (g_Trace)
						Log("score_max=%.3g\n", score_max);
#endif
					//save the rotation matrix
					for (k = 0; k < 3; k++)
						{
						t0[k] = t[k];
						u0[k][0] = u[k][0];
						u0[k][1] = u[k][1];
						u0[k][2] = u[k][2];
						}
					}

				//check if it converges            
				if (n_cut == ka)
					{
					for (k = 0; k < n_cut; k++)
						{
						if (i_ali[k] != k_ali[k]) break;
						}
					if (k == n_cut) break;
					}
				} //for iteration            

			if (i < iL_max)
				{
				i = i + simplify_step; //shift the fragment        
				if (i > iL_max) i = iL_max;  //do this to use the last missed fragment
				}
			else if (i >= iL_max) break;
			}//while(1)
			//end of one fragment
		}//for(i_init
	myfree(k_ali);
	myfree(L_ini);
	myfree(i_ali);
	return score_max;
	}

double TMA::TMscore8_search_standard(double** r1, double** r2,
	double** xtm, double** ytm, double** xt, int Lali,
	double t0[3], double u0[3][3], int simplify_step, int score_sum_method,
	double* Rcomm, double local_d0_search, double score_d8, double d0)
	{
	int i, m;
	double score_max, score, rmsd;
	const int kmax = Lali;
	//    int k_ali[kmax], ka, k;

	int ka, k;
	int* k_ali = myalloc(int, kmax);
	double t[3];
	double u[3][3];
	double d;

	//iterative parameters
	int n_it = 20;            //maximum number of iterations
	int n_init_max = 6; //maximum number of different fragment length 
	//    int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4   
	int* L_ini = myalloc(int, n_init_max);
	int L_ini_min = 4;
	if (Lali < L_ini_min) L_ini_min = Lali;

	int n_init = 0, i_init;
	for (i = 0; i < n_init_max - 1; i++)
		{
		n_init++;
		L_ini[i] = (int)(Lali / pow(2.0, (double)i));
		if (L_ini[i] <= L_ini_min)
			{
			L_ini[i] = L_ini_min;
			break;
			}
		}
	if (i == n_init_max - 1)
		{
		n_init++;
		L_ini[i] = L_ini_min;
		}

	score_max = -1;
	//find the maximum score starting from local structures superposition
//    int i_ali[kmax], n_cut;
	int n_cut;
	int* i_ali = myalloc(int, kmax);
	int L_frag; //fragment length
	int iL_max; //maximum starting postion for the fragment

	for (i_init = 0; i_init < n_init; i_init++)
		{
		L_frag = L_ini[i_init];
		iL_max = Lali - L_frag;

		i = 0;
		while (1)
			{
			//extract the fragment starting from position i 
			ka = 0;
			for (k = 0; k < L_frag; k++)
				{
				int kk = k + i;
				r1[k][0] = xtm[kk][0];
				r1[k][1] = xtm[kk][1];
				r1[k][2] = xtm[kk][2];

				r2[k][0] = ytm[kk][0];
				r2[k][1] = ytm[kk][1];
				r2[k][2] = ytm[kk][2];

				k_ali[ka] = kk;
				ka++;
				}
			//extract rotation matrix based on the fragment
			Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);
			if (simplify_step != 1)
				*Rcomm = 0;
			do_rotation(xtm, xt, Lali, t, u);

			//get subsegment of this fragment
			d = local_d0_search - 1;
			n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, &score,
				score_sum_method, score_d8, d0);

			if (score > score_max)
				{
				score_max = score;

				//save the rotation matrix
				for (k = 0; k < 3; k++)
					{
					t0[k] = t[k];
					u0[k][0] = u[k][0];
					u0[k][1] = u[k][1];
					u0[k][2] = u[k][2];
					}
				}

			//try to extend the alignment iteratively            
			d = local_d0_search + 1;
			for (int it = 0; it < n_it; it++)
				{
				ka = 0;
				for (k = 0; k < n_cut; k++)
					{
					m = i_ali[k];
					r1[k][0] = xtm[m][0];
					r1[k][1] = xtm[m][1];
					r1[k][2] = xtm[m][2];

					r2[k][0] = ytm[m][0];
					r2[k][1] = ytm[m][1];
					r2[k][2] = ytm[m][2];

					k_ali[ka] = m;
					ka++;
					}
				//extract rotation matrix based on the fragment                
				Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);
				do_rotation(xtm, xt, Lali, t, u);
				n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, &score,
					score_sum_method, score_d8, d0);
				if (score > score_max)
					{
					score_max = score;

					//save the rotation matrix
					for (k = 0; k < 3; k++)
						{
						t0[k] = t[k];
						u0[k][0] = u[k][0];
						u0[k][1] = u[k][1];
						u0[k][2] = u[k][2];
						}
					}

				//check if it converges            
				if (n_cut == ka)
					{
					for (k = 0; k < n_cut; k++)
						{
						if (i_ali[k] != k_ali[k]) break;
						}
					if (k == n_cut) break;
					}
				} //for iteration            

			if (i < iL_max)
				{
				i = i + simplify_step; //shift the fragment        
				if (i > iL_max) i = iL_max;  //do this to use the last missed fragment
				}
			else if (i >= iL_max) break;
			}//while(1)
			//end of one fragment
		}//for(i_init
	myfree(k_ali);
	myfree(L_ini);
	myfree(i_ali);
	return score_max;
	}

//Comprehensive TMscore search engine
// input:   two vector sets: x, y
//          an alignment invmap0[] between x and y
//          simplify_step: 1 or 40 or other integers
//          score_sum_method: 0 for score over all pairs
//                            8 for socre over the pairs with dist<score_d8
// output:  the best rotaion matrix t, u that results in highest TMscore
double TMA::detailed_search(double** r1, double** r2, double** xtm, double** ytm,
	double** xt, double** x, double** y, int xlen, int ylen,
	int invmap0[], double t[3], double u[3][3], int simplify_step,
	int score_sum_method, double local_d0_search, double Lnorm,
	double score_d8, double d0)
	{
	//x is model, y is template, try to superpose onto y
	int i, j, k;
	double tmscore;
	double rmsd;

	k = 0;
	for (i = 0; i < ylen; i++)
		{
		j = invmap0[i];
		if (j >= 0) //aligned
			{
			xtm[k][0] = x[j][0];
			xtm[k][1] = x[j][1];
			xtm[k][2] = x[j][2];

			ytm[k][0] = y[i][0];
			ytm[k][1] = y[i][1];
			ytm[k][2] = y[i][2];
			k++;
			}
		}

	//detailed search 40-->1
	tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u, simplify_step,
		score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0);
	return tmscore;
	}

double TMA::detailed_search_standard(
	double** r1, double** r2,
	double** xtm, double** ytm,
	double** xt, double** x, double** y,
	int /*xlen*/, int ylen,
	const int invmap0[],
	double t[3], double u[3][3],
	int simplify_step, int score_sum_method,
	double local_d0_search,
	const bool& bNormalize, double Lnorm, double score_d8, double d0)
	{
	//x is model, y is template, try to superpose onto y
	int i, j, k;
	double tmscore;
	double rmsd;

	k = 0;
	for (i = 0; i < ylen; i++)
		{
		j = invmap0[i];
		if (j >= 0) //aligned
			{
			xtm[k][0] = x[j][0];
			xtm[k][1] = x[j][1];
			xtm[k][2] = x[j][2];

			ytm[k][0] = y[i][0];
			ytm[k][1] = y[i][1];
			ytm[k][2] = y[i][2];
			k++;
			}
		}

	//detailed search 40-->1
	tmscore = TMscore8_search_standard(r1, r2, xtm, ytm, xt, k, t, u,
		simplify_step, score_sum_method, &rmsd, local_d0_search, score_d8, d0);
	if (bNormalize)// "-i", to use standard_TMscore, then bNormalize=true, else bNormalize=false; 
		tmscore = tmscore * k / Lnorm;

	return tmscore;
	}

//compute the score quickly in three iterations
double TMA::get_score_fast(double** r1, double** r2, double** xtm, double** ytm,
	double** x, double** y, int xlen, int ylen, int invmap[],
	double d0, double d0_search, double t[3], double u[3][3])
	{
	double rms, tmscore, tmscore1, tmscore2;
	int i, j, k;

	k = 0;
	for (j = 0; j < ylen; j++)
		{
		i = invmap[j];
		if (i >= 0)
			{
			r1[k][0] = x[i][0];
			r1[k][1] = x[i][1];
			r1[k][2] = x[i][2];

			r2[k][0] = y[j][0];
			r2[k][1] = y[j][1];
			r2[k][2] = y[j][2];

			xtm[k][0] = x[i][0];
			xtm[k][1] = x[i][1];
			xtm[k][2] = x[i][2];

			ytm[k][0] = y[j][0];
			ytm[k][1] = y[j][1];
			ytm[k][2] = y[j][2];

			k++;
			}
		else if (i != -1) Die("Wrong map!\n");
		}
	Kabsch(r1, r2, k, 1, &rms, t, u);

	//evaluate score   
	double di;
	const int len = k;
	//    double dis[len];    
	double* dis = myalloc(double, len);
	double d00 = d0_search;
	double d002 = d00 * d00;
	double d02 = d0 * d0;

	int n_ali = k;
	double xrot[3];
	tmscore = 0;
	for (k = 0; k < n_ali; k++)
		{
		transform(t, u, &xtm[k][0], xrot);
		di = dist(xrot, &ytm[k][0]);
		dis[k] = di;
		tmscore += 1 / (1 + di / d02);
		}

	//second iteration 
	double d002t = d002;
	while (1)
		{
		j = 0;
		for (k = 0; k < n_ali; k++)
			{
			if (dis[k] <= d002t)
				{
				r1[j][0] = xtm[k][0];
				r1[j][1] = xtm[k][1];
				r1[j][2] = xtm[k][2];

				r2[j][0] = ytm[k][0];
				r2[j][1] = ytm[k][1];
				r2[j][2] = ytm[k][2];

				j++;
				}
			}
		//there are not enough feasible pairs, relieve the threshold 
		if (j < 3 && n_ali>3) d002t += 0.5;
		else break;
		}

	if (n_ali != j)
		{
		Kabsch(r1, r2, j, 1, &rms, t, u);
		tmscore1 = 0;
		for (k = 0; k < n_ali; k++)
			{
			transform(t, u, &xtm[k][0], xrot);
			di = dist(xrot, &ytm[k][0]);
			dis[k] = di;
			tmscore1 += 1 / (1 + di / d02);
			}

		//third iteration
		d002t = d002 + 1;

		while (1)
			{
			j = 0;
			for (k = 0; k < n_ali; k++)
				{
				if (dis[k] <= d002t)
					{
					r1[j][0] = xtm[k][0];
					r1[j][1] = xtm[k][1];
					r1[j][2] = xtm[k][2];

					r2[j][0] = ytm[k][0];
					r2[j][1] = ytm[k][1];
					r2[j][2] = ytm[k][2];

					j++;
					}
				}
			//there are not enough feasible pairs, relieve the threshold 
			if (j < 3 && n_ali>3) d002t += 0.5;
			else break;
			}

		//evaluate the score
		Kabsch(r1, r2, j, 1, &rms, t, u);
		tmscore2 = 0;
		for (k = 0; k < n_ali; k++)
			{
			transform(t, u, &xtm[k][0], xrot);
			di = dist(xrot, &ytm[k][0]);
			tmscore2 += 1 / (1 + di / d02);
			}
		}
	else
		{
		tmscore1 = tmscore;
		tmscore2 = tmscore;
		}

	if (tmscore1 >= tmscore) tmscore = tmscore1;
	if (tmscore2 >= tmscore) tmscore = tmscore2;
	myfree(dis);
	return tmscore; // no need to normalize this score because it will not be used for latter scoring
	}

//perform gapless threading to find the best initial alignment
//input: x, y, xlen, ylen
//output: y2x0 stores the best alignment: e.g., 
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
double TMA::get_initial(double** r1, double** r2, double** xtm, double** ytm,
	double** x, double** y, int xlen, int ylen, int* y2x,
	double d0, double d0_search, double t[3], double u[3][3])
	{
	int min_len = min(xlen, ylen);
	if (min_len < 3) Die("Sequence is too short <3!\n");

	int min_ali = min_len / 2;              //minimum size of considered fragment 
	if (min_ali <= 5)  min_ali = 5;
	int n1, n2;
	n1 = -ylen + min_ali;
	n2 = xlen - min_ali;

	int i, j, k, k_best;
	double tmscore, tmscore_max = -1;

	k_best = n1;
	for (k = n1; k <= n2; k += 1)
		{
		//get the map
		for (j = 0; j < ylen; j++)
			{
			i = j + k;
			if (i >= 0 && i < xlen) y2x[j] = i;
			else y2x[j] = -1;
			}

		//evaluate the map quickly in three iterations
		//this is not real tmscore, it is used to evaluate the goodness of the initial alignment
		tmscore = get_score_fast(r1, r2, xtm, ytm,
			x, y, xlen, ylen, y2x, d0, d0_search, t, u);
		if (tmscore >= tmscore_max)
			{
			tmscore_max = tmscore;
			k_best = k;
			}
		}

	//extract the best map
	k = k_best;
	for (j = 0; j < ylen; j++)
		{
		i = j + k;
		if (i >= 0 && i < xlen) y2x[j] = i;
		else y2x[j] = -1;
		}

	return tmscore_max;
	}

void TMA::smooth(int* sec, int len)
	{
	int i, j;
	//smooth single  --x-- => -----
	for (i = 2; i < len - 2; i++)
		{
		if (sec[i] == 2 || sec[i] == 4)
			{
			j = sec[i];
			if (sec[i - 2] != j && sec[i - 1] != j && sec[i + 1] != j && sec[i + 2] != j)
				sec[i] = 1;
			}
		}

	//   smooth double 
	//   --xx-- => ------
	for (i = 0; i < len - 5; i++)
		{
		//helix
		if (sec[i] != 2 && sec[i + 1] != 2 && sec[i + 2] == 2 && sec[i + 3] == 2 &&
			sec[i + 4] != 2 && sec[i + 5] != 2)
			{
			sec[i + 2] = 1;
			sec[i + 3] = 1;
			}

		//beta
		if (sec[i] != 4 && sec[i + 1] != 4 && sec[i + 2] == 4 && sec[i + 3] == 4 &&
			sec[i + 4] != 4 && sec[i + 5] != 4)
			{
			sec[i + 2] = 1;
			sec[i + 3] = 1;
			}
		}

	//smooth connect
	for (i = 0; i < len - 2; i++)
		{
		if (sec[i] == 2 && sec[i + 1] != 2 && sec[i + 2] == 2) sec[i + 1] = 2;
		else if (sec[i] == 4 && sec[i + 1] != 4 && sec[i + 2] == 4) sec[i + 1] = 4;
		}
	}

char TMA::sec_str(double dis13, double dis14, double dis15,
	double dis24, double dis25, double dis35)
	{
	char s = 'C';

	double delta = 2.1;
	if (fabs(dis15 - 6.37) < delta && fabs(dis14 - 5.18) < delta &&
		fabs(dis25 - 5.18) < delta && fabs(dis13 - 5.45) < delta &&
		fabs(dis24 - 5.45) < delta && fabs(dis35 - 5.45) < delta)
		{
		s = 'H'; //helix                        
		return s;
		}

	delta = 1.42;
	if (fabs(dis15 - 13) < delta && fabs(dis14 - 10.4) < delta &&
		fabs(dis25 - 10.4) < delta && fabs(dis13 - 6.1) < delta &&
		fabs(dis24 - 6.1) < delta && fabs(dis35 - 6.1) < delta)
		{
		s = 'E'; //strand
		return s;
		}

	if (dis15 < 8) s = 'T'; //turn
	return s;
	}

///* secondary stucture assignment for protein:
// * 1->coil, 2->helix, 3->turn, 4->strand */
//void TMA::make_sec(double** x, int len, char* sec)
//	{
//	int j1, j2, j3, j4, j5;
//	double d13, d14, d15, d24, d25, d35;
//	for (int i = 0; i < len; i++)
//		{
//		sec[i] = 'C';
//		j1 = i - 2;
//		j2 = i - 1;
//		j3 = i;
//		j4 = i + 1;
//		j5 = i + 2;
//
//		if (j1 >= 0 && j5 < len)
//			{
//			d13 = sqrt(dist(x[j1], x[j3]));
//			d14 = sqrt(dist(x[j1], x[j4]));
//			d15 = sqrt(dist(x[j1], x[j5]));
//			d24 = sqrt(dist(x[j2], x[j4]));
//			d25 = sqrt(dist(x[j2], x[j5]));
//			d35 = sqrt(dist(x[j3], x[j5]));
//			sec[i] = sec_str(d13, d14, d15, d24, d25, d35);
//			}
//		}
//	sec[len] = 0;
//	}

//get initial alignment from secondary structure alignment
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g., 
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
//void TMA::get_initial_ss(bool** path, double** val,
//	const char* secx, const char* secy, int xlen, int ylen, int* y2x)
//	{
//	double gap_open = -1.0;
//	NWDP_TM3(path, val, secx, secy, xlen, ylen, gap_open, y2x);
//	}

// get_initial5 in TMalign fortran, get_initial_local in TMalign c by yangji
//get initial alignment of local structure superposition
//input: x, y, xlen, ylen
//output: y2x stores the best alignment: e.g., 
//y2x[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
bool TMA::get_initial5(double** r1, double** r2, double** xtm, double** ytm,
	bool** path, double** val,
	double** x, double** y, int xlen, int ylen, int* y2x,
	double d0, double d0_search, const double D0_MIN)
	{
	double GL, rmsd;
	double t[3];
	double u[3][3];

	double d01 = d0 + 1.5;
	if (d01 < D0_MIN) d01 = D0_MIN;
	double d02 = d01 * d01;

	double GLmax = 0;
	int aL = min(xlen, ylen);
	int* invmap = new int[ylen + 1];

	// jump on sequence1-------------->
	int n_jump1 = 0;
	if (xlen > 250)
		n_jump1 = 45;
	else if (xlen > 200)
		n_jump1 = 35;
	else if (xlen > 150)
		n_jump1 = 25;
	else
		n_jump1 = 15;
	if (n_jump1 > (xlen / 3))
		n_jump1 = xlen / 3;

	// jump on sequence2-------------->
	int n_jump2 = 0;
	if (ylen > 250)
		n_jump2 = 45;
	else if (ylen > 200)
		n_jump2 = 35;
	else if (ylen > 150)
		n_jump2 = 25;
	else
		n_jump2 = 15;
	if (n_jump2 > (ylen / 3))
		n_jump2 = ylen / 3;

	// fragment to superimpose-------------->
	int n_frag[2] = { 20, 100 };
	if (n_frag[0] > (aL / 3))
		n_frag[0] = aL / 3;
	if (n_frag[1] > (aL / 2))
		n_frag[1] = aL / 2;

	// start superimpose search-------------->
	bool flag = false;
	for (int i_frag = 0; i_frag < 2; i_frag++)
		{
		int m1 = xlen - n_frag[i_frag] + 1;
		int m2 = ylen - n_frag[i_frag] + 1;

		for (int i = 0; i < m1; i = i + n_jump1) //index starts from 0, different from FORTRAN
			{
			for (int j = 0; j < m2; j = j + n_jump2)
				{
				for (int k = 0; k < n_frag[i_frag]; k++) //fragment in y
					{
					r1[k][0] = x[k + i][0];
					r1[k][1] = x[k + i][1];
					r1[k][2] = x[k + i][2];

					r2[k][0] = y[k + j][0];
					r2[k][1] = y[k + j][1];
					r2[k][2] = y[k + j][2];
					}

				// superpose the two structures and rotate it
				Kabsch(r1, r2, n_frag[i_frag], 1, &rmsd, t, u);

				double gap_open = 0.0;
				NWDP_TM2(path, val, x, y, xlen, ylen,
					t, u, d02, gap_open, invmap);
				GL = get_score_fast(r1, r2, xtm, ytm, x, y, xlen, ylen,
					invmap, d0, d0_search, t, u);
				if (GL > GLmax)
					{
					GLmax = GL;
					for (int ii = 0; ii < ylen; ii++) y2x[ii] = invmap[ii];
					flag = true;
					}
				}
			}
		}

	delete[] invmap;
	return flag;
	}

void TMA::score_matrix_rmsd_sec(double** r1, double** r2, double** score,
	const char* secx, const char* secy, double** x, double** y,
	int xlen, int ylen, int* y2x, const double D0_MIN, double d0)
	{
	double t[3], u[3][3];
	double rmsd, dij;
	double d01 = d0 + 1.5;
	if (d01 < D0_MIN) d01 = D0_MIN;
	double d02 = d01 * d01;

	double xx[3];
	int i, k = 0;
	for (int j = 0; j < ylen; j++)
		{
		i = y2x[j];
		if (i >= 0)
			{
			r1[k][0] = x[i][0];
			r1[k][1] = x[i][1];
			r1[k][2] = x[i][2];

			r2[k][0] = y[j][0];
			r2[k][1] = y[j][1];
			r2[k][2] = y[j][2];

			k++;
			}
		}
	Kabsch(r1, r2, k, 1, &rmsd, t, u);


	for (int ii = 0; ii < xlen; ii++)
		{
		transform(t, u, &x[ii][0], xx);
		for (int jj = 0; jj < ylen; jj++)
			{
			dij = dist(xx, &y[jj][0]);
			if (secx[ii] == secy[jj])
				score[ii + 1][jj + 1] = 1.0 / (1 + dij / d02) + 0.5;
			else
				score[ii + 1][jj + 1] = 1.0 / (1 + dij / d02);
			}
		}
	}

////get initial alignment from secondary structure and previous alignments
////input: x, y, xlen, ylen
////output: y2x stores the best alignment: e.g., 
////y2x[j]=i means:
////the jth element in y is aligned to the ith element in x if i>=0 
////the jth element in y is aligned to a gap in x if i==-1
//void TMA::get_initial_ssplus(double** r1, double** r2, double** score, bool** path,
//	double** val, const char* secx, const char* secy, double** x, double** y,
//	int xlen, int ylen, int* y2x0, int* y2x, const double D0_MIN, double d0)
//	{
//	//create score matrix for DP
//	score_matrix_rmsd_sec(r1, r2, score, secx, secy, x, y, xlen, ylen,
//		y2x0, D0_MIN, d0);
//
//	double gap_open = -1.0;
//	NWDP_TM1(score, path, val, xlen, ylen, gap_open, y2x);
//	}

void TMA::find_max_frag(double** x, int len, int* start_max,
	int* end_max, double dcu0)
	{
	int r_min, fra_min = 4;           //minimum fragment for search
	int start;
	int Lfr_max = 0;

	r_min = (int)(len * 1.0 / 3.0); //minimum fragment, in case too small protein
	if (r_min > fra_min) r_min = fra_min;

	int inc = 0;
	double dcu0_cut = dcu0 * dcu0;;
	double dcu_cut = dcu0_cut;

	while (Lfr_max < r_min)
		{
		Lfr_max = 0;
		int j = 1;    //number of residues at nf-fragment
		start = 0;
		for (int i = 1; i < len; i++)
			{
			if (dist(x[i - 1], x[i]) < dcu_cut)
				{
				j++;

				if (i == (len - 1))
					{
					if (j > Lfr_max)
						{
						Lfr_max = j;
						*start_max = start;
						*end_max = i;
						}
					j = 1;
					}
				}
			else
				{
				if (j > Lfr_max)
					{
					Lfr_max = j;
					*start_max = start;
					*end_max = i - 1;
					}

				j = 1;
				start = i;
				}
			}// for i;

		if (Lfr_max < r_min)
			{
			inc++;
			double dinc = pow(1.1, (double)inc) * dcu0;
			dcu_cut = dinc * dinc;
			}
		}//while <;    
	}

//perform fragment gapless threading to find the best initial alignment
//input: x, y, xlen, ylen
//output: y2x0 stores the best alignment: e.g., 
//y2x0[j]=i means:
//the jth element in y is aligned to the ith element in x if i>=0 
//the jth element in y is aligned to a gap in x if i==-1
double TMA::get_initial_fgt(double** r1, double** r2, double** xtm, double** ytm,
	double** x, double** y, int xlen, int ylen,
	int* y2x, double d0, double d0_search,
	double dcu0, double t[3], double u[3][3])
	{
	int fra_min = 4;           //minimum fragment for search
	int fra_min1 = fra_min - 1;  //cutoff for shift, save time

	int xstart = 0, ystart = 0, xend = 0, yend = 0;

	find_max_frag(x, xlen, &xstart, &xend, dcu0);
	find_max_frag(y, ylen, &ystart, &yend, dcu0);


	int Lx = xend - xstart + 1;
	int Ly = yend - ystart + 1;
	int* ifr, * y2x_;
	int L_fr = min(Lx, Ly);
	ifr = new int[L_fr];
	y2x_ = new int[ylen + 1];

	//select what piece will be used. The original implement may cause 
	//asymetry, but only when xlen==ylen and Lx==Ly
	//if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1

	if (Lx < Ly || (Lx == Ly && xlen < ylen))
		{
		for (int i = 0; i < L_fr; i++) ifr[i] = xstart + i;
		}
	else if (Lx > Ly || (Lx == Ly && xlen > ylen))
		{
		for (int i = 0; i < L_fr; i++) ifr[i] = ystart + i;
		}
	else // solve asymetric for 1x5gA vs 2q7nA5
		{
		/* In this case, L0==xlen==ylen; L_fr==Lx==Ly */
		int L0 = xlen;
		double tmscore, tmscore_max = -1;
		int i, j, k;
		int n1, n2;
		int min_len;
		int min_ali;

		/* part 1, normalized by xlen */
		for (i = 0; i < L_fr; i++) ifr[i] = xstart + i;

		if (L_fr == L0)
			{
			n1 = (int)(L0 * 0.1); //my index starts from 0
			n2 = (int)(L0 * 0.89);
			j = 0;
			for (i = n1; i <= n2; i++)
				{
				ifr[j] = ifr[i];
				j++;
				}
			L_fr = j;
			}

		int L1 = L_fr;
		min_len = min(L1, ylen);
		min_ali = (int)(min_len / 2.5); //minimum size of considered fragment 
		if (min_ali <= fra_min1)  min_ali = fra_min1;
		n1 = -ylen + min_ali;
		n2 = L1 - min_ali;

		for (k = n1; k <= n2; k += 1)
			{
			//get the map
			for (j = 0; j < ylen; j++)
				{
				i = j + k;
				if (i >= 0 && i < L1) y2x_[j] = ifr[i];
				else             y2x_[j] = -1;
				}

			//evaluate the map quickly in three iterations
			tmscore = get_score_fast(r1, r2, xtm, ytm, x, y, xlen, ylen, y2x_,
				d0, d0_search, t, u);

			if (tmscore >= tmscore_max)
				{
				tmscore_max = tmscore;
				for (j = 0; j < ylen; j++) y2x[j] = y2x_[j];
				}
			}

		/* part 2, normalized by ylen */
		L_fr = Ly;
		for (i = 0; i < L_fr; i++) ifr[i] = ystart + i;

		if (L_fr == L0)
			{
			n1 = (int)(L0 * 0.1); //my index starts from 0
			n2 = (int)(L0 * 0.89);

			j = 0;
			for (i = n1; i <= n2; i++)
				{
				ifr[j] = ifr[i];
				j++;
				}
			L_fr = j;
			}

		int L2 = L_fr;
		min_len = min(xlen, L2);
		min_ali = (int)(min_len / 2.5); //minimum size of considered fragment 
		if (min_ali <= fra_min1)  min_ali = fra_min1;
		n1 = -L2 + min_ali;
		n2 = xlen - min_ali;

		for (k = n1; k <= n2; k++)
			{
			//get the map
			for (j = 0; j < ylen; j++) y2x_[j] = -1;

			for (j = 0; j < L2; j++)
				{
				i = j + k;
				if (i >= 0 && i < xlen) y2x_[ifr[j]] = i;
				}

			//evaluate the map quickly in three iterations
			tmscore = get_score_fast(r1, r2, xtm, ytm,
				x, y, xlen, ylen, y2x_, d0, d0_search, t, u);
			if (tmscore >= tmscore_max)
				{
				tmscore_max = tmscore;
				for (j = 0; j < ylen; j++) y2x[j] = y2x_[j];
				}
			}

		delete[] ifr;
		delete[] y2x_;
		return tmscore_max;
		}


	int L0 = min(xlen, ylen); //non-redundant to get_initial1
	if (L_fr == L0)
		{
		int n1 = (int)(L0 * 0.1); //my index starts from 0
		int n2 = (int)(L0 * 0.89);

		int j = 0;
		for (int i = n1; i <= n2; i++)
			{
			ifr[j] = ifr[i];
			j++;
			}
		L_fr = j;
		}


	//gapless threading for the extracted fragment
	double tmscore, tmscore_max = -1;

	if (Lx < Ly || (Lx == Ly && xlen <= ylen))
		{
		int L1 = L_fr;
		int min_len = min(L1, ylen);
		int min_ali = (int)(min_len / 2.5);              //minimum size of considered fragment 
		if (min_ali <= fra_min1)  min_ali = fra_min1;
		int n1, n2;
		n1 = -ylen + min_ali;
		n2 = L1 - min_ali;

		int i, j, k;
		for (k = n1; k <= n2; k += 1)
			{
			//get the map
			for (j = 0; j < ylen; j++)
				{
				i = j + k;
				if (i >= 0 && i < L1) y2x_[j] = ifr[i];
				else             y2x_[j] = -1;
				}

			//evaluate the map quickly in three iterations
			tmscore = get_score_fast(r1, r2, xtm, ytm, x, y, xlen, ylen, y2x_,
				d0, d0_search, t, u);

			if (tmscore >= tmscore_max)
				{
				tmscore_max = tmscore;
				for (j = 0; j < ylen; j++) y2x[j] = y2x_[j];
				}
			}
		}
	else
		{
		int L2 = L_fr;
		int min_len = min(xlen, L2);
		int min_ali = (int)(min_len / 2.5);              //minimum size of considered fragment 
		if (min_ali <= fra_min1)  min_ali = fra_min1;
		int n1, n2;
		n1 = -L2 + min_ali;
		n2 = xlen - min_ali;

		int i, j, k;

		for (k = n1; k <= n2; k++)
			{
			//get the map
			for (j = 0; j < ylen; j++) y2x_[j] = -1;

			for (j = 0; j < L2; j++)
				{
				i = j + k;
				if (i >= 0 && i < xlen) y2x_[ifr[j]] = i;
				}

			//evaluate the map quickly in three iterations
			tmscore = get_score_fast(r1, r2, xtm, ytm,
				x, y, xlen, ylen, y2x_, d0, d0_search, t, u);
			if (tmscore >= tmscore_max)
				{
				tmscore_max = tmscore;
				for (j = 0; j < ylen; j++) y2x[j] = y2x_[j];
				}
			}
		}


	delete[] ifr;
	delete[] y2x_;
	return tmscore_max;
	}

//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap
double TMA::DP_iter(double** r1, double** r2, double** xtm, double** ytm,
	double** xt, bool** path, double** val, double** x, double** y,
	int xlen, int ylen, double t[3], double u[3][3], int invmap0[],
	int g1, int g2, int iteration_max, double local_d0_search,
	double D0_MIN, double Lnorm, double d0, double score_d8)
	{
	double gap_open[2] = { -0.6, 0 };
	double rmsd;
	int* invmap = new int[ylen + 1];

	int iteration, i, j, k;
	double tmscore, tmscore_max, tmscore_old = 0;
	int score_sum_method = 8, simplify_step = 40;
	tmscore_max = -1;

	//double d01=d0+1.5;
	double d02 = d0 * d0;
	for (int g = g1; g < g2; g++)
		{
		for (iteration = 0; iteration < iteration_max; iteration++)
			{
			NWDP_TM2(path, val, x, y, xlen, ylen,
				t, u, d02, gap_open[g], invmap);

			k = 0;
			for (j = 0; j < ylen; j++)
				{
				i = invmap[j];

				if (i >= 0) //aligned
					{
					xtm[k][0] = x[i][0];
					xtm[k][1] = x[i][1];
					xtm[k][2] = x[i][2];

					ytm[k][0] = y[j][0];
					ytm[k][1] = y[j][1];
					ytm[k][2] = y[j][2];
					k++;
					}
				}

			tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
				simplify_step, score_sum_method, &rmsd, local_d0_search,
				Lnorm, score_d8, d0);


			if (tmscore > tmscore_max)
				{
				tmscore_max = tmscore;
				for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
				}

			if (iteration > 0)
				{
				if (fabs(tmscore_old - tmscore) < 0.000001) break;
				}
			tmscore_old = tmscore;
			}// for iteration           

		}//for gapopen


	delete[]invmap;
	return tmscore_max;
	}

double TMA::standard_TMscore(double** r1, double** r2, double** xtm, double** ytm,
	double** xt, double** x, double** y, int xlen, int ylen, int invmap[],
	int& L_ali, double& RMSD, double D0_MIN, double Lnorm, double d0,
	double d0_search, double score_d8, double t[3], double u[3][3],
	const int mol_type)
	{
	D0_MIN = 0.5;
	Lnorm = ylen;
	if (mol_type > 0) // RNA
		{
		if (Lnorm <= 11) d0 = 0.3;
		else if (Lnorm > 11 && Lnorm <= 15) d0 = 0.4;
		else if (Lnorm > 15 && Lnorm <= 19) d0 = 0.5;
		else if (Lnorm > 19 && Lnorm <= 23) d0 = 0.6;
		else if (Lnorm > 23 && Lnorm < 30)  d0 = 0.7;
		else d0 = (0.6 * pow((Lnorm * 1.0 - 0.5), 1.0 / 2) - 2.5);
		}
	else
		{
		if (Lnorm > 21) d0 = (1.24 * pow((Lnorm * 1.0 - 15), 1.0 / 3) - 1.8);
		else d0 = D0_MIN;
		if (d0 < D0_MIN) d0 = D0_MIN;
		}
	double d0_input = d0;// Scaled by seq_min

	double tmscore;// collected alined residues from invmap
	int n_al = 0;
	int i;
	for (int j = 0; j < ylen; j++)
		{
		i = invmap[j];
		if (i >= 0)
			{
			xtm[n_al][0] = x[i][0];
			xtm[n_al][1] = x[i][1];
			xtm[n_al][2] = x[i][2];

			ytm[n_al][0] = y[j][0];
			ytm[n_al][1] = y[j][1];
			ytm[n_al][2] = y[j][2];

			r1[n_al][0] = x[i][0];
			r1[n_al][1] = x[i][1];
			r1[n_al][2] = x[i][2];

			r2[n_al][0] = y[j][0];
			r2[n_al][1] = y[j][1];
			r2[n_al][2] = y[j][2];

			n_al++;
			}
		else if (i != -1) Die("Wrong map!\n");
		}
	L_ali = n_al;

	Kabsch(r1, r2, n_al, 0, &RMSD, t, u);
	RMSD = sqrt(RMSD / (1.0 * n_al));

	int temp_simplify_step = 1;
	int temp_score_sum_method = 0;
	d0_search = d0_input;
	double rms = 0.0;
	tmscore = TMscore8_search_standard(r1, r2, xtm, ytm, xt, n_al, t, u,
		temp_simplify_step, temp_score_sum_method, &rms, d0_input,
		score_d8, d0);
	tmscore = tmscore * n_al / (1.0 * Lnorm);

	return tmscore;
	}

/* copy the value of t and u into t0,u0 */
void TMA::copy_t_u(double t[3], double u[3][3], double t0[3], double u0[3][3])
	{
	int i, j;
	for (i = 0; i < 3; i++)
		{
		t0[i] = t[i];
		for (j = 0; j < 3; j++) u0[i][j] = u[i][j];
		}
	}

/* calculate approximate TM-score given rotation matrix */
double TMA::approx_TM(const int xlen, const int ylen, const int a_opt,
	double** xa, double** ya, double t[3], double u[3][3],
	const int invmap0[], const int mol_type)
	{
	double Lnorm_0 = ylen; // normalized by the second protein
	if (a_opt == -2 && xlen > ylen) Lnorm_0 = xlen;      // longer
	else if (a_opt == -1 && xlen < ylen) Lnorm_0 = xlen; // shorter
	else if (a_opt == 1) Lnorm_0 = (xlen + ylen) / 2.;     // average

	double D0_MIN;
	double Lnorm;
	double d0;
	double d0_search;
	parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
	double TMtmp = 0;
	double d;
	double xtmp[3] = { 0,0,0 };

	for (int i = 0, j = 0; j < ylen; j++)
		{
		i = invmap0[j];
		if (i >= 0)//aligned
			{
			transform(t, u, &xa[i][0], &xtmp[0]);
			d = sqrt(dist(&xtmp[0], &ya[j][0]));
			TMtmp += 1 / (1 + (d / d0) * (d / d0));
			//if (d <= score_d8) TMtmp+=1/(1+(d/d0)*(d/d0));
			}
		}
	TMtmp /= Lnorm_0;
	return TMtmp;
	}

void TMA::clean_up_after_approx_TM(int* invmap0, int* invmap,
	double** score, bool** path, double** val, double** xtm, double** ytm,
	double** xt, double** r1, double** r2, const int xlen, const int minlen)
	{
	delete[] invmap0;
	delete[] invmap;
	DeleteArray(&score, xlen + 1);
	DeleteArray(&path, xlen + 1);
	DeleteArray(&val, xlen + 1);
	DeleteArray(&xtm, minlen);
	DeleteArray(&ytm, minlen);
	DeleteArray(&xt, xlen);
	DeleteArray(&r1, minlen);
	DeleteArray(&r2, minlen);
	return;
	}

/* Entry function for TM-align. Return TM-score calculation status:
 * 0   - full TM-score calculation
 * 1   - terminated due to exception
 * 2-7 - pre-terminated due to low TM-score */
int TMA::TMalign_main(double** xa, double** ya,
	const char* seqx, const char* seqy,
	double& TM1, double& TM2,
	double& d0A, double& d0B,
	string& seqM, string& seqxA, string& seqyA,
	const int xlen, const int ylen)
	{
	double t0[3], u0[3][3];
	double d0_0, TM_0;
	double d0_out = 5.0;
	double rmsd0 = 0.0;
	double Liden = 0;
	int n_ali = 0;
	int n_ali8 = 0;

	TM2 = 0;

	////////////////////////////////////////////////////////
	double D0_MIN;        //for d0
	double Lnorm;         //normalization length
	double score_d8, d0, d0_search, dcu0;//for TMscore search
	double t[3], u[3][3]; //Kabsch translation vector and rotation matrix
	double** score;       // Input score table for dynamic programming
	bool** path;        // for dynamic programming  
	double** val;         // for dynamic programming  
	double** xtm, ** ytm;  // for TMscore search engine
	double** xt;          //for saving the superposed version of r_1 or xtm
	double** r1, ** r2;    // for Kabsch rotation

	const int a_opt = 0;
	const int mol_type = -2; // > 0 is RNA

	/***********************/
	/* allocate memory     */
	/***********************/
	int minlen = min(xlen, ylen);
	NewArray(&score, xlen + 1, ylen + 1);
	NewArray(&path, xlen + 1, ylen + 1);
	NewArray(&val, xlen + 1, ylen + 1);
	NewArray(&xtm, minlen, 3);
	NewArray(&ytm, minlen, 3);
	NewArray(&xt, xlen, 3);
	NewArray(&r1, minlen, 3);
	NewArray(&r2, minlen, 3);

	/***********************/
	/*    parameter set    */
	/***********************/
	parameter_set4search(xlen, ylen, D0_MIN, Lnorm,
		score_d8, d0, d0_search, dcu0);
	int simplify_step = 40; //for similified search engine
	int score_sum_method = 8;  //for scoring method, whether only sum over pairs with dis<score_d8

	int i;
	int* invmap0 = new int[ylen + 1];
	int* invmap = new int[ylen + 1];
	double TM, TMmax = -1;
	for (i = 0; i < ylen; i++) invmap0[i] = -1;

	double ddcc = 0.4;
	if (Lnorm <= 40) ddcc = 0.1;   //Lnorm was setted in parameter_set4search
	double local_d0_search = d0_search;

#if 0
	/******************************************************/
	/*    get initial alignment with gapless threading    */
	/******************************************************/
	get_initial(r1, r2, xtm, ytm, xa, ya, xlen, ylen, invmap0, d0,
		d0_search, t, u);

	TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap0,
		t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
		score_d8, d0);
	if (TM > TMmax) TMmax = TM;
	if (TMcut > 0) copy_t_u(t, u, t0, u0);

	//run dynamic programing iteratively to find the best alignment
	TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa, ya, xlen, ylen,
		t, u, invmap, 0, 2, 30, local_d0_search,
		D0_MIN, Lnorm, d0, score_d8);
	if (TM > TMmax)
		{
		TMmax = TM;
		for (int i = 0; i < ylen; i++) invmap0[i] = invmap[i];
		if (TMcut > 0) copy_t_u(t, u, t0, u0);
		}

	if (TMcut > 0) // pre-terminate if TM-score is too low
		{
		double TMtmp = approx_TM(xlen, ylen, a_opt,
			xa, ya, t0, u0, invmap0, mol_type);

		if (TMtmp < 0.5 * TMcut)
			{
			TM1 = TM2 = TM3 = TM4 = TM5 = TMtmp;
			clean_up_after_approx_TM(invmap0, invmap, score, path, val,
				xtm, ytm, xt, r1, r2, xlen, minlen);
			return 2;
			}
		}

	/************************************************************/
	/*    get initial alignment based on secondary structure    */
	/************************************************************/
	get_initial_ss(path, val, secx, secy, xlen, ylen, invmap);
	TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap,
		t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
		score_d8, d0);

	if (TM > TMmax)
		{
		TMmax = TM;
		for (int i = 0; i < ylen; i++) invmap0[i] = invmap[i];
		if (TMcut > 0) copy_t_u(t, u, t0, u0);
		}
	if (TM > TMmax * 0.2)
		{
		TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa, ya,
			xlen, ylen, t, u, invmap, 0, 2, 30,
			local_d0_search, D0_MIN, Lnorm, d0, score_d8);
		if (TM > TMmax)
			{
			TMmax = TM;
			for (int i = 0; i < ylen; i++) invmap0[i] = invmap[i];
			if (TMcut > 0) copy_t_u(t, u, t0, u0);
			}
		}

	if (TMcut > 0) // pre-terminate if TM-score is too low
		{
		double TMtmp = approx_TM(xlen, ylen, a_opt,
			xa, ya, t0, u0, invmap0, mol_type);

		if (TMtmp < 0.52 * TMcut)
			{
			TM1 = TM2 = TM3 = TM4 = TM5 = TMtmp;
			clean_up_after_approx_TM(invmap0, invmap, score, path, val,
				xtm, ytm, xt, r1, r2, xlen, minlen);
			return 3;
			}
		}

#endif // 0
	/************************************************************/
	/*    get initial alignment based on local superposition    */
	/************************************************************/
	//=initial5 in original TM-align
	bool Initial5Ok = get_initial5(r1, r2, xtm, ytm, path, val, xa, ya,
		xlen, ylen, invmap, d0, d0_search, D0_MIN);
	if (!Initial5Ok)
		return 2;

	TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
		invmap, t, u, simplify_step, score_sum_method,
		local_d0_search, Lnorm, score_d8, d0);
	if (TM > TMmax)
		{
		TMmax = TM;
		for (int i = 0; i < ylen; i++) invmap0[i] = invmap[i];
		}
	if (TM > TMmax * ddcc)
		{
		TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa, ya,
			xlen, ylen, t, u, invmap, 0, 2, 2, local_d0_search,
			D0_MIN, Lnorm, d0, score_d8);
		if (TM > TMmax)
			{
			TMmax = TM;
			for (int i = 0; i < ylen; i++) invmap0[i] = invmap[i];
			}
		}

#if 0
	/********************************************************************/
	/* get initial alignment by local superposition+secondary structure */
	/********************************************************************/
	//=initial3 in original TM-align
	get_initial_ssplus(r1, r2, score, path, val, secx, secy, xa, ya,
		xlen, ylen, invmap0, invmap, D0_MIN, d0);
	TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap,
		t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
		score_d8, d0);
	if (TM > TMmax)
		{
		TMmax = TM;
		for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
		if (TMcut > 0) copy_t_u(t, u, t0, u0);
		}
	if (TM > TMmax * ddcc)
		{
		TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa, ya,
			xlen, ylen, t, u, invmap, 0, 2, 30,
			local_d0_search, D0_MIN, Lnorm, d0, score_d8);
		if (TM > TMmax)
			{
			TMmax = TM;
			for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
			if (TMcut > 0) copy_t_u(t, u, t0, u0);
			}
		}

	if (TMcut > 0) // pre-terminate if TM-score is too low
		{
		double TMtmp = approx_TM(xlen, ylen, a_opt,
			xa, ya, t0, u0, invmap0, mol_type);

		if (TMtmp < 0.56 * TMcut)
			{
			TM1 = TM2 = TM3 = TM4 = TM5 = TMtmp;
			clean_up_after_approx_TM(invmap0, invmap, score, path, val,
				xtm, ytm, xt, r1, r2, xlen, minlen);
			return 5;
			}
		}

	/*******************************************************************/
	/*    get initial alignment based on fragment gapless threading    */
	/*******************************************************************/
	//=initial4 in original TM-align
	get_initial_fgt(r1, r2, xtm, ytm, xa, ya, xlen, ylen,
		invmap, d0, d0_search, dcu0, t, u);
	TM = detailed_search(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen, invmap,
		t, u, simplify_step, score_sum_method, local_d0_search, Lnorm,
		score_d8, d0);
	if (TM > TMmax)
		{
		TMmax = TM;
		for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
		if (TMcut > 0) copy_t_u(t, u, t0, u0);
		}
	if (TM > TMmax * ddcc)
		{
		TM = DP_iter(r1, r2, xtm, ytm, xt, path, val, xa, ya,
			xlen, ylen, t, u, invmap, 1, 2, 2, local_d0_search, D0_MIN,
			Lnorm, d0, score_d8);
		if (TM > TMmax)
			{
			TMmax = TM;
			for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
			if (TMcut > 0) copy_t_u(t, u, t0, u0);
			}
		}

	if (TMcut > 0) // pre-terminate if TM-score is too low
		{
		double TMtmp = approx_TM(xlen, ylen, a_opt,
			xa, ya, t0, u0, invmap0, mol_type);

		if (TMtmp < 0.58 * TMcut)
			{
			TM1 = TM2 = TM3 = TM4 = TM5 = TMtmp;
			clean_up_after_approx_TM(invmap0, invmap, score, path, val,
				xtm, ytm, xt, r1, r2, xlen, minlen);
			return 6;
			}
		}
#endif

	//*******************************************************************//
	//    The alignment will not be changed any more in the following    //
	//*******************************************************************//
	//check if the initial alignment is generated approriately
	bool flag = false;
	for (i = 0; i < ylen; i++)
		{
		if (invmap0[i] >= 0)
			{
			flag = true;
			break;
			}
		}
	if (!flag)
	// ========= NO ALIGNMENT FOUND =========
		return 1;

	//********************************************************************//
	//    Detailed TMscore search engine --> prepare for final TMscore    //
	//********************************************************************//
	//run detailed TMscore search engine for the best alignment, and
	//extract the best rotation matrix (t, u) for the best alginment
	simplify_step = 1;
	score_sum_method = 8;
	TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
		invmap0, t, u, simplify_step, score_sum_method, local_d0_search,
		false, Lnorm, score_d8, d0);

	//select pairs with dis<d8 for final TMscore computation and output alignment
	int k = 0;
	int* m1, * m2;
	double d;
	m1 = new int[xlen]; //alignd index in x
	m2 = new int[ylen]; //alignd index in y
	do_rotation(xa, xt, xlen, t, u);
	k = 0;

	for (int j = 0; j < ylen; j++)
		{
		i = invmap0[j];
		if (i >= 0)//aligned
			{
			n_ali++;
			d = sqrt(dist(&xt[i][0], &ya[j][0]));
			if (d <= score_d8)
				{
				m1[k] = i;
				m2[k] = j;

				xtm[k][0] = xa[i][0];
				xtm[k][1] = xa[i][1];
				xtm[k][2] = xa[i][2];

				ytm[k][0] = ya[j][0];
				ytm[k][1] = ya[j][1];
				ytm[k][2] = ya[j][2];

				r1[k][0] = xt[i][0];
				r1[k][1] = xt[i][1];
				r1[k][2] = xt[i][2];
				r2[k][0] = ya[j][0];
				r2[k][1] = ya[j][1];
				r2[k][2] = ya[j][2];

				k++;
				}
			}
		}
	n_ali8 = k;

	Kabsch(r1, r2, n_ali8, 0, &rmsd0, t, u);

// rmsd0 is used for final output, only recalculate rmsd0, not t & u
	rmsd0 = sqrt(rmsd0 / n_ali8);

	//****************************************//
	//              Final TMscore             //
	//    Please set parameters for output    //
	//****************************************//
	double rmsd;
	simplify_step = 1;
	score_sum_method = 0;
	double Lnorm_0 = ylen;

	//normalized by length of structure A
	parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
	d0A = d0;
	d0_0 = d0A;
	local_d0_search = d0_search;
	TM1 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0, simplify_step,
		score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0);
	TM_0 = TM1;
#if 1
	//normalized by length of structure B
	parameter_set4final(xlen + 0.0, D0_MIN, Lnorm, d0, d0_search, mol_type);
	d0B = d0;
	local_d0_search = d0_search;
	TM2 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t, u, simplify_step,
		score_sum_method, &rmsd, local_d0_search, Lnorm, score_d8, d0);
#endif
	/* derive alignment from superposition */
	int ali_len = xlen + ylen; //maximum length of alignment
	seqxA.assign(ali_len, '-');
	seqM.assign(ali_len, ' ');
	seqyA.assign(ali_len, '-');

	//do_rotation(xa, xt, xlen, t, u);
	do_rotation(xa, xt, xlen, t0, u0);

#if 0
	{
	Log("\n");
	Log("Rotation:\n");
	for (int i = 0; i < xlen; ++i)
		{
		Log("[%4u]  %6.2f  %6.2f  %6.2f  =>   %6.2f  %6.2f  %6.2f\n",
		  i, xa[i][0], xa[i][1], xa[i][2],
		  xt[i][0], xt[i][1], xt[i][2]);
		}
	}
#endif

	int kk = 0, i_old = 0, j_old = 0;
	d = 0;
	for (int k = 0; k < n_ali8; k++)
		{
		for (int i = i_old; i < m1[k]; i++)
			{
			//align x to gap
			seqxA[kk] = seqx[i];
			seqyA[kk] = '-';
			seqM[kk] = ' ';
			kk++;
			}

		for (int j = j_old; j < m2[k]; j++)
			{
			//align y to gap
			seqxA[kk] = '-';
			seqyA[kk] = seqy[j];
			seqM[kk] = ' ';
			kk++;
			}

		seqxA[kk] = seqx[m1[k]];
		seqyA[kk] = seqy[m2[k]];
		Liden += (seqxA[kk] == seqyA[kk]);
		d = sqrt(dist(&xt[m1[k]][0], &ya[m2[k]][0]));
		if (d < d0_out) seqM[kk] = ':';
		else         seqM[kk] = '.';
		kk++;
		i_old = m1[k] + 1;
		j_old = m2[k] + 1;
		}

	//tail
	for (int i = i_old; i < xlen; i++)
		{
		//align x to gap
		seqxA[kk] = seqx[i];
		seqyA[kk] = '-';
		seqM[kk] = ' ';
		kk++;
		}
	for (int j = j_old; j < ylen; j++)
		{
		//align y to gap
		seqxA[kk] = '-';
		seqyA[kk] = seqy[j];
		seqM[kk] = ' ';
		kk++;
		}
	seqxA = seqxA.substr(0, kk);
	seqyA = seqyA.substr(0, kk);
	seqM = seqM.substr(0, kk);

	/* free memory */
	clean_up_after_approx_TM(invmap0, invmap, score, path, val,
		xtm, ytm, xt, r1, r2, xlen, minlen);
	delete[] m1;
	delete[] m2;
	return 0; // zero for no exception
	}

void TMA::WriteAln(FILE *f,
	const int xlen, const int ylen,
	const double TM1, const double TM2,
	double d0A, double d0B,
	const string& seqM, const string& seqxA, const string& seqyA) const
	{
	if (f == 0)
		return;

	const uint ColCount = SIZE(seqxA);
	asserta(SIZE(seqyA) == ColCount);
	asserta(SIZE(seqM) == ColCount);

	fprintf(f, "\n");
	fprintf(m_faln, "____________________________\n");
	fprintf(f, ">%s\n", m_Q->m_Label.c_str());
	fprintf(f, ">%s\n", m_R->m_Label.c_str());
	fprintf(f, "\n");

	const char *Q = seqxA.c_str();
	const char *R = seqyA.c_str();
	const char *M = seqM.c_str();

	uint ColStart = 0;
	for (;;)
		{
		uint ColEnd = ColStart + 100;
		if (ColEnd >= ColCount)
			ColEnd = ColCount;
		asserta(ColEnd > ColStart);
		uint n = ColEnd - ColStart;

		fprintf(f, "Q  %*.*s\n", n, n, Q + ColStart);
		fprintf(f, "   %*.*s\n", n, n, M + ColStart);
		fprintf(f, "R  %*.*s\n", n, n, R + ColStart);
		fprintf(f, "\n");

		ColStart = ColEnd;
		if (ColStart >= ColCount)
			break;
		}

	fprintf(f, "\n");
//	fprintf(f, "TM-score= %5.3f (QL=%d, dQ=%.2f)\n", TM2, xlen, d0B);
	fprintf(f, "TM-score= %5.3f (RL=%d, dR=%.2f)\n", TM1, ylen, d0A);
	}

/**************************************************************************
Implemetation of Kabsch algoritm for finding the best rotation matrix
---------------------------------------------------------------------------
x     - x(i,m) are coordinates of atom m in set x            (input)
y     - y(i,m) are coordinates of atom m in set y            (input)
n     - n is number of atom pairs                            (input)
mode  - 0: calculate rms only                                (input)
  1:calculate u,t only                                (takes medium)
  2:calculate rms,u,t                                 (takes longer)
rms   - sum of w*(ux+t-y)**2 over all atom pairs            (output)
u    - u(i,j) is   rotation  matrix for best superposition  (output)
t    - t(i)   is translation vector for best superposition  (output)
**************************************************************************/
bool TMA::Kabsch(double** x, double** y, int n, int mode, double* rms,
	double t[3], double u[3][3])
	{
	int i, j, m, m1, l, k;
	double e0, rms1, d, h, g;
	double cth, sth, sqrth, p, det, sigma;
	double xc[3], yc[3];
	double a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
	double sqrt3 = 1.73205080756888, tol = 0.01;
	int ip[] = { 0, 1, 3, 1, 2, 4, 3, 4, 5 };
	int ip2312[] = { 1, 2, 0, 1 };

	int a_failed = 0, b_failed = 0;
	double epsilon = 0.00000001;

	//initializtation
	*rms = 0;
	rms1 = 0;
	e0 = 0;
	double c1[3], c2[3];
	double s1[3], s2[3];
	double sx[3], sy[3], sz[3];
	for (i = 0; i < 3; i++)
		{
		s1[i] = 0.0;
		s2[i] = 0.0;

		sx[i] = 0.0;
		sy[i] = 0.0;
		sz[i] = 0.0;
		}

	for (i = 0; i < 3; i++)
		{
		xc[i] = 0.0;
		yc[i] = 0.0;
		t[i] = 0.0;
		for (j = 0; j < 3; j++)
			{
			u[i][j] = 0.0;
			r[i][j] = 0.0;
			a[i][j] = 0.0;
			if (i == j)
				{
				u[i][j] = 1.0;
				a[i][j] = 1.0;
				}
			}
		}

	if (n < 1) return false;

	//compute centers for vector sets x, y
	for (i = 0; i < n; i++)
		{
		for (j = 0; j < 3; j++)
			{
			c1[j] = x[i][j];
			c2[j] = y[i][j];

			s1[j] += c1[j];
			s2[j] += c2[j];
			}

		for (j = 0; j < 3; j++)
			{
			sx[j] += c1[0] * c2[j];
			sy[j] += c1[1] * c2[j];
			sz[j] += c1[2] * c2[j];
			}
		}
	for (i = 0; i < 3; i++)
		{
		xc[i] = s1[i] / n;
		yc[i] = s2[i] / n;
		}
	if (mode == 2 || mode == 0)
		for (int mm = 0; mm < n; mm++)
			for (int nn = 0; nn < 3; nn++)
				e0 += (x[mm][nn] - xc[nn]) * (x[mm][nn] - xc[nn]) +
				(y[mm][nn] - yc[nn]) * (y[mm][nn] - yc[nn]);
	for (j = 0; j < 3; j++)
		{
		r[j][0] = sx[j] - s1[0] * s2[j] / n;
		r[j][1] = sy[j] - s1[1] * s2[j] / n;
		r[j][2] = sz[j] - s1[2] * s2[j] / n;
		}

	//compute determinat of matrix r
	det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])\
		- r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])\
		+ r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
	sigma = det;

	//compute tras(r)*r
	m = 0;
	for (j = 0; j < 3; j++)
		{
		for (i = 0; i <= j; i++)
			{
			rr[m] = r[0][i] * r[0][j] + r[1][i] * r[1][j] + r[2][i] * r[2][j];
			m++;
			}
		}

	double spur = (rr[0] + rr[2] + rr[5]) / 3.0;
	double cof = (((((rr[2] * rr[5] - rr[4] * rr[4]) + rr[0] * rr[5])\
		- rr[3] * rr[3]) + rr[0] * rr[2]) - rr[1] * rr[1]) / 3.0;
	det = det * det;

	for (i = 0; i < 3; i++) e[i] = spur;

	if (spur > 0)
		{
		d = spur * spur;
		h = d - cof;
		g = (spur * cof - det) / 2.0 - spur * h;

		if (h > 0)
			{
			sqrth = sqrt(h);
			d = h * h * h - g * g;
			if (d < 0.0) d = 0.0;
			d = atan2(sqrt(d), -g) / 3.0;
			cth = sqrth * cos(d);
			sth = sqrth * sqrt3 * sin(d);
			e[0] = (spur + cth) + cth;
			e[1] = (spur - cth) + sth;
			e[2] = (spur - cth) - sth;

			if (mode != 0)
				{//compute a                
				for (l = 0; l < 3; l = l + 2)
					{
					d = e[l];
					ss[0] = (d - rr[2]) * (d - rr[5]) - rr[4] * rr[4];
					ss[1] = (d - rr[5]) * rr[1] + rr[3] * rr[4];
					ss[2] = (d - rr[0]) * (d - rr[5]) - rr[3] * rr[3];
					ss[3] = (d - rr[2]) * rr[3] + rr[1] * rr[4];
					ss[4] = (d - rr[0]) * rr[4] + rr[1] * rr[3];
					ss[5] = (d - rr[0]) * (d - rr[2]) - rr[1] * rr[1];

					if (fabs(ss[0]) <= epsilon) ss[0] = 0.0;
					if (fabs(ss[1]) <= epsilon) ss[1] = 0.0;
					if (fabs(ss[2]) <= epsilon) ss[2] = 0.0;
					if (fabs(ss[3]) <= epsilon) ss[3] = 0.0;
					if (fabs(ss[4]) <= epsilon) ss[4] = 0.0;
					if (fabs(ss[5]) <= epsilon) ss[5] = 0.0;

					if (fabs(ss[0]) >= fabs(ss[2]))
						{
						j = 0;
						if (fabs(ss[0]) < fabs(ss[5])) j = 2;
						}
					else if (fabs(ss[2]) >= fabs(ss[5])) j = 1;
					else j = 2;

					d = 0.0;
					j = 3 * j;
					for (i = 0; i < 3; i++)
						{
						k = ip[i + j];
						a[i][l] = ss[k];
						d = d + ss[k] * ss[k];
						}


					//if( d > 0.0 ) d = 1.0 / sqrt(d);
					if (d > epsilon) d = 1.0 / sqrt(d);
					else d = 0.0;
					for (i = 0; i < 3; i++) a[i][l] = a[i][l] * d;
					}//for l

				d = a[0][0] * a[0][2] + a[1][0] * a[1][2] + a[2][0] * a[2][2];
				if ((e[0] - e[1]) > (e[1] - e[2]))
					{
					m1 = 2;
					m = 0;
					}
				else
					{
					m1 = 0;
					m = 2;
					}
				p = 0;
				for (i = 0; i < 3; i++)
					{
					a[i][m1] = a[i][m1] - d * a[i][m];
					p = p + a[i][m1] * a[i][m1];
					}
				if (p <= tol)
					{
					p = 1.0;
					for (i = 0; i < 3; i++)
						{
						if (p < fabs(a[i][m])) continue;
						p = fabs(a[i][m]);
						j = i;
						}
					k = ip2312[j];
					l = ip2312[j + 1];
					p = sqrt(a[k][m] * a[k][m] + a[l][m] * a[l][m]);
					if (p > tol)
						{
						a[j][m1] = 0.0;
						a[k][m1] = -a[l][m] / p;
						a[l][m1] = a[k][m] / p;
						}
					else a_failed = 1;
					}//if p<=tol
				else
					{
					p = 1.0 / sqrt(p);
					for (i = 0; i < 3; i++) a[i][m1] = a[i][m1] * p;
					}//else p<=tol  
				if (a_failed != 1)
					{
					a[0][1] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
					a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];
					a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
					}
				}//if(mode!=0)       
			}//h>0

			//compute b anyway
		if (mode != 0 && a_failed != 1)//a is computed correctly
			{
			//compute b
			for (l = 0; l < 2; l++)
				{
				d = 0.0;
				for (i = 0; i < 3; i++)
					{
					b[i][l] = r[i][0] * a[0][l] +
						r[i][1] * a[1][l] + r[i][2] * a[2][l];
					d = d + b[i][l] * b[i][l];
					}
				//if( d > 0 ) d = 1.0 / sqrt(d);
				if (d > epsilon) d = 1.0 / sqrt(d);
				else d = 0.0;
				for (i = 0; i < 3; i++) b[i][l] = b[i][l] * d;
				}
			d = b[0][0] * b[0][1] + b[1][0] * b[1][1] + b[2][0] * b[2][1];
			p = 0.0;

			for (i = 0; i < 3; i++)
				{
				b[i][1] = b[i][1] - d * b[i][0];
				p += b[i][1] * b[i][1];
				}

			if (p <= tol)
				{
				p = 1.0;
				for (i = 0; i < 3; i++)
					{
					if (p < fabs(b[i][0])) continue;
					p = fabs(b[i][0]);
					j = i;
					}
				k = ip2312[j];
				l = ip2312[j + 1];
				p = sqrt(b[k][0] * b[k][0] + b[l][0] * b[l][0]);
				if (p > tol)
					{
					b[j][1] = 0.0;
					b[k][1] = -b[l][0] / p;
					b[l][1] = b[k][0] / p;
					}
				else b_failed = 1;
				}//if( p <= tol )
			else
				{
				p = 1.0 / sqrt(p);
				for (i = 0; i < 3; i++) b[i][1] = b[i][1] * p;
				}
			if (b_failed != 1)
				{
				b[0][2] = b[1][0] * b[2][1] - b[1][1] * b[2][0];
				b[1][2] = b[2][0] * b[0][1] - b[2][1] * b[0][0];
				b[2][2] = b[0][0] * b[1][1] - b[0][1] * b[1][0];
				//compute u
				for (i = 0; i < 3; i++)
					for (j = 0; j < 3; j++)
						u[i][j] = b[i][0] * a[j][0] +
						b[i][1] * a[j][1] + b[i][2] * a[j][2];
				}

			//compute t
			for (i = 0; i < 3; i++)
				t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) -
				u[i][2] * xc[2];
			}//if(mode!=0 && a_failed!=1)
		}//spur>0
	else //just compute t and errors
		{
		//compute t
		for (i = 0; i < 3; i++)
			t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) -
			u[i][2] * xc[2];
		}//else spur>0 

		//compute rms
	for (i = 0; i < 3; i++)
		{
		if (e[i] < 0) e[i] = 0;
		e[i] = sqrt(e[i]);
		}
	d = e[2];
	if (sigma < 0.0) d = -d;
	d = (d + e[1]) + e[0];

	if (mode == 2 || mode == 0)
		{
		rms1 = (e0 - d) - d;
		if (rms1 < 0.0) rms1 = 0.0;
		}

	*rms = rms1;
	return true;
	}

/* Input: score[1:len1, 1:len2], and gap_open
 * Output: j2i[1:len2] \in {1:len1} U {-1}
 * path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical */
void TMA::NWDP_TM1(double** score, bool** path, double** val,
	int len1, int len2, double gap_open, int j2i[])
	{
	int i, j;
	double h, v, d;

	//initialization
	for (i = 0; i <= len1; i++)
		{
		val[i][0] = 0;
		//val[i][0]=i*gap_open;
		path[i][0] = false; //not from diagonal
		}

	for (j = 0; j <= len2; j++)
		{
		val[0][j] = 0;
		//val[0][j]=j*gap_open;
		path[0][j] = false; //not from diagonal
		j2i[j] = -1;    //all are not aligned, only use j2i[1:len2]
		}


	//decide matrix and path
	for (i = 1; i <= len1; i++)
		{
		for (j = 1; j <= len2; j++)
			{
			d = val[i - 1][j - 1] + score[i][j]; //diagonal

			//symbol insertion in horizontal (= a gap in vertical)
			h = val[i - 1][j];
			if (path[i - 1][j]) h += gap_open; //aligned in last position

			//symbol insertion in vertical
			v = val[i][j - 1];
			if (path[i][j - 1]) v += gap_open; //aligned in last position


			if (d >= h && d >= v)
				{
				path[i][j] = true; //from diagonal
				val[i][j] = d;
				}
			else
				{
				path[i][j] = false; //from horizontal
				if (v >= h) val[i][j] = v;
				else val[i][j] = h;
				}
			} //for i
		} //for j

		//trace back to extract the alignment
	i = len1;
	j = len2;
	while (i > 0 && j > 0)
		{
		if (path[i][j]) //from diagonal
			{
			j2i[j - 1] = i - 1;
			i--;
			j--;
			}
		else
			{
			h = val[i - 1][j];
			if (path[i - 1][j]) h += gap_open;

			v = val[i][j - 1];
			if (path[i][j - 1]) v += gap_open;

			if (v >= h) j--;
			else i--;
			}
		}
	}

/* Input: vectors x, y, rotation matrix t, u, scale factor d02, and gap_open
 * Output: j2i[1:len2] \in {1:len1} U {-1}
 * path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical */
void TMA::NWDP_TM2(bool** path, double** val, double** x, double** y,
	int len1, int len2, double t[3], double u[3][3],
	double d02, double gap_open, int j2i[])
	{
	int i, j;
	double h, v, d;

	//initialization. use old val[i][0] and val[0][j] initialization
	//to minimize difference from TMalign fortran version
	for (i = 0; i <= len1; i++)
		{
		val[i][0] = 0;
		//val[i][0]=i*gap_open;
		path[i][0] = false; //not from diagonal
		}

	for (j = 0; j <= len2; j++)
		{
		val[0][j] = 0;
		//val[0][j]=j*gap_open;
		path[0][j] = false; //not from diagonal
		j2i[j] = -1;    //all are not aligned, only use j2i[1:len2]
		}
	double xx[3], dij;


	//decide matrix and path
	for (i = 1; i <= len1; i++)
		{
		transform(t, u, &x[i - 1][0], xx);
		for (j = 1; j <= len2; j++)
			{
			dij = dist(xx, &y[j - 1][0]);
			d = val[i - 1][j - 1] + 1.0 / (1 + dij / d02);

			//symbol insertion in horizontal (= a gap in vertical)
			h = val[i - 1][j];
			if (path[i - 1][j]) h += gap_open; //aligned in last position

			//symbol insertion in vertical
			v = val[i][j - 1];
			if (path[i][j - 1]) v += gap_open; //aligned in last position


			if (d >= h && d >= v)
				{
				path[i][j] = true; //from diagonal
				val[i][j] = d;
				}
			else
				{
				path[i][j] = false; //from horizontal
				if (v >= h) val[i][j] = v;
				else val[i][j] = h;
				}
			} //for i
		} //for j

		//trace back to extract the alignment
	i = len1;
	j = len2;
	while (i > 0 && j > 0)
		{
		if (path[i][j]) //from diagonal
			{
			j2i[j - 1] = i - 1;
			i--;
			j--;
			}
		else
			{
			h = val[i - 1][j];
			if (path[i - 1][j]) h += gap_open;

			v = val[i][j - 1];
			if (path[i][j - 1]) v += gap_open;

			if (v >= h) j--;
			else i--;
			}
		}
	}

/* This is the same as the previous NWDP_TM, except for the lack of rotation
 * Input: vectors x, y, scale factor d02, and gap_open
 * Output: j2i[1:len2] \in {1:len1} U {-1}
 * path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical */
void TMA::NWDP_SE(bool** path, double** val, double** x, double** y,
	int len1, int len2, double d02, double gap_open, int j2i[])
	{
	int i, j;
	double h, v, d;

	for (i = 0; i <= len1; i++)
		{
		val[i][0] = 0;
		path[i][0] = false; //not from diagonal
		}

	for (j = 0; j <= len2; j++)
		{
		val[0][j] = 0;
		path[0][j] = false; //not from diagonal
		j2i[j] = -1;    //all are not aligned, only use j2i[1:len2]
		}
	double dij;

	//decide matrix and path
	for (i = 1; i <= len1; i++)
		{
		for (j = 1; j <= len2; j++)
			{
			dij = dist(&x[i - 1][0], &y[j - 1][0]);
			d = val[i - 1][j - 1] + 1.0 / (1 + dij / d02);

			//symbol insertion in horizontal (= a gap in vertical)
			h = val[i - 1][j];
			if (path[i - 1][j]) h += gap_open; //aligned in last position

			//symbol insertion in vertical
			v = val[i][j - 1];
			if (path[i][j - 1]) v += gap_open; //aligned in last position


			if (d >= h && d >= v)
				{
				path[i][j] = true; //from diagonal
				val[i][j] = d;
				}
			else
				{
				path[i][j] = false; //from horizontal
				if (v >= h) val[i][j] = v;
				else val[i][j] = h;
				}
			} //for i
		} //for j

		//trace back to extract the alignment
	i = len1;
	j = len2;
	while (i > 0 && j > 0)
		{
		if (path[i][j]) //from diagonal
			{
			j2i[j - 1] = i - 1;
			i--;
			j--;
			}
		else
			{
			h = val[i - 1][j];
			if (path[i - 1][j]) h += gap_open;

			v = val[i][j - 1];
			if (path[i][j - 1]) v += gap_open;

			if (v >= h) j--;
			else i--;
			}
		}
	}

/* +ss
 * Input: secondary structure secx, secy, and gap_open
 * Output: j2i[1:len2] \in {1:len1} U {-1}
 * path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical */
void TMA::NWDP_TM3(bool** path, double** val, const char* secx, const char* secy,
	const int len1, const int len2, const double gap_open, int j2i[])
	{
	int i, j;
	double h, v, d;

	//initialization
	for (i = 0; i <= len1; i++)
		{
		val[i][0] = 0;
		//val[i][0]=i*gap_open;
		path[i][0] = false; //not from diagonal
		}

	for (j = 0; j <= len2; j++)
		{
		val[0][j] = 0;
		//val[0][j]=j*gap_open;
		path[0][j] = false; //not from diagonal
		j2i[j] = -1;    //all are not aligned, only use j2i[1:len2]
		}

	//decide matrix and path
	for (i = 1; i <= len1; i++)
		{
		for (j = 1; j <= len2; j++)
			{
			d = val[i - 1][j - 1] + 1.0 * (secx[i - 1] == secy[j - 1]);

			//symbol insertion in horizontal (= a gap in vertical)
			h = val[i - 1][j];
			if (path[i - 1][j]) h += gap_open; //aligned in last position

			//symbol insertion in vertical
			v = val[i][j - 1];
			if (path[i][j - 1]) v += gap_open; //aligned in last position

			if (d >= h && d >= v)
				{
				path[i][j] = true; //from diagonal
				val[i][j] = d;
				}
			else
				{
				path[i][j] = false; //from horizontal
				if (v >= h) val[i][j] = v;
				else val[i][j] = h;
				}
			} //for i
		} //for j

		//trace back to extract the alignment
	i = len1;
	j = len2;
	while (i > 0 && j > 0)
		{
		if (path[i][j]) //from diagonal
			{
			j2i[j - 1] = i - 1;
			i--;
			j--;
			}
		else
			{
			h = val[i - 1][j];
			if (path[i - 1][j]) h += gap_open;

			v = val[i][j - 1];
			if (path[i][j - 1]) v += gap_open;

			if (v >= h) j--;
			else i--;
			}
		}
	}

void TMA::parameter_set4search(const int xlen, const int ylen,
	double& D0_MIN, double& Lnorm,
	double& score_d8, double& d0, double& d0_search, double& dcu0)
	{
	//parameter initilization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
	D0_MIN = 0.5;
	dcu0 = 4.25;                       //update 3.85-->4.25

	Lnorm = min(xlen, ylen);        //normaliz TMscore by this in searching
	if (Lnorm <= 19)                    //update 15-->19
		d0 = 0.168;                   //update 0.5-->0.168
	else d0 = (1.24 * pow((Lnorm * 1.0 - 15), 1.0 / 3) - 1.8);
	D0_MIN = d0 + 0.8;              //this should be moved to above
	d0 = D0_MIN;                  //update: best for search    

	d0_search = d0;
	if (d0_search > 8)   d0_search = 8;
	if (d0_search < 4.5) d0_search = 4.5;

	score_d8 = 1.5 * pow(Lnorm * 1.0, 0.3) + 3.5; //remove pairs with dis>d8 during search & final
	}

void TMA::parameter_set4final(const double len, double& D0_MIN, double& Lnorm,
	double& d0, double& d0_search, const int mol_type)
	{
	D0_MIN = 0.5;

	Lnorm = len;            //normaliz TMscore by this in searching
	if (Lnorm <= 21) d0 = 0.5;
	else d0 = (1.24 * pow((Lnorm * 1.0 - 15), 1.0 / 3) - 1.8);
	if (d0 < D0_MIN) d0 = D0_MIN;
	d0_search = d0;
	if (d0_search > 8)   d0_search = 8;
	if (d0_search < 4.5) d0_search = 4.5;
	}

void TMA::parameter_set4scale(const int len, const double d_s, double& Lnorm,
	double& d0, double& d0_search)
	{
	d0 = d_s;
	Lnorm = len;            //normaliz TMscore by this in searching
	d0_search = d0;
	if (d0_search > 8)   d0_search = 8;
	if (d0_search < 4.5) d0_search = 4.5;
	}

//     1, collect those residues with dis<d;
//     2, calculate TMscore
int TMA::score_fun8( double **xa, double **ya, int n_ali, double d, int i_ali[],
    double *score1, int score_sum_method, const double Lnorm, 
    const double score_d8, const double d0)
{
    double score_sum=0, di;
    double d_tmp=d*d;
    double d02=d0*d0;
    double score_d8_cut = score_d8*score_d8;
    
    int i, n_cut, inc=0;

#if 0
	if (g_Trace)
		Log("score_fun8()\n");
#endif
    while(1)
    {
        n_cut=0;
        score_sum=0;
        for(i=0; i<n_ali; i++)
        {
            di = dist(xa[i], ya[i]);
            if(di<d_tmp)
            {
                i_ali[n_cut]=i;
                n_cut++;
            }
			else
				{
#if 0
				if (g_Trace)
					Log("Discard dist di=%.3g\n", di);
#endif
				}
            if(score_sum_method==8)
            {                
                if(di<=score_d8_cut) score_sum += 1/(1+di/d02);
            }
            else score_sum += 1/(1+di/d02);
        }
        //there are not enough feasible pairs, reliefe the threshold         
        if(n_cut<3 && n_ali>3)
        {
            inc++;
            double dinc=(d+inc*0.5);
            d_tmp = dinc * dinc;
        }
        else break;
    }  

    *score1=score_sum/Lnorm;
    return n_cut;
}

int TMA::score_fun8_standard(double** xa, double** ya, int n_ali, double d,
	int i_ali[], double* score1, int score_sum_method,
	double score_d8, double d0)
	{
	double score_sum = 0, di;
	double d_tmp = d * d;
	double d02 = d0 * d0;
	double score_d8_cut = score_d8 * score_d8;

	int i, n_cut, inc = 0;
	while (1)
		{
		n_cut = 0;
		score_sum = 0;
		for (i = 0; i < n_ali; i++)
			{
			di = dist(xa[i], ya[i]);
			if (di < d_tmp)
				{
				i_ali[n_cut] = i;
				n_cut++;
				}
			if (score_sum_method == 8)
				{
				if (di <= score_d8_cut) score_sum += 1 / (1 + di / d02);
				}
			else
				{
				score_sum += 1 / (1 + di / d02);
				}
			}
		//there are not enough feasible pairs, reliefe the threshold         
		if (n_cut < 3 && n_ali>3)
			{
			inc++;
			double dinc = (d + inc * 0.5);
			d_tmp = dinc * dinc;
			}
		else break;
		}

	*score1 = score_sum / n_ali;
	return n_cut;
	}

uint TMA::ReadCal(const string &FileName, char *Seq, double **a)
	{
	FILE *f = OpenStdioFile(FileName);
	string LabelLine;
	bool Ok = ReadLineStdioFile(f, LabelLine);
	asserta(Ok);

	string Line;
	vector<string> Fields;
	uint n = 0;
	while (ReadLineStdioFile(f, Line))
		{
	/***
	>102l
	M       43.619  -1.924  8.869
	N       40.445  -0.876  10.670
	I       38.254  2.240   11.220
	F       40.340  3.621   14.036
	***/
		Split(Line, Fields, '\t');
		if (Fields.size() != 4 || Fields[0].size() != 1)
			Die("Invalid .cal record '%s'", Line.c_str());

		char aa = Fields[0][0];
		Seq[n] = aa;

		double X = StrToFloat(Fields[1]);
		double Y = StrToFloat(Fields[2]);
		double Z = StrToFloat(Fields[3]);

		a[n][0] = X;
		a[n][1] = Y;
		a[n][2] = Z;
		++n;
		}
	return n;
	}

double TMA::AlignChains(const PDBChain &Q, const PDBChain &R)
	{
	m_Q = &Q;
	m_R = &R;

	const uint QL = Q.GetSeqLength();
	const uint RL = R.GetSeqLength();

	char *seqx = new char[QL+1];
	char *seqy = new char[RL+1];

	memcpy(seqx, Q.m_Seq.c_str(), QL+1);
	memcpy(seqy, R.m_Seq.c_str(), RL+1);

	double **xa = 0;
	double **ya = 0;
	NewArray(&xa, QL, 3);
	NewArray(&ya, RL, 3);

	int xlen = (int) QL;
	int ylen = (int) RL;

	for (int i = 0; i < xlen; ++i)
		{
		xa[i][0] = Q.m_Xs[i];
		xa[i][1] = Q.m_Ys[i];
		xa[i][2] = Q.m_Zs[i];
		}

	for (int i = 0; i < ylen; ++i)
		{
		ya[i][0] = R.m_Xs[i];
		ya[i][1] = R.m_Ys[i];
		ya[i][2] = R.m_Zs[i];
		}

	double d0A = DBL_MAX;
	double d0B = DBL_MAX;
	string seqM;
	int iResult = TMalign_main(xa, ya, seqx, seqy, m_TM1, m_TM2, d0A, d0B,
	  seqM, m_QRow, m_RRow, xlen, ylen);
	m_Lock.lock();
	{
	if (iResult == 0)
		{
		WriteAln(m_faln, xlen, ylen, m_TM1, m_TM2, d0A, d0B, seqM, m_QRow, m_RRow);

		SeqToFasta(m_ffasta2, m_Q->m_Label, m_QRow);
		SeqToFasta(m_ffasta2, m_R->m_Label, m_RRow);
		}
	else
		{
		if (m_faln != 0)
			{
			fprintf(m_faln, "\n");
			fprintf(m_faln, "____________________________\n");
			fprintf(m_faln, ">%s\n", m_Q->m_Label.c_str());
			fprintf(m_faln, "no alignment to >%s\n", m_R->m_Label.c_str());
			}
		}
	m_Lock.unlock();
	}

	DeleteArray(&xa, QL);
	DeleteArray(&ya, RL);
	delete[] seqx;
	delete[] seqy;
	if (iResult != 0)
		return 0;
	return m_TM1;
	}
int TMA::TMalign_main_score(
	const string &rowa, const string &rowb,
	double** xa, double** ya,
	const char* seqx, const char* secy,
	double& TM1, double& TM2,
	double& d0A, double& d0B,
	string& seqM, string& seqxA, string& seqyA,
	const int xlen, const int ylen)
	{
	double D0_MIN;        //for d0
	double Lnorm;         //normalization length
	double score_d8, d0, d0_search, dcu0;//for TMscore search
	double t[3], u[3][3]; //Kabsch translation vector and rotation matrix
	double** xtm, ** ytm;  // for TMscore search engine
	double** xt;          //for saving the superposed version of r_1 or xtm
	double** r1, ** r2;    // for Kabsch rotation

	/***********************/
	/* allocate memory     */
	/***********************/
	int minlen = min(xlen, ylen);
	//NewArray(&val, xlen+1, ylen+1);
	NewArray(&xtm, minlen, 3);
	NewArray(&ytm, minlen, 3);
	NewArray(&xt, xlen, 3);
	NewArray(&r1, minlen, 3);
	NewArray(&r2, minlen, 3);

	/***********************/
	/*    parameter set    */
	/***********************/
	parameter_set4search(xlen, ylen, D0_MIN, Lnorm,
		score_d8, d0, d0_search, dcu0);
	if (0)
		{
		Log("parameter_set4search:\n");
		Log("   D0_MIN  %.3g\n", D0_MIN);
		Log("   Lnorm   %.3g\n", Lnorm);
		Log("score_d8   %.3g\n", score_d8);
		Log("      d0   %.3g\n", d0);
		Log("d0_search  %.3g\n", d0_search);
		Log("    dcu0   %.3g\n", dcu0);
		}

	int i;
	int* invmap0 = new int[ylen + 1];
	int* invmap = new int[ylen + 1];
	double TM, TMmax = -1;
	for (i = 0; i < ylen; i++) invmap0[i] = -1;

	double ddcc = 0.4;
	if (Lnorm <= 40) ddcc = 0.1;   //Lnorm was setted in parameter_set4search
	double local_d0_search = d0_search;

	for (int j = 0; j < ylen; j++)// Set aligned position to be "-1"
		invmap[j] = -1;

	int i1 = -1;// in C version, index starts from zero, not from one
	int i2 = -1;
	int L1 = (int) rowa.size();
	int L2 = (int) rowb.size();
	int L = min(L1, L2);// Get positions for aligned residues
	for (int kk1 = 0; kk1 < L; kk1++)
		{
		if (rowa[kk1] != '-') i1++;
		if (rowb[kk1] != '-')
			{
			i2++;
			if (i2 >= ylen || i1 >= xlen) kk1 = L;
			else if (rowa[kk1] != '-') invmap[i2] = i1;
			}
		}

	//--------------- 2. Align proteins from original alignment
	double prevD0_MIN = D0_MIN;// stored for later use
	double prevLnorm = Lnorm;
	double prevd0 = d0;
	double rmsd_ali = DBL_MAX;
	int L_ali = INT_MAX;
	double TM_ali = standard_TMscore(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
		invmap, L_ali, rmsd_ali, D0_MIN, Lnorm, d0, d0_search, score_d8,
		t, u, -2);
	D0_MIN = prevD0_MIN;
	Lnorm = prevLnorm;
	d0 = prevd0;
	int simplify_step = 40;
	int score_sum_method = 8;
	TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
		invmap, t, u, simplify_step, score_sum_method,
		local_d0_search, true, Lnorm, score_d8, d0);
	if (TM > TMmax)
		{
		TMmax = TM;
		for (i = 0; i < ylen; i++) invmap0[i] = invmap[i];
		}

	//check if the initial alignment is generated approriately
	bool flag = false;
	for (i = 0; i < ylen; i++)
		{
		if (invmap0[i] >= 0)
			{
			flag = true;
			break;
			}
		}
	if (!flag)
		Die("No alignment");

	simplify_step = 40;
	score_sum_method = 8;
	TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa, ya, xlen, ylen,
		invmap0, t, u, simplify_step, score_sum_method, local_d0_search,
		false, Lnorm, score_d8, d0);

	//select pairs with dis<d8 for final TMscore computation and output alignment
	int k = 0;
	int* m1, * m2;
	double d;
	m1 = new int[xlen]; //alignd index in x
	m2 = new int[ylen]; //alignd index in y
	do_rotation(xa, xt, xlen, t, u);
	k = 0;
	int n_ali = 0;
	for (int j = 0; j < ylen; j++)
		{
		i = invmap0[j];
		if (i >= 0)//aligned
			{
			n_ali++;
			d = sqrt(dist(&xt[i][0], &ya[j][0]));
			m1[k] = i;
			m2[k] = j;

			xtm[k][0] = xa[i][0];
			xtm[k][1] = xa[i][1];
			xtm[k][2] = xa[i][2];

			ytm[k][0] = ya[j][0];
			ytm[k][1] = ya[j][1];
			ytm[k][2] = ya[j][2];

			r1[k][0] = xt[i][0];
			r1[k][1] = xt[i][1];
			r1[k][2] = xt[i][2];
			r2[k][0] = ya[j][0];
			r2[k][1] = ya[j][1];
			r2[k][2] = ya[j][2];

			k++;
			}
		}
	int n_ali8 = k;

	double rmsd;
	double Lnorm_0 = ylen;

	//normalized by length of structure A
	parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, -2);
	d0A = d0;
	double d0_0 = d0A;
	local_d0_search = d0_search;
	double t0[3], u0[3][3];
	{
	const int simplify_step = 1;
	const int score_sum_method = 0;
	TM1 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
		simplify_step, score_sum_method, &rmsd, local_d0_search,
		Lnorm, score_d8, d0);
	}
	double TM_0 = TM1;

	//normalized by length of structure B
	parameter_set4final(xlen + 0.0, D0_MIN, Lnorm, d0, d0_search, -2);
	d0B = d0;
	local_d0_search = d0_search;
	{
	const int simplify_step = 1;
	const int score_sum_method = 0;
	TM2 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t, u,
		simplify_step, score_sum_method, &rmsd, local_d0_search,
		Lnorm, score_d8, d0);
	}

	return 0; // zero for no exception
	}

double TMA::CalcTMScore(const PDBChain &Q, const PDBChain &R,
  const string &RowQ, const string &RowR)
	{
	m_Q = &Q;
	m_R = &R;

	const uint QL = Q.GetSeqLength();
	const uint RL = R.GetSeqLength();

	char *seqx = new char[QL+1];
	char *seqy = new char[RL+1];

	memcpy(seqx, Q.m_Seq.c_str(), QL+1);
	memcpy(seqy, R.m_Seq.c_str(), RL+1);

	double **xa = 0;
	double **ya = 0;
	NewArray(&xa, QL, 3);
	NewArray(&ya, RL, 3);

	int xlen = (int) QL;
	int ylen = (int) RL;

	for (int i = 0; i < xlen; ++i)
		{
		xa[i][0] = Q.m_Xs[i];
		xa[i][1] = Q.m_Ys[i];
		xa[i][2] = Q.m_Zs[i];
		}

	for (int i = 0; i < ylen; ++i)
		{
		ya[i][0] = R.m_Xs[i];
		ya[i][1] = R.m_Ys[i];
		ya[i][2] = R.m_Zs[i];
		}

	double d0A = DBL_MAX;
	double d0B = DBL_MAX;
	string seqM;
	int iResult = TMalign_main_score(
		RowQ, RowR,
	  xa, ya, seqx, seqy, m_TM1, m_TM2, d0A, d0B,
	  seqM, m_QRow, m_RRow, xlen, ylen);

	return m_TM1;
	}
