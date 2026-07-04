#include "myutils.h"
#include "hitdata.h"
#include "flat_chain.h"
#include "kabsch.h"
#include "abcxyz.h"
#include <cmath>
#include <algorithm>

uint hitdata::m_tm_iter_max_iters = 4;
uint hitdata::m_tm_iter_min_pairs = 20;
double hitdata::m_tm_iter_min_frac = 0.5;
double hitdata::m_tm_iter_min_improve = 0.001;
double hitdata::m_tm_iter_d0_mult = 1.0;

static double tm_d0(double L)
	{
	if (L <= 15)
		return 0.5;
	if (L <= 21)
		return 0.494*pow(L - 15, 1.0/3.0) - 0.276;
	return 1.24*pow(L - 15, 1.0/3.0) - 1.8;
	}

static double tm_mean_L(const flat_chain_t *query, const flat_chain_t *target)
	{
	return (double(query->m_L) + double(target->m_L))/2.0;
	}

static uint tm_count_path_matches(const char *path, uint ncol)
	{
	uint M = 0;
	for (uint col = 0; col < ncol; ++col)
		if (path[col] == 'M')
			++M;
	return M;
	}

static uint tm_fill_path_pairs(
	hitdata &hit,
	uint M,
	double *x_data,
	double *y_data)
	{
	uint posQ = hit.qlo;
	uint posT = hit.tlo;
	uint m = 0;
	for (uint col = 0; col < hit.ncol; ++col)
		{
		const char c = hit.path[col];
		if (c == 'M')
			{
			asserta(posQ < hit.query->m_L);
			asserta(posT < hit.target->m_L);
			double *xp = x_data + m*3;
			double *yp = y_data + m*3;
			float xq, yq, zq;
			float xt, yt, zt;
			hit.query->get_coords(posQ, xq, yq, zq);
			hit.target->get_coords(posT, xt, yt, zt);
			xp[0] = xq;
			xp[1] = yq;
			xp[2] = zq;
			yp[0] = xt;
			yp[1] = yt;
			yp[2] = zt;
			++m;
			++posQ;
			++posT;
			}
		else if (c == 'D')
			++posT;
		else if (c == 'I')
			++posQ;
		else
			asserta(false);
		}
	asserta(m == M);
	return M;
	}

static uint tm_count_kept(uint M, const uint8_t *keep)
	{
	uint n = 0;
	for (uint i = 0; i < M; ++i)
		if (keep[i])
			++n;
	return n;
	}

static bool tm_keep_equal(uint M, const uint8_t *a, const uint8_t *b)
	{
	for (uint i = 0; i < M; ++i)
		if (a[i] != b[i])
			return false;
	return true;
	}

static void tm_kabsch_kept(
	hitdata &hit,
	uint M,
	const uint8_t *keep,
	const double *x_data,
	const double *y_data,
	uint nkeep,
	double t[3],
	double u[3][3])
	{
	uint k = 0;
	for (uint i = 0; i < M; ++i)
		{
		if (!keep[i])
			continue;
		hit.tm_x_ptrs_buf[k] = (double *) (x_data + i*3);
		hit.tm_y_ptrs_buf[k] = (double *) (y_data + i*3);
		++k;
		}
	asserta(k == nkeep);
	Kabsch(hit.tm_x_ptrs_buf, hit.tm_y_ptrs_buf, int(nkeep), t, u);
	}

static double tm_score_kept(
	uint M,
	const uint8_t *keep,
	const double *x_data,
	const double *y_data,
	const double t[3],
	const double u[3][3],
	double inv_d0_sq,
	double L)
	{
	double score_sum = 0;
	for (uint i = 0; i < M; ++i)
		{
		if (!keep[i])
			continue;
		double xt[3];
		const double *xi = x_data + i*3;
		const double *yi = y_data + i*3;
		transform(t, u, xi, xt);
		const double dx = xt[0] - yi[0];
		const double dy = xt[1] - yi[1];
		const double dz = xt[2] - yi[2];
		const double di_sq = dx*dx + dy*dy + dz*dz;
		score_sum += 1.0/(1.0 + di_sq*inv_d0_sq);
		}
	return score_sum/L;
	}

static void tm_distances_kept(
	uint M,
	const uint8_t *keep,
	const double *x_data,
	const double *y_data,
	const double t[3],
	const double u[3][3],
	double *dist_buf)
	{
	for (uint i = 0; i < M; ++i)
		{
		if (!keep[i])
			continue;
		double xt[3];
		const double *xi = x_data + i*3;
		const double *yi = y_data + i*3;
		transform(t, u, xi, xt);
		const double dx = xt[0] - yi[0];
		const double dy = xt[1] - yi[1];
		const double dz = xt[2] - yi[2];
		dist_buf[i] = sqrt(dx*dx + dy*dy + dz*dz);
		}
	}

struct tm_pair_dist
	{
	uint idx;
	double dist;
	};

static bool tm_pair_dist_less(const tm_pair_dist &a, const tm_pair_dist &b)
	{
	return a.dist < b.dist;
	}

static void tm_trim_keep(
	uint M,
	const uint8_t *keep,
	const double *dist_buf,
	double dist_thresh,
	uint min_keep,
	uint8_t *new_keep)
	{
	for (uint i = 0; i < M; ++i)
		new_keep[i] = 0;

	tm_pair_dist *pairs = (tm_pair_dist *) myalloc(tm_pair_dist, M);
	uint np = 0;
	for (uint i = 0; i < M; ++i)
		{
		if (!keep[i])
			continue;
		tm_pair_dist &pd = pairs[np++];
		pd.idx = i;
		pd.dist = dist_buf[i];
		if (dist_buf[i] <= dist_thresh)
			new_keep[i] = 1;
		}

	const uint nnew = tm_count_kept(M, new_keep);
	if (nnew >= min_keep)
		{
		myfree(pairs);
		return;
		}

	if (np == 0)
		{
		myfree(pairs);
		return;
		}

	uint mk = min_keep;
	if (mk > np)
		mk = np;

	std::sort(pairs, pairs + np, tm_pair_dist_less);
	for (uint k = 0; k < mk; ++k)
		new_keep[pairs[k].idx] = 1;
	myfree(pairs);
	}

double hitdata::calc_tm()
	{
	asserta(query);
	asserta(target);
	asserta(query->m_xyz);
	asserta(target->m_xyz);
	if (path == 0 || ncol == 0)
		return 0;

	const uint M = tm_count_path_matches(path, ncol);
	if (M == 0)
		return 0;

	tm_coord_ensure(M);
	double *x_data = tm_x_data(M);
	double *y_data = tm_y_data(M);
	tm_fill_path_pairs(*this, M, x_data, y_data);

	double t[3];
	double u[3][3];
	for (uint i = 0; i < M; ++i)
		{
		tm_x_ptrs_buf[i] = x_data + i*3;
		tm_y_ptrs_buf[i] = y_data + i*3;
		}
	Kabsch(tm_x_ptrs_buf, tm_y_ptrs_buf, int(M), t, u);

	const double L = tm_mean_L(query, target);
	const double d0 = tm_d0(L);
	const double inv_d0_sq = 1.0/(d0*d0);

	for (uint i = 0; i < M; ++i)
		tm_keep_buf[i] = 1;
	return tm_score_kept(M, tm_keep_buf, x_data, y_data, t, u, inv_d0_sq, L);
	}

double hitdata::calc_tm_iterate()
	{
	asserta(query);
	asserta(target);
	asserta(query->m_xyz);
	asserta(target->m_xyz);
	if (path == 0 || ncol == 0)
		return 0;

	const uint M = tm_count_path_matches(path, ncol);
	if (M < 3)
		return 0;

	tm_coord_ensure(M);
	double *x_data = tm_x_data(M);
	double *y_data = tm_y_data(M);
	tm_fill_path_pairs(*this, M, x_data, y_data);

	const double L = tm_mean_L(query, target);
	const double d0 = tm_d0(L);
	const double dist_thresh = d0*m_tm_iter_d0_mult;
	const double inv_d0_sq = 1.0/(d0*d0);

	uint min_keep = m_tm_iter_min_pairs;
	const uint min_frac_keep = uint(ceil(m_tm_iter_min_frac*double(M)));
	if (min_frac_keep > min_keep)
		min_keep = min_frac_keep;
	if (min_keep > M)
		min_keep = M;
	if (min_keep < 3)
		min_keep = 3;

	for (uint i = 0; i < M; ++i)
		tm_keep_buf[i] = 1;

	double t[3];
	double u[3][3];
	double tm = 0;
	double prev_tm = -1;
	for (uint iter = 0; iter < m_tm_iter_max_iters; ++iter)
		{
		const uint nkeep = tm_count_kept(M, tm_keep_buf);
		if (nkeep < 3)
			return 0;

		tm_kabsch_kept(*this, M, tm_keep_buf, x_data, y_data, nkeep, t, u);
		tm = tm_score_kept(M, tm_keep_buf, x_data, y_data, t, u, inv_d0_sq, L);

		if (iter > 0 && tm - prev_tm < m_tm_iter_min_improve)
			break;

		tm_distances_kept(M, tm_keep_buf, x_data, y_data, t, u, tm_dist_buf);

		for (uint i = 0; i < M; ++i)
			tm_new_keep_buf[i] = 0;
		tm_trim_keep(M, tm_keep_buf, tm_dist_buf, dist_thresh,
			min_keep, tm_new_keep_buf);

		if (tm_keep_equal(M, tm_keep_buf, tm_new_keep_buf))
			break;

		memcpy(tm_keep_buf, tm_new_keep_buf, M);
		prev_tm = tm;
		}

	const uint nfinal = tm_count_kept(M, tm_keep_buf);
	if (nfinal >= 3)
		{
		tm_kabsch_kept(*this, M, tm_keep_buf, x_data, y_data, nfinal, t, u);
		tm = tm_score_kept(M, tm_keep_buf, x_data, y_data, t, u, inv_d0_sq, L);
		}

	return tm;
	}
