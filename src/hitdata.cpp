#include "myutils.h"
#include "flat_chain.h"
#include "flat_params.h"
#include "hitdata.h"

void hitdata::cigar_free()
	{
	if (cigar_heap != 0)
		{
		myfree(cigar_heap);
		cigar_heap = 0;
		cigar_heap_cap = 0;
		}
	}

bool hitdata::cigar_ensure(uint need)
	{
	const uint used = cigar_len;
	const uint cap = cigar_heap ? cigar_heap_cap : CIGAR_BUFSIZE;
	if (used + need <= cap)
		return true;

	const uint new_cap = used + need;
	char *nb = (char *) myalloc(char, new_cap);
	const char *src = cigar_heap ? cigar_heap : cigar_buf;
	memcpy(nb, src, used);
	cigar_free();
	cigar_heap = nb;
	cigar_heap_cap = new_cap;
	return true;
	}

void hitdata::cigar_put_uint(uint n)
	{
	if (n < 10)
		{
		cigar_ensure(1);
		char *p = cigar_heap ? cigar_heap : cigar_buf;
		p[cigar_len++] = char('0' + n);
		return;
		}

	char tmp[10];
	char *t = tmp;
	do
		{
		*t++ = char('0' + (n % 10));
		n /= 10;
		}
	while (n > 0);
	const uint nd = uint(t - tmp);
	cigar_ensure(nd);
	char *p = cigar_heap ? cigar_heap : cigar_buf;
	p += cigar_len;
	while (t > tmp)
		*p++ = *--t;
	cigar_len += nd;
	}

void hitdata::cigar_put_op(uint n, char op)
	{
	cigar_put_uint(n);
	cigar_ensure(1);
	char *p = cigar_heap ? cigar_heap : cigar_buf;
	p[cigar_len++] = op;
	}

void hitdata::tm_coord_free()
	{
	if (tm_x_heap != 0)
		{
		myfree(tm_x_heap);
		tm_x_heap = 0;
		}
	if (tm_y_heap != 0)
		{
		myfree(tm_y_heap);
		tm_y_heap = 0;
		}
	tm_heap_cap = 0;
	}

bool hitdata::tm_coord_ensure(uint npairs)
	{
	if (npairs <= TM_BUFSIZE)
		return true;

	if (npairs <= tm_heap_cap)
		return true;

	const uint new_cap = npairs;
	double *nx = (double *) myalloc(double, new_cap*3);
	double *ny = (double *) myalloc(double, new_cap*3);
	tm_coord_free();
	tm_x_heap = nx;
	tm_y_heap = ny;
	tm_heap_cap = new_cap;
	return true;
	}

double *hitdata::tm_x_data(uint npairs)
	{
	if (npairs <= TM_BUFSIZE)
		return tm_x_buf;
	asserta(tm_x_heap != 0);
	return tm_x_heap;
	}

double *hitdata::tm_y_data(uint npairs)
	{
	if (npairs <= TM_BUFSIZE)
		return tm_y_buf;
	asserta(tm_y_heap != 0);
	return tm_y_heap;
	}

void hitdata::fill(const flat_params &params)
	{
	ids = 0;
	diffs = 0;
	gaps = 0;
	cigar_free();
	cigar_len = 0;

	const char *qaa = query->m_aa->m_data;
	const char *taa = target->m_aa->m_data;
	const uint LQ = query->m_L;
	const uint LT = target->m_L;
	uint qpos = qlo;
	uint tpos = tlo;

	if (ncol > 0)
		{
		char lastc = path[0];
		asserta(lastc == 'M');
		uint cigaropn = 0;
		if (qlo > 0)
			cigar_put_op(qlo, 'S');
		for (uint i = 0; i < ncol; ++i)
			{
			const char c = path[i];
			if (c == lastc)
				++cigaropn;
			else
				{
				cigar_put_op(cigaropn, lastc);
				cigaropn = 1;
				}

			switch (c)
				{
			case 'M':
				{
				asserta(qpos < LQ);
				asserta(tpos < LT);
				const char q = qaa[qpos];
				const char t = taa[tpos];
				if (toupper(q) == toupper(t))
					++ids;
				else
					++diffs;
				++qpos;
				++tpos;
				break;
				}

			case 'D':
				++gaps;
				++tpos;
				break;

			case 'I':
				++gaps;
				++qpos;
				break;

			default:
				asserta(false);
				}
			lastc = c;
			}
		cigar_put_op(cigaropn, lastc);
		}

	qhi = qpos - 1;
	thi = tpos - 1;
	pvalue = calc_pvalue(TS, params.m_pvm);
	TM = float(calc_tm_iterate());
	}
