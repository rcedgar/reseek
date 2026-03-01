#pragma once

#include "flat_enum.h"

extern atomic<int64_t> g_flat_creates[FE_N];
extern atomic<int64_t> g_flat_destroys[FE_N];
extern atomic<int64_t> g_flat_bytes[FE_N];

template<typename T, FE fe>
class flat_base
	{
public:
	atomic<int> m_refcount;
	T *m_data;
	uint32_t m_size;
	const char *m_srcfile;
	int m_srcline;

protected:
	flat_base(uint32_t n)
		{
		m_size = n;
		m_data = (T*) aligned_malloc(m_size*sizeof(T));
		g_flat_bytes[fe] += n*sizeof(T);
		++g_flat_creates[fe];
		m_refcount = 1;
		m_srcfile = 0;
		m_srcline = 0;
		}

	flat_base(uint32_t n, const char *srcfile, int srcline) :
		flat_base(n)
		{
		m_srcfile = srcfile;
		m_srcline = srcline;
		}

	~flat_base()
		{
		asserta(m_refcount == 0);
		++g_flat_destroys[fe];
		g_flat_bytes[fe] -= m_size*sizeof(T);
		aligned_free(m_data);
		m_data = 0;
		}

public:
	template<class Derived>
	void base_release(Derived *& p)
		{
		int r = --p->m_refcount;
		asserta(r >= 0);
		if (r == 0)
			delete p;
		p = 0;
		}
	};

class museq_t : public flat_base<uint8_t, FE_museq>
	{
public:
	museq_t(uint32_t n) : flat_base<uint8_t, FE_museq>(n)
		{ }
	museq_t(uint32_t n, const char *srcfile, int srcline) :
		flat_base<uint8_t, FE_museq>(n, srcfile, srcline)
		{ }
 	};

class chainaa_t : public flat_base<char, FE_chainaa>
	{
public:
	chainaa_t(uint32_t n) : flat_base<char, FE_chainaa>(n)
		{ }
	chainaa_t(uint32_t n, const char *srcfile, int srcline) :
		flat_base<char, FE_chainaa>(n, srcfile, srcline)
		{ }
 	};

template<typename T, FE fe>
class flat_mx : public flat_base<T, fe>
	{
public:
	uint32_t m_rows;
	uint32_t m_cols;

protected:
	flat_mx(int32_t rows, int32_t cols) :
		flat_base<T, fe>(rows*cols)
		{
		m_rows = rows;
		m_cols = cols;
		}

	flat_mx(int32_t rows, int32_t cols, const char *srcfile, int srcline) :
		flat_mx(rows, cols)
		{
		this->m_srcfile = srcfile;
		this->m_srcline = srcline;
		}

public:
	void set(uint32_t i, uint32_t j, T value)
		{
		assert(i < m_rows);
		assert(j < m_cols);
		this->m_data[i*m_cols + j] = value;
		}
	};

class chainxyz_t : public flat_mx<uint16_t, FE_chainxyz>
	{
public:
	chainxyz_t(uint32_t L) :
		flat_mx<uint16_t, FE_chainxyz>(L, 3)
		{ }
	chainxyz_t(uint32_t L, const char *srcfile, int srcline) :
		flat_mx<uint16_t, FE_chainxyz>(L, 3, srcfile, srcline)
		{ }
	};

#define create_museq(n)		new museq_t((n), __FILE__, __LINE__);
#define create_chainaa(n)	new chainaa_t((n), __FILE__, __LINE__);
#define create_chainxyz(n)	new chainxyz_t((n), __FILE__, __LINE__);

#define release(p)	p->base_release((p))
#define release0(p)	(p ? p->base_release((p)) : (void) 0)

void log_flat_stats(const string &msg = "");
