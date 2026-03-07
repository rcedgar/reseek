#pragma once

/***
Action			Code						Comments
------			----						--------
Create object	ptr_museq = create_museq(n)	initial refcount=1
Refcopy object	up(ptr_museq)				++refcount
Release object	down(ptr_museq)				deletes when refcount=1
				down0(ptr_museq)			allows nullptr (e.g. d'tor)

Function returning an object pointer p increments refcount.
Caller ensures that down(p) is called.

Simple lifetime
---------------
Effectively new ... delete inside function body.
	ptr_museq = create_museq(n);
	// ...
	down(ptr_museq);

Returning pointer to object from create_xxx()
---------------------------------------------
Same idiom as create_xxx().
	ptr_museq = my_chaq->get_museq();
	// ... same function ...
	down(ptr_museq);	// nullptr not allowed

Returning object to higher caller
---------------------------------
Reference count is not changed, now it its higher
caller's responsibility to down().
	ptr_museq = ...;	// create_xxx() or call to lower object
	// ... arbitrary code ...
	return ptr_museq

Member pointer
--------------
	m_ptr_museq = nullptr;	// c'tor
	// ... arbitary code ...
	m_ptr_museq = create_xxx(); OR lower_obj->get_museq(); // anywhere
	// ... arbitary code ...
	down0(m_ptr_museq);		// d'tor, allows nullptr if never set

Returning member pointer
------------------------
Increment refcount, caller will decrement when no longer needed.
	up(m_ptr_museq);	// increment refcount immediately before return
	return m_ptr_museq;
***/

#include "flat_enum.h"
#include "flat_dist_types.h"

extern const ic_t sid2ic[65536];

// Global atomics shared by all threads
// Simpler, faster and smaller compared to ObjMgr
// where one object per thread.
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

template<typename T, FE fe>
class flat_vec : public flat_base<T, fe>
	{
public:
	flat_vec(uint32_t n) : flat_base<T, fe>(n)
		{ }
	flat_vec(uint32_t n, const char *srcfile, int srcline) :
		flat_base<T, fe>(n, srcfile, srcline)
		{ }
	};

template<typename T, FE fe>
class flat_mx : public flat_base<T, fe>
	{
public:
	uint32_t m_rows;
	uint32_t m_cols;

public:
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
	T get(uint32_t i, uint32_t j) const
		{
		assert(i < m_rows);
		assert(j < m_cols);
		return this->m_data[i*m_cols + j];
		}
	};

class chainxyz_t : public flat_mx<uint16_t, FE_chainxyz>
	{
public:
	chainxyz_t(uint32_t L) :
		flat_mx<uint16_t, FE_chainxyz>(3, L)
		{ }
	chainxyz_t(uint32_t L, const char *srcfile, int srcline) :
		flat_mx<uint16_t, FE_chainxyz>(3, L, srcfile, srcline)
		{ }
	};

using chainaa_t = flat_vec<char, FE_chainaa>;
using museq_t = flat_vec<uint8_t, FE_museq>;
using featseq_t = flat_vec<uint8_t, FE_museq>;
using ss3_t = flat_vec<char, FE_ss3>;
using nnvec_t = flat_vec<uint16_t, FE_nnvec>;
using floatvec_t = flat_vec<float, FE_floatvec>;

using chaindistmx_t = flat_vec<sid_t, FE_chaindistmx>;
using megaprof_t = flat_mx<uint8_t, FE_chaindistmx>;

#define create_floatvec(n)	new floatvec_t((n), __FILE__, __LINE__);
#define create_museq(n)		new museq_t((n), __FILE__, __LINE__);
#define create_featseq(n)	new featseq_t((n), __FILE__, __LINE__);
#define create_nnvec(n)		new nnvec_t((n), __FILE__, __LINE__);
#define create_ss3(n)		new ss3_t((n), __FILE__, __LINE__);
#define create_chainaa(n)	new chainaa_t((n), __FILE__, __LINE__);
#define create_chainxyz(n)	new chainxyz_t((n), __FILE__, __LINE__);

#define create_chaindistmx(L, M)	new chaindistmx_t((L)*(M), __FILE__, __LINE__);
#define create_megaprof(nfeat, L)	new megaprof_t((nfeat), (L), __FILE__, __LINE__);

#define down(p)		p->base_release((p))
#define down0(p)	(p ? p->base_release((p)) : (void) 0)

void log_flat_stats(const string &msg = "");
