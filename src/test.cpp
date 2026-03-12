#include "myutils.h"
#include "flat_enum.h"

using sid_t = int;
using ic_t = int;

template<typename T, FE fe>
class flat_base
	{
private:
	atomic<int> m_refcount;

public:
	T* m_data;
	uint32_t m_size;
#if TRACK_SRC
	FE m_fe = fe;
	const char *m_srcfile = 0;
	int m_srcline = 0;
	list<void *>::iterator m_list_iter;
#endif

protected:
	flat_base()
		{
		m_size = 0;
		m_data = (T*) aligned_malloc(m_size*sizeof(T));
		m_refcount = 1;
#if TRACK_SRC
		m_srcfile = 0;
		m_srcline = 0;
#endif
		}

	~flat_base()
		{
		asserta(m_refcount == 0);
		aligned_free(m_data);
		m_data = 0;
#if TRACK_SRC
		if (m_srcfile)
			{
			g_flat_obj_list_lock.lock();
			g_flat_obj_list.erase(m_list_iter);
			g_flat_obj_list_lock.unlock();
			}
#endif
		}

public:
	void add_ref()
		{
		++m_refcount;
		}

	void release_ref()
		{
		assert(m_refcount > 0);
		--m_refcount;
		}

	int get_refcount() const
		{
		return m_refcount;
		}

	void falloc(uint32_t n)
		{
		asserta(m_size == 0);
		m_size = n;
		m_data = (T*) aligned_malloc(m_size*sizeof(T));
		g_flat_bytes[fe] += n*sizeof(T);
#if TRACE_SRC
		m_srcfile = srcfile;
		m_srcline = srcline;
#endif
		}
	};

template<typename T, FE fe>
class flat_vec : public flat_base<T, fe>
	{
protected:
	flat_vec() : flat_base<T, fe>() {} // explicit do-nothing
	};

template<typename T, FE fe>
class flat_mx : public flat_base<T, fe>
	{
public:
	uint32_t m_rows = 0;
	uint32_t m_cols = 0;

protected:
	flat_mx() : flat_base()
		{
		m_rows = 0;
		m_cols = 0;
		}

	flat_mx(uint32_t rows, uint32_t cols) : flat_base()
		{
		falloc2(rows, cols)
		}
	};

class chaindistmx_t : public flat_mx<ic_t, FE_chaindistmx>
	{
protected:
	static chaindistmx_t *new_chaindistmx()
		{
		chaindistmx_t *p = new chaindistmx_t;
		return p;
		}
	};
