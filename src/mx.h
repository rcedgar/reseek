#ifndef mx_h
#define mx_h

#include <mutex>

#define ONE_BUFFER	1

template<class T> class Mx
	{
public:
	uint m_RowCount = 0;
	uint m_ColCount = 0;
	T **m_Data = 0;
#if ONE_BUFFER
	T *m_Buffer = 0;
#endif

private:
	Mx(const Mx &rhs);
	Mx(Mx &rhs);

public:
	Mx()
		{
		m_RowCount = 0;
		m_ColCount = 0;
		m_Data = 0;
		m_Buffer = 0;
		}

	~Mx()
		{
		FreeData();
		}

	void Alloc(unsigned RowCount, unsigned ColCount, const char *fn, int linenr)
		{
		asserta(RowCount > 0 && ColCount > 0);
		FreeData();
		m_Data = myalloc(T *, RowCount);
#if ONE_BUFFER
		uint n = RowCount*ColCount;
		m_Buffer = myalloc(T, n);
		for (uint i = 0; i < RowCount; ++i)
			m_Data[i] = m_Buffer + i*ColCount;
#else
		for (unsigned i = 0; i < RowCount; ++i)
			m_Data[i] = myalloc(T, ColCount);
#endif
		m_RowCount = RowCount;
		m_ColCount = ColCount;
		}

	void Clear()
		{
		FreeData();
		}

	void FreeData()
		{
		if (m_RowCount == 0)
			{
			asserta(m_ColCount == 0);
			asserta(m_Data == 0);
#if ONE_BUFFER
			asserta(m_Buffer == 0);
#endif
			return;
			}
#if ONE_BUFFER
		myfree(m_Buffer);
#else
		for (unsigned i = 0; i < m_RowCount; ++i)
			myfree(m_Data[i]);
#endif
		myfree(m_Data);

		m_Data = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		m_Buffer = 0;
		}

	T **GetData()
		{
		return m_Data;
		}

	T Get(unsigned i, unsigned j) const
		{
		assert(i < m_RowCount);
		assert(j < m_ColCount);
		return m_Data[i][j];
		}

	void Put(unsigned i, unsigned j, T x) const
		{
		assert(i < m_RowCount);
		assert(j < m_ColCount);
		m_Data[i][j] = x;
		}

	const T *const *const GetData() const
		{
		return (const T *const *) m_Data;
		}

	//void Copy(const Mx<T> &rhs)
	//	{
	//	Alloc("Copy", rhs.m_RowCount, rhs.m_ColCount, rhs.m_SeqDB, rhs.m_IdA, rhs.m_IdB);
	//	const T * const *Data = rhs.GetData();
	//	for (unsigned i = 0; i < m_RowCount; ++i)
	//		for (unsigned j = 0; j < m_ColCount; ++j)
	//			m_Data[i][j] = Data[i][j];
	//	}

	void Assign(T v)
		{
		for (unsigned i = 0; i < m_RowCount; ++i)
			for (unsigned j = 0; j < m_ColCount; ++j)
				m_Data[i][j] = v;
		}

	void Init(T v)
		{
		for (unsigned i = 0; i < m_RowCount; ++i)
			for (unsigned j = 0; j < m_ColCount; ++j)
				m_Data[i][j] = v;
		}
	};

void LogMxAllocCounts(const char *s);

#endif // mx_h
