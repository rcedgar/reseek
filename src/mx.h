#ifndef mx_h
#define mx_h

#include <mutex>

extern uint g_MxAllocCount;
extern uint g_MxFreeCount;
extern mutex g_MxLock;

template<class T> class Mx
	{
public:
	uint m_RowCount = 0;
	uint m_ColCount = 0;
	T **m_Data = 0;

	~Mx()
		{
		FreeData();
		}

	void Alloc(unsigned RowCount, unsigned ColCount)
		{
		FreeData();
		g_MxLock.lock();
		++g_MxAllocCount;
		g_MxLock.unlock();
		m_Data = myalloc(T *, RowCount);
		for (unsigned i = 0; i < RowCount; ++i)
			m_Data[i] = myalloc(T, ColCount);

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
			return;
			}
		g_MxLock.lock();
		++g_MxFreeCount;
		g_MxLock.unlock();
		for (unsigned i = 0; i < m_RowCount; ++i)
			myfree(m_Data[i]);
		myfree(m_Data);

		m_Data = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		}

	T **GetData()
		{
		return (T **) m_Data;
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
