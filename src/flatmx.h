#pragma once

// Matrix allocated as single buffer to enable
// single malloc with no initialization, this
// should be more efficient that STL vector
// of vectors which requires multiple mallocs
// and unnecessary initializtion.
template<class T>
class FlatMx
	{
public:
	size_t m_RowCount = 0;
	size_t m_ColCount = 0;
	size_t m_Size = 0;
	T *m_Buffer = 0;

private:
	FlatMx();	// disable default c'tor

public:
	FlatMx(size_t RowCount, size_t ColCount)
		{
		m_RowCount = RowCount;
		m_ColCount = ColCount;
		m_Size = RowCount*ColCount*sizeof(T);
		m_Buffer = (T *) malloc(m_Size);
		}

	~FlatMx()
		{
		assert(m_Buffer != 0);
		_chkmem();
		free(m_Buffer);
		_chkmem();
		m_Buffer = 0;
		}

	size_t Idx(uint i, uint j) const
		{
		assert(i < m_RowCount && j < m_ColCount);
		size_t Idx = i*m_ColCount + j;
		assert(Idx < m_Size);
		return Idx;
		}

	void Set(uint i, uint j, T x)
		{
		m_Buffer[Idx(i, j)] = x;
		}

	T Get(uint i, uint j) const
		{
		return m_Buffer[Idx(i, j)];
		}
	};
