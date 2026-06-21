#pragma once

#include "myutils.h"

template<class T, uint MaxSize> class MaxBuff
	{
public:
	uint m_Size = 0;
	T *m_Data = 0;
	bool m_Filled = false;

public:
	~MaxBuff() { Free(); }

	void Free()
		{
		myfree(m_Data);
		m_Size = 0;
		m_Data = 0;
		m_Filled = false;
		}

	T *Alloc(uint Size)
		{
		m_Filled = false;
		if (Size <= m_Size && m_Size <= MaxSize)
			return m_Data;

		Free();
		m_Size = max(Size, MaxSize);
		m_Data = myalloc(T, m_Size);
		return m_Data;
		}

	void Filled()
		{
		m_Filled = true;
		}
	};
