#pragma once

class diag
	{
public:
	int m_LQ = 0;
	int m_LT = 0;

	diag(int LQ, int LT)
		{
		assert(LQ > 0 && LT > 0);
		m_LQ = LQ;
		m_LT = LT;
		}

	int getd(int i, int j) const
		{
		return (m_LQ + j) - i - 1;
		}

	int getmind() const
		{
		return 0;
		}

	int getmaxd() const
		{
		return (m_LQ + m_LT) - 2;
		}

	int getj(int d, int i) const
		{
		int j = (d + i + 1) - m_LQ ;
		assert(j >= 0 && j < m_LT);
		return j;
		}

	int geti(int d, int j) const
		{
		int i = (m_LQ + j) - 1 - d;
		assert(i >= 0 && i < m_LQ);
		return i;
		}

	int getmini(int d) const
		{
		int mini = m_LQ - d - 1;
		if (mini < 0)
			mini = 0;
		assert(mini < m_LQ);
		return mini;
		}

	int getminj(int d) const
		{
		int minj = d + 1 - m_LQ;
		if (minj < 0)
			minj = 0;
		assert(minj < m_LT);
		return minj;
		}

	int getmaxi(int d) const
		{
		int maxi = m_LQ + m_LT - d - 2;
		if (maxi >= m_LQ)
			maxi = m_LQ - 1;
		return maxi;
		}

	int getmaxj(int d) const
		{
		int maxj = d;
		if (maxj >= m_LT)
			maxj = m_LT - 1;
		return maxj;
		}
	};
