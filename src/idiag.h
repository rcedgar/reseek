#pragma once

// Foldseek-compatible diagonals
// see $src/foldseek_rce/notes/diagonal_score.docx
class idiag
	{
public:
	int m_LQ = 0;
	int m_LT = 0;

	idiag(int LQ, int LT)
		{
		assert(LQ > 0 && LT > 0);
		m_LQ = LQ;
		m_LT = LT;
		}

	int getd(int i, int j) const
		{
		return i - j;
		}

	int getmind() const
		{
		return -(m_LT - 1);
		}

	int getmaxd() const
		{
		return m_LQ - 1;
		}

	int getj(int d, int i) const
		{
		int j = i - d;
		assert(j >= 0 && j < m_LT);
		return j;
		}

	int geti(int d, int j) const
		{
		int i = d + j;
		assert(i >= 0 && i < m_LQ);
		return i;
		}

	int getmini(int d) const
		{
		if (d < 0)
			return 0;
		assert(d < m_LQ);
		return d;
		}

	int getminj(int d) const
		{
		if (d > 0)
			return 0;
		assert(-d < m_LT);
		return -d;
		}

	int getmaxi(int d) const
		{
		int maxi = min(m_LQ - 1, d + m_LT-1);
		assert(maxi >= 0 && maxi < m_LQ);
		return maxi;
		}

	int getmaxj(int d) const
		{
		int maxj = min(m_LT - 1, m_LQ - 1 - d);
		assert(maxj >= 0 && maxj < m_LT);
		return maxj;
		}

	int getlen(int d) const
		{
		int maxi = getmaxi(d);
		int mini = getmini(d);
		int len = maxi - mini + 1;
		assert(len > 0);
		assert(len <= min(m_LQ, m_LT));
#if DEBUG
		int maxj = getmaxj(d);
		int minj = getminj(d);
		int len2 = maxj - minj + 1;
		assert(len2 == len);
#endif
		return len;
		}
	};
