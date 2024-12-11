struct BPData
	{
	uint Pos;
	bool IsLo;
	uint Index;

	void LogMe() const
		{
		Log("BP%s Pos %u Ix %u", (IsLo ? "lo" : "hi"), Pos, Index);
		}
	};

class Chainer
	{
public:
	Chainer()
		{
		m_BPs = 0;
		m_ChainScores = 0;
		m_TB = 0;
		}
	~Chainer()
		{
		Clear();
		}

public:
	BPData *m_BPs = 0;
	float *m_ChainScores = 0;
	uint *m_TB = 0;

public:
	float Chain(const vector<uint> &Los, const vector<uint> &His,
	  const vector<float> &Scores, vector<uint> &Idxs);
	void Clear()
		{
		if (m_BPs != 0) { myfree(m_BPs); m_BPs = 0; }
		if (m_ChainScores != 0) { myfree(m_ChainScores); m_ChainScores = 0; }
		if (m_TB != 0) { myfree(m_TB); }
		}


	static float GetChainScore(const uint *Los, const uint *His, const float *Scores,
	  uint N, const vector<uint> &Idxs);
	};