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
	BPData *m_BPs = 0;
	float *m_ChainScores = 0;
	uint *m_TB = 0;

public:
	void Chain(const uint *Los, const uint *His,
	  const float *Scores, uint N, vector<uint> &Idxs);
	void Clear();

	static float GetChainScore(const uint *Los, const uint *His, const float *Scores,
	  uint N, const vector<uint> &Idxs);
	};