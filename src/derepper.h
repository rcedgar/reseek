#pragma once

class Derepper
	{
public:
	const vector<vector<byte> > *m_vs = 0;
	vector<uint> m_IdxToCluster;
	vector<vector<uint> > m_ClusterToMemberIdxs;
	vector<vector<uint> > m_HashTable;
	uint m_N = UINT_MAX;
	uint m_K = UINT_MAX;
	uint m_HashSlots = UINT_MAX;

public:
	void Clear()
		{
		m_HashTable.clear();
		m_IdxToCluster.clear();
		m_ClusterToMemberIdxs.clear();
		m_N = UINT_MAX;
		m_HashSlots = UINT_MAX;
		}

	uint GetN() const { return SIZE(*m_vs); }
	void Run(const vector<vector<byte> > &vs, uint K);
	const vector<byte> &Get(uint Idx) const;
	uint GetHash(uint Idx) const;
	void Search(uint Idx);
	bool Eq(uint Idx1, uint Idx2) const;
	void Validate() const;
	void LogMe() const;
	void GetSizeOrder(vector<uint> &Order) const;
	};
