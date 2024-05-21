#include "myutils.h"
#include "derepper.h"
#include "murmur.h"
#include "sort.h"

const vector<byte> &Derepper::Get(uint Idx) const
	{
	const vector<vector<byte> > &vs = *m_vs;
	asserta(Idx < SIZE(vs));
	return vs[Idx];
	}

uint Derepper::GetHash(uint Idx) const
	{
	const vector<byte> &v = Get(Idx);
	uint mur = murmur3_32(v.data(), SIZE(v));
	uint h = mur%m_HashSlots;
	return h;
	}

bool Derepper::Eq(uint Idx1, uint Idx2) const
	{
	const vector<byte> &v1 = Get(Idx1);
	const vector<byte> &v2 = Get(Idx2);
	asserta(SIZE(v1) == m_K && SIZE(v2) == m_K);
	for (uint i = 0; i < m_K; ++i)
		{
		if (v1[i] != v2[i])
			return false;
		}
	return true;
	}

void Derepper::Search(uint Idx)
	{
	uint h = GetHash(Idx);
	const vector<byte> &v = Get(Idx);
	for (uint i = 0; i < m_HashSlots; ++i)
		{
		vector<uint> &Row = m_HashTable[h];
		if (Row.empty())
			{
			uint Cluster = SIZE(m_ClusterToMemberIdxs);
			Row.push_back(Idx);
			m_IdxToCluster.push_back(Cluster);

			vector<uint> Members;
			Members.push_back(Idx);
			m_ClusterToMemberIdxs.push_back(Members);
			return;
			}
		for (uint j = 0; j < SIZE(Row); ++j)
			{
			uint Idx2 = Row[j];
			if (Eq(Idx, Idx2))
				{
				asserta(Idx2 < SIZE(m_IdxToCluster));
				uint Cluster = m_IdxToCluster[Idx2];
				m_IdxToCluster.push_back(Cluster);
				m_ClusterToMemberIdxs[Cluster].push_back(Idx);
				return;
				}
			}
		h = (h + 1)%m_HashSlots;
		}
	}

void Derepper::Validate() const
	{
	const uint ClusterCount = SIZE(m_ClusterToMemberIdxs);

	for (uint Idx = 0; Idx < m_N; ++Idx)
		{
		uint Cluster = m_IdxToCluster[Idx];
		asserta(Cluster < ClusterCount);
		const vector<uint> &MemberIdxs = m_ClusterToMemberIdxs[Cluster];
		bool Found = false;
		for (uint i = 0; i < SIZE(MemberIdxs); ++i)
			{
			if (MemberIdxs[i] == Idx)
				{
				Found = true;
				break;
				}
			}
		asserta(Found);
		}

	vector<bool> Done(m_N);
	for (uint Cluster = 0; Cluster < ClusterCount; ++Cluster)
		{
		const vector<uint> &MemberIdxs = m_ClusterToMemberIdxs[Cluster];
		asserta(!MemberIdxs.empty());
		for (uint i = 0; i < SIZE(MemberIdxs); ++i)
			{
			uint Idx = MemberIdxs[i];
			asserta(!Done[Idx]);
			Done[Idx] = true;
			}

		for (uint i = 0; i < SIZE(MemberIdxs); ++i)
			{
			uint Idxi = MemberIdxs[i];
			for (uint j = i; j < SIZE(MemberIdxs); ++j)
				{
				uint Idxj = MemberIdxs[j];
				asserta(Eq(Idxi, Idxj));
				}
			}
		}

	for (uint Idx = 0; Idx < m_N; ++Idx)
		asserta(Done[Idx]);
	}

void Derepper::Run(const vector<vector<byte> > &vs, uint K)
	{
	Clear();
	m_vs = &vs;
	m_K = K;
	m_N = GetN();
	m_HashSlots = 3*m_N  + 1;
	m_HashTable.resize(m_HashSlots);

	for (uint Idx = 0; Idx < m_N; ++Idx)
		Search(Idx);
	}

void Derepper::LogMe() const
	{
	const vector<vector<byte> > &vs = *m_vs;
	const uint ClusterCount = SIZE(m_ClusterToMemberIdxs);
	Log("%u clusters\n", ClusterCount);
	for (uint Cluster = 0; Cluster < ClusterCount; ++Cluster)
		{
		const vector<uint> &MemberIdxs = m_ClusterToMemberIdxs[Cluster];
		Log("%7u  [%7u]  ", Cluster, SIZE(MemberIdxs));
		uint Idx = MemberIdxs[0];
		for (uint i = 0; i < m_K; ++i)
			Log(" %3u", vs[Idx][i]);
		Log("\n");
		}
	}

void Derepper::GetSizeOrder(vector<uint> &Order) const
	{
	Order.clear();
	const uint ClusterCount = SIZE(m_ClusterToMemberIdxs);
	vector<uint> Sizes;
	for (uint i = 0; i < ClusterCount; ++i)
		{
		uint Size = SIZE(m_ClusterToMemberIdxs[i]);
		Sizes.push_back(Size);
		}
	Order.resize(ClusterCount);
	QuickSortOrderDesc(Sizes.data(), ClusterCount, Order.data());
	}

void cmd_test_derep()
	{
	vector<vector<byte> > vs;

	for (uint i = 0; i < 10000; ++i)
		{
		vector<byte> v;
		v.push_back(randu32()%4);
		v.push_back(randu32()%4);
		v.push_back(randu32()%4);
		v.push_back(randu32()%4);
		vs.push_back(v);
		}

	Derepper D;
	D.Run(vs, 4);
	D.LogMe();
	D.Validate();
	vector<uint> Order;
	D.GetSizeOrder(Order);
	for (uint i = 0; i < min(SIZE(Order), 32u); ++i)
		{
		uint Cluster = Order[i];
		uint Size = SIZE(D.m_ClusterToMemberIdxs[Cluster]);
		}
	}
