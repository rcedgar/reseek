#include "myutils.h"
#include "kmeans.h"

float KMeans::GetEuclideanDist(const vector<float> &v1, const vector<float> &v2)
	{
	const uint n = SIZE(v1);
	asserta(SIZE(v2) == n);
	float Sum2 = 0;
	for (uint i = 0; i < n; ++i)
		{
		float diff = v1[i] - v2[i];
		Sum2 += diff*diff;
		}
	return sqrt(Sum2);
	}

uint KMeans::GetBestFitClusterIdx(const vector<float> &v, float *ptrMinDist) const
	{
	asserta(SIZE(v) == m_M);
	float MinDist = FLT_MAX;
	uint BestCluster = UINT_MAX;
	for (uint k = 0; k < m_K; ++k)
		{
		float d = m_GetDist(v, m_Centroids[k]);
		if (k == 0 || d < MinDist)
			{
			BestCluster = k;
			MinDist = d;
			}
		}
	if (ptrMinDist != 0)
		*ptrMinDist = MinDist;
	return BestCluster;
	}

void KMeans::SetCentroidsGivenCurrentClusterAssignments()
	{
	m_Sizes.clear();
	m_Sums.clear();
	m_Sizes.resize(m_K);
	m_Sums.resize(m_K);
	for (uint k = 0; k < m_K; ++k)
		m_Sums[k].resize(m_M, 0);

	asserta(SIZE(m_ClusterIdxs) == m_N);
	const vector<vector<float> > &vs = *m_ptrvs;
	for (uint i = 0; i < m_N; ++i)
		{
		uint ClusterIdx = m_ClusterIdxs[i];
		asserta(ClusterIdx < m_K);
		m_Sizes[ClusterIdx] += 1;
		const vector<float> &v = vs[i];
		asserta(SIZE(v) == m_M);
		for (uint j = 0; j < m_M; ++j)
			{
			float x = v[j];
			m_Sums[ClusterIdx][j] += x;
			}
		}

	for (uint k = 0; k < m_K; ++k)
		{
		uint Count = m_Sizes[k];
		for (uint j = 0; j < m_M; ++j)
			{
			if (Count == 0)
				m_Centroids[k][j] = 0;
			else
				m_Centroids[k][j] = m_Sums[k][j]/Count;
			}
		}
	}

uint KMeans::AssignClustersGivenCurrentCentroids(float &TotalDist)
	{
	TotalDist = 0;
	m_Sizes.clear();
	m_Sizes.resize(m_K, 0);
	uint ChangeCount = 0;
	asserta(SIZE(m_ClusterIdxs) == m_N);
	const vector<vector<float> > &vs = *m_ptrvs;
	for (uint i = 0; i < m_N; ++i)
		{
		const vector<float> &v = vs[i];
		uint OldCluster = m_ClusterIdxs[i];
		float d = 0;
		uint NewCluster = GetBestFitClusterIdx(v, &d);
		TotalDist += d;
		m_Sizes[NewCluster] += 1;
		if (NewCluster != OldCluster)
			{
			m_ClusterIdxs[i] = NewCluster;
			++ChangeCount;
			}
		}
	return ChangeCount;
	}

void KMeans::Run(const vector<vector<float> > &vs, uint K, uint MaxIters)
	{
	m_N = SIZE(vs);
	asserta(m_N > 0);
	m_M = SIZE(vs[0]);
	for (uint i = 0; i < m_N; ++i)
		asserta(SIZE(vs[i]) == m_M);
	m_K = K;
	m_ptrvs = &vs;
	m_Centroids.clear();
	m_Centroids.resize(K, vector<float>(m_M, 0));
	vector<uint> ClusterIdxs;
	vector<uint> Sizes;
	m_ClusterIdxs.clear();
	m_ClusterIdxs.reserve(m_N);
	for (uint i = 0; i < m_N; ++i)
		m_ClusterIdxs.push_back(randu32()%K);

	m_MinChangeCount = UINT_MAX;
	m_MinTotalDist = FLT_MAX;
	float TotalDist0 = FLT_MAX;
	for (uint Iter = 0; Iter < MaxIters; ++Iter)
		{
		SetCentroidsGivenCurrentClusterAssignments();

		float TotalDist = 0;
		uint ChangeCount = AssignClustersGivenCurrentCentroids(TotalDist);
		if (Iter == 0)
			{
			TotalDist0 = TotalDist;
			if (TotalDist0 == 0)
				TotalDist0 = 1;
			}
		m_MinTotalDist = min(m_MinTotalDist, TotalDist);
		if(ChangeCount < m_MinChangeCount)
			{
			m_MinChangeCount = ChangeCount;
			m_BestCentroids = m_Centroids;
			}
		double Err = TotalDist/TotalDist0;
		if (ChangeCount == 0)
			{
			Progress("Iter %5u   changes %7u  err %.3g\n", Iter, ChangeCount, Err);
			break;
			}
		Progress("Iter %5u   changes %7u  err %.3g\r", Iter, ChangeCount, Err);
		Log("Iter %u, %u changes, err %.3g\n", Iter, ChangeCount, Err);
		}
	}

void KMeans::WriteCentroids(FILE *f) const
	{
	if (f == 0)
		return;
	string CmdLine;
	GetCmdLine(CmdLine);
	fprintf(f, "# %s\n", CmdLine.c_str());
	fprintf(f, "# [%s]\n", g_GitVer);
	fprintf(f, "K=%u", m_K);
	fprintf(f, "\tM=%u", m_M);
	fprintf(f, "\n");
	for(uint k = 0; k < m_K; ++k)
		{
		fprintf(f, "k=%u", k);
		fprintf(f, "\tsize=%u", m_Sizes[k]);
		for(uint j = 0; j < m_M; ++j)
			fprintf(f, "\t%.6g", m_Centroids[k][j]);
		fprintf(f, "\n");
		}
	}
