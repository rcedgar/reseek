#pragma once

typedef float KMEANS_GET_DIST_FUNC(
	const vector<float> &v1, const vector<float> &v2);

class KMeans
	{
public:
	uint m_K = UINT_MAX;	// Number of clusters
	uint m_M = UINT_MAX;	// Number of dimensions
	uint m_N = UINT_MAX;	// Number of input vectors
	const vector<vector<float> > *m_ptrvs = 0;
	KMEANS_GET_DIST_FUNC *m_GetDist = GetEuclideanDist;
	uint m_MinChangeCount = UINT_MAX;

// Vectors of size N
	vector<uint> m_ClusterIdxs;

// Vectors of size K
	vector<uint> m_Sizes;

// Arrays of size K x M
	vector<vector<float> > m_Centroids;
	vector<vector<float> > m_BestCentroids;
	vector<vector<float> > m_Sums;

public:
	void Run(const vector<vector<float> > &vs, uint K, uint MaxIters);
	void SetCentroidsGivenCurrentClusterAssignments();
	uint GetBestFitClusterIdx(const vector<float> &v) const;
	uint AssignClustersGivenCurrentCentroids();

public:
	static float GetEuclideanDist(
		const vector<float> &v1, const vector<float> &v2);
	};
