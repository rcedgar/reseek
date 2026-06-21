#include "myutils.h"
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>

/**
 * Calculates a Latin Hypercube Sample.
 * @param M The number of samples to generate.
 * @param N The number of variables (dimensions).
 * @param mins vector size N containing minimum values for each variable.
 * @param maxs vector of size N containing maximum values for each variable.
 * @return result is a vector of vectors (size M x N) containing the sampled values.
 */
void latin_hypercube(
	const vector<double> &mins,
	const vector<double> &maxs,
	size_t M,
	vector<vector<double> > &result)
	{
	const size_t N = mins.size();
	asserta(maxs.size() == N);
	// result[sample_index][variable_index]
	result.resize(M, vector<double>(N));

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(0.0, 1.0);

	for (int j = 0; j < N; ++j) {
		// 1. Create a list of interval idxs [0, 1, ..., M-1]
		vector<int> idxs(M);
		iota(idxs.begin(), idxs.end(), 0);

		// 2. Shuffle idxs to ensure the "Latin" property (random permutation)
		shuffle(idxs.begin(), idxs.end(), gen);

		double range = maxs[j] - mins[j];
		double interval_width = range / M;

		for (int i = 0; i < M; ++i) {
			// 3. Pick a random point within the specific interval
			double lower_bound = mins[j] + (idxs[i] * interval_width);
			double random_offset = dis(gen) * interval_width;

			result[i][j] = lower_bound + random_offset;
			}
		}
	}