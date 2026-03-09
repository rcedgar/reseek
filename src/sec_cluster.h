#pragma once

#include "flat_distmx.h"
#include "chaq.h"

// Cluster subset of local distance
//	matrix by k-meeans clustering
class sec_cluster
	{
public:
	uint m_K = 0;				// number of clusters for K-means
	uint m_N = 0;				// number of residues, size of m_vs
	uint m_D = 0;				// dimension of feature vector, length of m_i/jvalues
	uint m_M = 0;				// band width for distance matrix (e.g. 64)
	int m_w = 0;				// band width for sec (e.g. 3), max index in m_i/jvalues
	const int* m_off1s = 0;		// +/- offsets from position
	const int* m_off2s = 0;		// +/- offsets from position
	uint* m_cluster_idxs = 0;	// current cluster assignments
	sid_t *m_means = 0;			// flat matrix of current means size m_K x m_D
	sid_t *m_vs;				// flat matrix of feature vectors size m_N x m_D
	uint* m_cluster_sizes = 0;	// cluster sizes
	const vector<flat_chain *> *m_chains = 0;

public:
	void log_params() const
		{
		Log("K %u, N %u, D %u, M %u, w %u\n", m_K, m_N, m_D, m_M, m_w);
		Log("off1s[%u] =", m_D);
		for (uint i = 0; i < m_D; ++i)
			{
			if (i > 0)
				Log(",");
			Log(" %2d", m_off1s[i]);
			}
		Log("\n");

		Log("off2s[%u] =", m_D);
		for (uint i = 0; i < m_D; ++i)
			{
			if (i > 0)
				Log(",");
			Log(" %2d", m_off2s[i]);
			}
		Log("\n");
		}

	void log_v(const sid_t *v) const
		{
		for (uint i = 0; i < m_D; ++i)
			Log(" %5u(%4.1f)", v[i], sid2dist(v[i]));
		Log("\n");
		}

	void log_means() const
		{
		Log("\nmeans:\n");
		for (uint cluster_idx = 0; cluster_idx < m_K; ++cluster_idx)
			{
			Log("%3u [%7u] ", cluster_idx, m_cluster_sizes[cluster_idx]);
			log_v(m_means + cluster_idx*m_D);
			}
		}

	void log_head_vs(uint n=10) const
		{
		Log("\nhead_vs(%u):\n", n);
		for (uint residue_idx = 0; residue_idx < min(n, m_N); ++residue_idx)
			{
			Log("[%3u] ", residue_idx);
			log_v(m_vs + residue_idx*m_D);
			}
		}

	void log_random_vs(uint n=100) const
		{
		Log("\nrandom_vs(%u):\n", n);
		for (uint i = 0; i < n; ++i)
			{
			uint residue_idx = randu32()%m_N;
			Log("[%7u] ", residue_idx);
			log_v(m_vs + residue_idx*m_D);
			}
		}

	void logme() const
		{
		log_params();
		log_random_vs();
		log_means();
		}

	// Adjacent residues in the backbone should have distance ~3.81 A
	// 5162 / 184690 bad backbones (3%)
	bool check_backbone(uint chain_idx, const sid_t* distmx, int pos, int L)
		{
		const sid_t backbone_sid = dist2sid(3.81f);
		const sid_t min_backbone_sid = backbone_sid - 10;
		const sid_t max_backbone_sid = backbone_sid + 10;
		assert(pos >= m_w && pos + m_w < L);
		for (int i = -int(m_w); i < int(m_w); ++i)
			{
			int ifirst_pos = int(pos)+i;
			assert(ifirst_pos >= 0 && ifirst_pos + 1 < L);
			uint32_t first_pos = uint32_t(ifirst_pos);
			uint k = banded_ij_to_k(m_M, first_pos, first_pos+1);
			sid_t sid = distmx[k];
			if (sid < min_backbone_sid || sid > max_backbone_sid)
				return false;
			}
		return true;
		}

	void get_v(uint chain_idx, const sid_t* distmx, int pos, int L, sid_t* v)
		{
		assert(pos >= m_w && pos + m_w < L);
		for (uint m = 0; m < m_D; ++m)
			{
			int off1 = m_off1s[m];
			int off2 = m_off2s[m];
			uint k = banded_ij_to_k(m_M, pos+off1, pos+off2);
			sid_t sid = distmx[k];
#if DEBUG
			assert(chain_idx < size(*m_chains));
			const flat_chain *chain = (*m_chains)[chain_idx];
			float dist_slow = chain->slow_float_dist(pos+off1, pos+off2);
			float dist_sid = sid2dist(sid);
			float diff = fabs(dist_slow - dist_sid);
			if (diff > 0.5f)
				Die("get_v(%s, pos1=%u, pos2=%u) sid=%u dist_sid=%.1f dist_slow=%.1f",
					chain->m_label.c_str(), pos+off1, pos+off2, sid, dist_sid, dist_slow);
#endif
			//_chkmem();//@@
			v[m] = sid;
			//_chkmem();//@@//FAILS HERE
			}
		}

	// Euclidean squared distance (no need to sqrt)
	uint32_t get_dist(const sid_t* v1, const sid_t* v2) const
		{
		uint32_t sum2 = 0;
		for (uint m = 0; m < m_D; ++m)
			{
			uint32_t diff = int32_t(v1[m]) - int32_t(v2[m]);
			sum2 += diff*diff;
			}
		return sum2;
		}

	// When looking for best match can give up early
	uint32_t get_dist_early_quit(const sid_t* v1, const sid_t* v2,
		uint32_t smallest_so_far) const
		{
		uint32_t sum2 = 0;
		for (uint m = 0; m < m_D; ++m)
			{
			uint32_t diff = int32_t(v1[m]) - int32_t(v2[m]);
			sum2 += diff*diff;
			if (sum2 >= smallest_so_far)
				return sum2;
			}
		return sum2;
		}
	
	uint assign_cluster(const sid_t* v) const
		{
		uint best_cluster = UINT32_MAX;
		uint32_t min_dist = UINT32_MAX;
		for (uint cluster_idx = 0; cluster_idx < m_K; ++cluster_idx)
			{
			uint32_t d = get_dist_early_quit(v, m_means + cluster_idx*m_D, min_dist);
			if (d < min_dist)
				{
				min_dist = d;
				best_cluster = cluster_idx;
				}
			}
		assert(best_cluster != UINT_MAX);
		return best_cluster;
		}

	void assign_random_means()
		{
		for (uint cluster_idx = 0; cluster_idx < m_K; ++cluster_idx)
			{
			uint residue_idx = randu32()%m_N;
			memcpy(m_means + cluster_idx*m_D, m_vs + residue_idx*m_D, m_D*sizeof(sid_t));
			//@@TODO
			{
			Log("\n");
			Log("Random mean %u:\n", cluster_idx);
			log_v(m_means + cluster_idx*m_D);
			}
			}
		}

	uint assign_clusters()
		{
		if (m_cluster_sizes != 0)
			myfree(m_cluster_sizes);
		m_cluster_sizes = myalloc(uint, m_K);
		zero_array(m_cluster_sizes, m_K);
		uint nrchanges = 0;
		for (uint residue_idx = 0; residue_idx < m_N; ++residue_idx)
			{
			uint old_cluster_idx = m_cluster_idxs[residue_idx];
			uint new_cluster_idx = assign_cluster(m_vs + residue_idx*m_D);
			if (new_cluster_idx != old_cluster_idx)
				{
				++nrchanges;
				m_cluster_idxs[residue_idx] = new_cluster_idx;
				}
			++m_cluster_sizes[new_cluster_idx];
			}
		return nrchanges;
		}

	sid_t get_random_value(uint d) const
		{
		assert(d < m_D);
		uint residue_idx = randu32()%m_N;
		return m_vs[residue_idx*m_D + d];
		}

	uint calc_means()
		{
		uint zero_count = 0;
		uint n = m_K*m_D;
		uint64_t *sums = myalloc(uint64_t, n);
		uint32_t *residue_counts = myalloc(uint32_t, m_K);
		zero_array(sums, n);
		zero_array(residue_counts, m_K);
#if DEBUG
		uint32_t *check_counts = myalloc(uint32_t, n);
		zero_array(check_counts, n);
#endif

		for (uint residue_idx = 0; residue_idx < m_N; ++residue_idx)
			{
			uint cluster_idx = m_cluster_idxs[residue_idx];
			++residue_counts[cluster_idx];
			for (uint d = 0; d < m_D; ++d)
				{
				sums[cluster_idx*m_D + d] += m_vs[residue_idx*m_D + d];
#if DEBUG
				check_counts[cluster_idx*m_D + d] += 1;
#endif
				}
			}

		uint sum_residue_count = 0;
		for (uint cluster_idx = 0; cluster_idx < m_K; ++cluster_idx)
			{
			uint residue_count = residue_counts[cluster_idx];
			sum_residue_count += residue_count;
			if (residue_count == 0)
				{
				++zero_count;
				uint random_residue_idx = randu32()%m_N;
				for (uint d = 0; d < m_D; ++d)
					m_means[cluster_idx*m_D + d] = m_vs[random_residue_idx*m_D + d];
				continue;
				}

			for (uint d = 0; d < m_D; ++d)
				{
#if DEBUG
				uint32_t check_count = check_counts[cluster_idx*m_D + d];
				assert(check_count == residue_count);
#endif
				uint64_t mean64 = sums[cluster_idx*m_D + d]/residue_count;
				sid_t mean = sid_t(mean64);
				asserta(uint64_t(mean) == mean64);
				m_means[cluster_idx*m_D + d] = mean;
				}
			}
		assert(sum_residue_count == m_N);
		myfree(sums);
		myfree(residue_counts);
		return zero_count;
		}

	void init(
		const vector<int> &off1s,
		const vector<int> &off2s)
		{
		m_D = SIZE(off1s);
		asserta(SIZE(off2s) == m_D);
		m_off1s = off1s.data();
		m_off2s = off2s.data();
		m_w = 0;
		for (uint i = 0; i < m_D; ++i)
			{
			m_w = max(m_w, abs(off1s[i]));
			m_w = max(m_w, abs(off2s[i]));
			}
		}

	void set_vs(const vector<flat_chain *> &chains)
		{
		m_N = 0;
		m_chains = &chains;
		asserta(m_D > 0);
		const uint nrchains = SIZE(chains);

		uint total_length = 0;
		for (uint chainidx = 0; chainidx < nrchains; ++chainidx)
			{
			uint L = chains[chainidx]->get_length();
			total_length += L;
			}

		// Will skip some residues, total_length is > m_N
		m_vs = myalloc(sid_t, m_D*total_length);
#if DEBUG
		memset(m_vs, 0xff, m_D*total_length*sizeof(sid_t));
#endif

		chaq c;
		uint residue_idx = 0;
		uint bad_backbones = 0;
		for (uint chain_idx = 0; chain_idx < nrchains; ++chain_idx)
			{
			const flat_chain *chain = chains[chain_idx];
			const int L = (int) chain->get_length();

			c.init(chain);
			const chaindistmx_t* dm = c.get_distmx(m_M);
			const sid_t *distmx = dm->m_data;
			for (int pos = m_w; pos < L - m_w; ++pos)
				{
				bool ok = check_backbone(chain_idx, distmx, pos, L);
				if (!ok)
					{
					++bad_backbones;
					continue;
					}
				get_v(chain_idx, distmx, pos, L, m_vs + m_D*residue_idx++);
				}
			}
		m_N = residue_idx;
		ProgressLog("%u / %u bad backbones\n", bad_backbones, m_N);
		}
	};
