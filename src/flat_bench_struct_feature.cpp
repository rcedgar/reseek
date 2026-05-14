#include "myutils.h"
#include "flat_bench_struct_feature.h"
#include "alpha.h"
#include "sort.h"
#include "triangle.h"
#include "flat_helpers.h"
#include "flat_aligner.h"

#define SHOW_PROGRESS 1

void flat_bench_struct_feature::ThreadBody_All(uint ThreadIdx)
	{
	asserta(!m_nu_filter);
	const uint NQ = SIZE(m_Labels);
	const uint PairCount = triangle_get_K(NQ);
	const uint nfeat = flat_alphas::get_nfeat();
	flat_aligner fa;
	fa.alloc();
	for (;;)
		{
		uint DomIdxT = m_NextQueryIdx++;
		if (DomIdxT >= NQ)
			return;

		const string &labelT = m_fp.get_label(DomIdxT);
		const uint8_t *profT = m_fp.get_profile(DomIdxT);
		const uint LT = m_fp.get_length(DomIdxT);
		fa.cacheT(labelT, profT, 0, LT);

		// Includes self-score for santify checking and because
		//   triangle*() functions include diagonal
		for (uint DomIdxQ = DomIdxT; DomIdxQ < NQ; ++DomIdxQ)
			{
			const string &labelQ = m_fp.get_label(DomIdxQ);
			const uint8_t *profQ = m_fp.get_profile(DomIdxQ);
			const uint LQ = m_fp.get_length(DomIdxQ);
			fa.alignQ(labelQ, profQ, 0, LQ);
			float Score = get_feature_value(DomIdxQ, DomIdxT, fa);
			if (isnan(Score))
				Die("isnan(%s,%s)", fa.m_labelQ.c_str(), fa.m_labelT.c_str());
			asserta(!isinf(Score));
			uint PairIdx = triangle_ij_to_k(DomIdxT, DomIdxQ, NQ);
			uint progress_count = m_aligned_pair_count++;
#if SHOW_PROGRESS
			if (ThreadIdx == 0 && progress_count%1000 == 0)
				ProgressStep(progress_count, PairCount, "Aligning");
#endif
			AppendHit(DomIdxT, DomIdxQ, Score);
			}
		}
	}

void flat_bench_struct_feature::ThreadBody_Dope(uint ThreadIdx)
	{
	assert(m_look);
	const uint ndom = m_look->get_ndom();
	const uint nfeat = flat_alphas::get_nfeat();
	flat_aligner fa;
	fa.alloc();
	uint CurrentDomIdxT = UINT_MAX;
	for (;;)
		{
		uint dopeidx = m_NextDopeIdx++;
		if (dopeidx >= m_dope_nhit)
			{
			fa.freemem();
			return;
			}
#if SHOW_PROGRESS
		if (ThreadIdx == 0 && dopeidx%1000 == 0)
			ProgressStep(dopeidx, m_dope_nhit, "Aligning");
#endif
		uint k = m_dope_ks[dopeidx];
		uint DomIdxQ, DomIdxT;
		triangle_k_to_ij(k, ndom, DomIdxT, DomIdxQ);

		if (DomIdxT == CurrentDomIdxT)
			++m_ncachehits;
		else
			{
			++m_ncachemisses;
			const string &labelT = m_fp.get_label(DomIdxT);
			const uint8_t *profT = m_fp.get_profile(DomIdxT);
			const uint LT = m_fp.get_length(DomIdxT);
			fa.cacheT(labelT, profT, 0, LT);
			CurrentDomIdxT = DomIdxT;
			}

		const string &labelQ = m_fp.get_label(DomIdxQ);
		const uint8_t *profQ = m_fp.get_profile(DomIdxQ);
		const uint LQ = m_fp.get_length(DomIdxQ);
		fa.alignQ(labelQ, profQ, 0, LQ);
		++m_aligned_pair_count;
		float Score = get_feature_value(DomIdxQ, DomIdxT, fa);
		asserta(!isnan(Score));
		asserta(!isinf(Score));
		uint PairIdx = triangle_ij_to_k(DomIdxT, DomIdxQ, ndom);
		AppendHit(DomIdxT, DomIdxQ, Score);
		}
	}

float flat_bench_struct_feature::get_feature_value(uint idxQ, uint idxT,
	const flat_aligner &fa) const
	{
	const string &feat = opt(feature);
	if (feat == "lddt")
		return get_lddt(idxQ, idxT, fa);
	else if (feat == "lddtpow")
		return get_lddtpow(idxQ, idxT, fa);
	else if (feat == "dali")
		return get_dali(idxQ, idxT, fa);
	else if (feat == "dalix")
		return get_dalix(idxQ, idxT, fa);
	else if (feat == "entropy")
		return get_entropy(idxQ, idxT, fa);
	Die("feat");
	return 0;
	}

float flat_bench_struct_feature::get_entropy(uint idxQ, uint idxT,
	const flat_aligner &fa) const
	{
	const uint nfeat = flat_alphas::get_nfeat();
	uint fi = UINT_MAX;
	for (uint i = 0; i < nfeat; ++i)
		{
		if (flat_alphas::m_alpha_names[i] == "sec32")
			{
			fi = i;
			break;
			}
		}
	if (fi == UINT_MAX)
		Die("entropy needs sec32");

	string path;
	fa.get_path_str(path);

	const uint8_t *profQ = m_fp.get_profile(idxQ);
	const uint8_t *profT = m_fp.get_profile(idxT);

	float H = flat_get_entropy(
		fa.m_labelQ, fa.m_labelT,
		path, fa.m_loQ, fa.m_LQ, fa.m_loT, fa.m_LT,
		profQ, profT, nfeat, fi);

	return H;
	}

float flat_bench_struct_feature::get_dali(uint idxQ, uint idxT,
	const flat_aligner &fa) const
	{
	asserta(idxQ < m_distmxs.size());
	asserta(idxT < m_distmxs.size());

	const sid_t *distmxQ = m_distmxs[idxQ];
	const sid_t *distmxT = m_distmxs[idxT];

	uint loQ = fa.m_loQ;
	uint loT = fa.m_loT;
	uint LQ = fa.m_LQ;
	uint LT = fa.m_LT;

	string path;
	uint nmatch = fa.get_path_str(path);
	uint ncol = uint(path.size());
	float dali = flat_get_dali(
		fa.m_labelQ, fa.m_labelT,
		path, loQ, LQ, loT, LT,
		distmxQ, distmxT);
	return dali;
	}

float flat_bench_struct_feature::get_dalix(uint idxQ, uint idxT,
	const flat_aligner &fa) const
	{
	asserta(idxQ < m_distmxs.size());
	asserta(idxT < m_distmxs.size());

	const sid_t *distmxQ = m_distmxs[idxQ];
	const sid_t *distmxT = m_distmxs[idxT];

	uint loQ = fa.m_loQ;
	uint loT = fa.m_loT;
	uint LQ = fa.m_LQ;
	uint LT = fa.m_LT;

	string path;
	uint nmatch = fa.get_path_str(path);
	uint ncol = uint(path.size());
	float *colscores = myalloc(float, nmatch);
	float dali = flat_get_dalix(
		fa.m_labelQ, fa.m_labelT,
		path, loQ, LQ, loT, LT,
		distmxQ, distmxT, colscores);
	myfree(colscores);
	return dali;
	}

float flat_bench_struct_feature::get_lddt(uint idxQ, uint idxT,
	const flat_aligner &fa) const
	{
	asserta(idxQ < m_distmxs.size());
	asserta(idxT < m_distmxs.size());

	const sid_t *distmxQ = m_distmxs[idxQ];
	const sid_t *distmxT = m_distmxs[idxT];

	uint loQ = fa.m_loQ;
	uint loT = fa.m_loT;
	uint LQ = fa.m_LQ;
	uint LT = fa.m_LT;

	float lddt = flat_getlddt_muscle_some_floats4(
		fa, distmxQ, distmxT);
	return lddt;
	}

float flat_bench_struct_feature::get_lddtpow(uint idxQ, uint idxT,
	const flat_aligner &fa) const
	{
	asserta(idxQ < m_distmxs.size());
	asserta(idxT < m_distmxs.size());

	const sid_t *distmxQ = m_distmxs[idxQ];
	const sid_t *distmxT = m_distmxs[idxT];

	uint loQ = fa.m_loQ;
	uint loT = fa.m_loT;
	uint LQ = fa.m_LQ;
	uint LT = fa.m_LT;

	string path;
	uint nmatch = fa.get_path_str(path);
	uint ncol = uint(path.size());
	float lddt = flat_getlddt_muscle_some_floats4(
		fa, distmxQ, distmxT);
	float maxL = max(LT, LQ) - 20.0f;
	if (maxL < 80)
		maxL = 80;
	float score = lddt*nmatch*2.0f/powf(maxL, 0.5);
	return score;
	}

void flat_bench_struct_feature::set_distmxs(uint M)
	{
	assert(m_look);
	const uint ndom = m_look->get_ndom();
	asserta(m_chains.size() == ndom);
	m_distmxs.clear();
	m_distmxs.reserve(ndom);
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		flat_chain_t *chain = m_chains[domidx];
		uint L = chain->get_length();
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain, distmx);
		m_distmxs.push_back(distmx);
		}
	}

void flat_bench_struct_feature::read_chains(const string &fn)
	{
	asserta(m_look);
	unordered_map<string, uint> label2idx;
	vector<flat_chain_t *> tmp_chains;
	read_flat_chains_idx_trunclabel(fn, tmp_chains, label2idx);
	m_chains.clear();
	uint ndom = m_look->get_ndom();
	m_chains.resize(ndom);
	for (unordered_map<string, uint>::const_iterator iter =
		m_look->m_dom2idx.begin();
		iter != m_look->m_dom2idx.end();
		++iter)
		{
		const string &label = iter->first;
		uint domidx = iter->second;
		unordered_map<string, uint>::const_iterator iter2 =
			label2idx.find(label);
		asserta(iter2 != label2idx.end());
		uint idx = iter2->second;
		flat_chain_t *chain = tmp_chains[idx];
		asserta(m_chains[domidx] == 0);
		m_chains[domidx] = chain;
		}
	}

void cmd_flat_bench_struct_feature()
	{
	asserta(optset_lookup);
	asserta(optset_fapattern);
	asserta(optset_mxpattern);
	asserta(optset_feature);
	asserta(optset_input);	// chains

	asserta(!optset_spec);
	asserta(!optset_varstr);

	const string &VarStr = g_Arg1;

	flat_bench_struct_feature FB;
	FB.ReadLookup(opt(lookup));
	FB.read_chains(opt(input));
	FB.set_distmxs(flat_params::m_distmx_bandwidth);
	if (optset_dope)
		FB.ReadDope(opt(dope));

	vector<string> param_names;
	vector<float> param_values;
	parse_varstr(VarStr, param_names, param_values);

	vector<string> feature_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		feature_names, weights,
		scalar_names, scalar_values);

	Die("TODO");
	//flat_alphas::load_alphas_obsolete(feature_names, opt(mxpattern));
	//FB.load_profiles_fapattern(opt(fapattern));
	FB.UpdateParamsFromVarStr(VarStr);
	FB.LogParams();
	FB.Alloc();

	uint ThreadCount = GetRequestedThreadCount();
	FB.Search(ThreadCount, false, opt(dope), UINT_MAX);
	FB.SetScoreOrder();
	FB.Bench();
	FB.WriteHits(opt(output), true);
	}
