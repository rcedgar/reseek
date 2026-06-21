#if 0
#include "myutils.h"
#include "seqdb.h"
#include "flat_chain.h"
#include "flat_params.h"
#include "lookup.h"
#include "alpha.h"
#include "entropy.h"

static uint get_unaligned_length(const string &seq)
	{
	uint n = 0;
	for (size_t i = 0; i < seq.size(); ++i)
		if (!isgap(seq[i]))
			++n;
	return n;
	}

void entropy::parse_fa2(SeqDB &DB, const bool is_tp)
	{
	DB.SetLabelToIndex();
	const uint nseq = DB.GetSeqCount();
	asserta(nseq%2 == 0);
	const uint npair = nseq/2;

	vector<string> flds;
	for (uint pairidx = 0; pairidx < npair; ++pairidx)
		{
		uint seqidxq = 2*pairidx;
		uint seqidxt = seqidxq + 1;

		string labelq = DB.GetLabel(seqidxq);
		const string &labelt = DB.GetLabel(seqidxt);

		Split(labelq, flds, ' ');
		labelq = flds[0];

		map<string, uint>::const_iterator iterq = m_label2idx.find(labelq);
		map<string, uint>::const_iterator itert = m_label2idx.find(labelt);
		asserta(iterq != m_label2idx.end());
		asserta(itert != m_label2idx.end());
		uint profidxq = iterq->second;
		uint profidxt = itert->second;

		const vector<vector<uint8_t> > &profq = m_profiles[profidxq];
		const vector<vector<uint8_t> > &proft = m_profiles[profidxt];
		const uint proflq = uint(profq[0].size());
		const uint proflt = uint(proft[0].size());

		const string &rowq = DB.GetSeq(seqidxq);
		const string &rowt = DB.GetSeq(seqidxt);
		const uint ncol = uint(rowq.size());
		asserta(rowt.size() == ncol);

		const uint ulq = get_unaligned_length(rowq);
		const uint ult = get_unaligned_length(rowt);

		const uint LQ = m_seqlengths[profidxq];
		const uint LT = m_seqlengths[profidxt];

		asserta(ulq == proflq);
		asserta(ult == proflt);

		asserta(ulq == LQ);
		asserta(ult == LT);

		vector<uint> posqs;
		vector<uint> posts;
		posqs.reserve(ncol);
		posts.reserve(ncol);
		uint posq = 0;
		uint post = 0;
		for (uint col = 0; col < ncol; ++col)
			{
			char q = rowq[col];
			char t = rowt[col];
			if (isupper(q) && isupper(t))
				{
				posqs.push_back(posq);
				posts.push_back(post);
				++m_total_col_count;
				}
			if (!isgap(q))
				++posq;
			if (!isgap(t))
				++post;
			}

		m_profidxqs.push_back(profidxq);
		m_profidxts.push_back(profidxt);
		m_posvecq.push_back(posqs);
		m_posvect.push_back(posts);
		m_pair_is_tp_vec.push_back(is_tp);
		}
	}

void entropy::load_fa2s(
	const string &tpfa2fn,
	const string &fpfa2fn)
	{
	m_TPDB.FromFasta(opt(fasta2_tp), true);
	m_FPDB.FromFasta(opt(fasta2_fp), true);

	parse_fa2(m_TPDB, true);
	parse_fa2(m_FPDB, false);
	}

float entropy::calc_col_score(
	const vector<vector<uint8_t> > &profileq, uint posq,
	const vector<vector<uint8_t> > &profilet, uint post) const
	{
	const size_t n = m_feature_subset.size();
	float total_score = 0;
	for (uint i = 0; i < n; ++i)
		{
		uint fi = m_feature_subset[i];
		assert(fi < profileq.size());
		assert(fi < profilet.size());
		assert(posq < profileq[fi].size());
		assert(post < profilet[fi].size());

		uint8_t codeq = profileq[fi][posq];
		uint8_t codet = profilet[fi][post];

		uint AS = m_alpha_sizes[fi];
		assert(codeq < AS);
		assert(codet < AS);
		float score = m_weighted_logoddsvec[fi][AS*codeq + codet];
		assert(score >= MIN_SANE_SCORE && score <= MAX_SANE_SCORE);
		total_score += score;
		}
	return total_score;
	}

void entropy::calc_col_scores()
	{
	m_scores.clear();
	m_is_tps.clear();
	m_scores.reserve(m_total_col_count);
	m_is_tps.reserve(m_total_col_count);

	const uint npair = uint(m_profidxqs.size());
	assert(m_profidxts.size() == npair);
	assert(m_posvecq.size() == npair);
	assert(m_posvect.size() == npair);
	assert(m_pair_is_tp_vec.size() == npair);

	for (uint pairidx = 0; pairidx < npair; ++pairidx)
		{
//		ProgressStep(pairidx, npair, "col scores");
		const uint profidxq = m_profidxqs[pairidx];
		const uint profidxt = m_profidxts[pairidx];
		assert(profidxq < m_profiles.size());
		assert(profidxt < m_profiles.size());

		const vector<vector<uint8_t> > &profq = m_profiles[profidxq];
		const vector<vector<uint8_t> > &proft = m_profiles[profidxt];

		const vector<uint> &posvecq = m_posvecq[pairidx];
		const vector<uint> &posvect = m_posvect[pairidx];
		const bool is_tp = m_pair_is_tp_vec[pairidx];
		const uint ncol = uint(posvecq.size());
		assert(uint(posvect.size()) == ncol);
		for (uint col = 0; col < ncol; ++col)
			{
			uint posq = posvecq[col];
			uint post = posvect[col];
			float score = calc_col_score(profq, posq, proft, post);
			m_scores.push_back(score);
			m_is_tps.push_back(is_tp);
			}
		}
	asserta(m_scores.size() == m_total_col_count);
	asserta(m_is_tps.size() == m_total_col_count);
	}

void entropy::set_logodds_subset(
	const vector<uint> &fis,
	const vector<float> &weights)
	{
	uint nfeat = uint(m_alpha_sizes.size());
	uint subset_size = uint(fis.size());
	asserta(weights.size() == subset_size);
	m_feature_subset = fis;
	for (uint i = 0; i < subset_size; ++i)
		asserta(m_feature_subset[i] < nfeat);

	asserta(m_feature_names.size() == nfeat);
	if (m_weights == 0)
		m_weights = myalloc(float, nfeat);
	memset(m_weights, 0, nfeat*sizeof(float));
	for (uint i = 0; i < subset_size; ++i)
		m_weights[m_feature_subset[i]] = weights[i];

	float sumw = 0;
	for (uint i = 0; i < subset_size; ++i)
		sumw += m_weights[m_feature_subset[i]];

	asserta(sumw > 1e-6);
	for (uint i = 0; i < subset_size; ++i)
		m_weights[m_feature_subset[i]] /= sumw;

	for (uint i = 0; i < subset_size; ++i)
		{
		uint fi = m_feature_subset[i];
		asserta(fi < nfeat);
		uint AS = m_alpha_sizes[fi];
		uint N = AS*AS;
		for (uint k = 0; k < N; ++k)
			{
			float uwscore = m_unweighted_logoddsvec[fi][k];
			assert(uwscore >= MIN_SANE_SCORE && uwscore <= MAX_SANE_SCORE);

			float wscore = uwscore*m_weights[fi];
			assert(wscore >= MIN_SANE_SCORE && wscore <= MAX_SANE_SCORE);

			m_weighted_logoddsvec[fi][k] = wscore;
			}
		}

	for (uint i = 0; i < subset_size; ++i)
		{
		uint fi = m_feature_subset[i];
		uint AS = m_alpha_sizes[fi];
		const float *low = m_weighted_logoddsvec[fi];
		const float *lou = m_unweighted_logoddsvec[fi];
		for (uint k = 0; k < AS*AS; ++k)
			{
			float wscore = low[k];
			float uscore = lou[k];
			asserta(wscore >= MIN_SANE_SCORE && wscore <= MAX_SANE_SCORE);
			asserta(uscore >= MIN_SANE_SCORE && uscore <= MAX_SANE_SCORE);
			}
		}
	}

float entropy::roc_auc(const vector<float>& scores,
	const vector<bool>& is_tp) const
	{
	const size_t N = scores.size();
	assert(is_tp.size() == N);
	assert(N > 0);

	size_t n_tp = 0;
	for (size_t i = 0; i < N; ++i)
		n_tp += (is_tp[i] ? 1 : 0);

	const size_t n_fp = N - n_tp;

	assert(n_tp > 0);
	assert(n_fp > 0);

	vector<size_t> order(N);
	for (size_t i = 0; i < N; ++i)
		order[i] = i;

	sort(order.begin(), order.end(),
		[&](size_t a, size_t b)
		{
		return scores[a] < scores[b];
		});

	float tp_rank_sum = 0.0;

	size_t i = 0;
	while (i < N)
		{
		size_t j = i + 1;
		const float s = scores[order[i]];
		while (j < N && scores[order[j]] == s)
			++j;

		// Tied block [i, j), 0-based positions in sorted order.
		// Average 1-based rank = ((i+1) + j) / 2.
		const float avg_rank = 0.5f * float(i + 1 + j);

		size_t tp_in_tie = 0;
		for (size_t k = i; k < j; ++k)
			tp_in_tie += (is_tp[order[k]] ? 1 : 0);

		tp_rank_sum += avg_rank * float(tp_in_tie);
		i = j;
		}

	const float npos = float(n_tp);
	const float nneg = float(n_fp);

	const float auc =
		(tp_rank_sum - npos * (npos + 1.0f) * 0.5f) / (npos * nneg);

	return auc;
	}

void cmd_train_fa2auc()
	{
	asserta(optset_lookup);
	asserta(optset_fasta2_tp);
	asserta(optset_fasta2_fp);

	entropy E;

	const string &filesfn = g_Arg1;
	vector<string> lines;
	ReadLinesFromFile(filesfn, lines);

	const size_t nfeat = lines.size();
	vector<string> fafns;
	vector<string> logoddsfns;
	vector<string> flds;
	vector<string> feature_names;
	for (size_t fi = 0; fi < nfeat; ++fi)
		{
		Split(lines[fi], flds, '\t');
		asserta(flds.size() == 2);

		const string &fafn = flds[0];
		const string &logoddsfn = flds[1];

		string fa_feature_name, logodds_feature_name;
		GetStemName(fafn, fa_feature_name);
		GetStemName(logoddsfn, logodds_feature_name);
		asserta(fa_feature_name == logodds_feature_name);
		feature_names.push_back(fa_feature_name);
		ProgressLog("[%3d] %s\n", fi, fa_feature_name.c_str());

		fafns.push_back(fafn);
		logoddsfns.push_back(logoddsfn);
		}

	E.load_profiles(fafns);
	E.read_logoddsvec(logoddsfns);
	E.load_fa2s(opt(fasta2_tp), opt(fasta2_fp));

	for (uint i = 0; i < E.m_feature_names.size(); ++i)
		{
		vector<uint> fis;
		vector<float> weights;
		const string &name = E.m_feature_names[i];

		fis.push_back(i);
		weights.push_back(1.0f);

		E.set_logodds_subset(fis, weights);
		E.calc_col_scores();

		float AUC = E.roc_auc(E.m_scores, E.m_is_tps);
		ProgressLog("name=%s AUC=%.4f\n", name.c_str(), AUC);
		}
	}

#endif
