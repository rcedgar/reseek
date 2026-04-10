#include "myutils.h"
#include "landmarks.h"
#include "flat_chain.h"
#include "chaq.h"
#include "seqdb.h"

static const uint M = 32;
static const uint MINL = 80;
static const uint MAXL = 1000;

static uint32_t n_ok_flanks = 0, n_ok_span = 0, n_ok_curv = 0, n_emit = 0;

static void find_landmark_candidates(
	const vector<xyz_t> &ca,
	const vector<sstype_t> &ss,
	const landmark_params_t &opt,
	vector<landmark_candidate_t> &out)
	{
	out.clear();
	const uint32_t n = (uint32_t) ca.size();
	if (n != (uint32_t) ss.size())
		return;
	if (n < 7)
		return;
	out.reserve(n);

	// Edge validity between i and i+1.
	vector<uint8_t> ok_edge(n > 0 ? n - 1 : 0, 1);
	if (opt.suppress_across_chain_breaks)
		{
		for (uint32_t i = 0; i + 1 < n; ++i)
			{
			double d = lm_dist(ca[i], ca[i+1]);
			ok_edge[i] = (d <= opt.max_adjacent_ca_dist ? 1 : 0);
			}
		}

	// Local geometry.
	vector<double> curv_deg(n, 0.0);        // angle(i-1,i,i+1)
	vector<double> tors_deg(n, 0.0);        // dihedral(i-1,i,i+1,i+2), stored at i
	vector<double> tors_delta_deg(n, 0.0);

	for (uint32_t i = 1; i + 1 < n; ++i)
		{
		if (opt.suppress_across_chain_breaks)
			{
			if (!ok_edge[i-1] || !ok_edge[i])
				continue;
			}
		curv_deg[i] = lm_angle_deg(ca[i-1], ca[i], ca[i+1]);
		}

	for (uint32_t i = 1; i + 2 < n; ++i)
		{
		if (opt.suppress_across_chain_breaks)
			{
			if (!ok_edge[i-1] || !ok_edge[i] || !ok_edge[i+1])
				continue;
			}
		tors_deg[i] = lm_dihedral_deg(ca[i-1], ca[i], ca[i+1], ca[i+2]);
		}

	for (uint32_t i = 2; i + 1 < n; ++i)
		tors_delta_deg[i] = lm_absdiff(tors_deg[i], tors_deg[i-1]);

	// Run lengths of SS states.
	vector<uint32_t> run_left(n, 1), run_right(n, 1);
	for (uint32_t i = 0; i < n; ++i)
		{
		run_left[i] = lm_same_state_run_left(ss, i);
		run_right[i] = lm_same_state_run_right(ss, i);
		}

	// Score tracks.
	vector<double> span_score(n, 0.0);
	vector<double> curvature_score(n, 0.0);
	vector<double> torsion_score(n, 0.0);
	vector<double> strand_hairpin_score(n, 0.0);
	vector<double> helix_kink_score(n, 0.0);
	vector<double> tr_hc_score(n, 0.0);
	vector<double> tr_ch_score(n, 0.0);
	vector<double> tr_ec_score(n, 0.0);
	vector<double> tr_ce_score(n, 0.0);
	vector<double> compact_score(n, 0.0);
	vector<double> nonlocal_score(n, 0.0);

	vector<double> best_span_dist(n, 1e30);
	vector<uint32_t> best_span_k(n, 0);

	double span_feature_threshold = opt.span_close_threshold;
	if (opt.enable_strand_hairpin_turn &&
		opt.strand_hairpin_span_threshold > span_feature_threshold)
		span_feature_threshold = opt.strand_hairpin_span_threshold;

	for (uint32_t i = 0; i < n; ++i)
		{
		// Precompute best span-closure feature for any detector that needs it.
		if (opt.enable_span_closure || opt.enable_strand_hairpin_turn)
			{
			double best_score = -1e30;
			double best_d = 1e30;
			uint32_t best_k = 0;

			for (uint32_t k = opt.span_k_min; k <= opt.span_k_max; ++k)
				{
				if (i < k || i + k >= n)
					continue;

				if (opt.suppress_across_chain_breaks)
					{
					if (lm_has_chain_break(ok_edge, i-k, i) ||
						lm_has_chain_break(ok_edge, i, i+k))
						continue;
					}

				double d = lm_dist(ca[i-k], ca[i+k]);

				if (d < best_d)
					{
					best_d = d;
					best_k = k;
					}

				double s = opt.span_close_score_bias - d;
				if (s > best_score)
					best_score = s;
				}

			best_span_dist[i] = best_d;
			best_span_k[i] = best_k;

			if (best_d <= span_feature_threshold)
				span_score[i] = opt.span_close_score_bias - best_d;
			}

		// Generic curvature anomaly.
		if (opt.enable_high_curvature && i >= 1 && i + 1 < n)
			{
			if (curv_deg[i] >= opt.curvature_min_deg)
				curvature_score[i] = curv_deg[i] - opt.curvature_min_deg;
			}

		// Generic torsion anomaly.
		if (opt.enable_torsion_flip && i >= 2 && i + 2 < n)
			{
			double s1 = fabs(tors_deg[i]) - opt.torsion_abs_min_deg;
			double s2 = tors_delta_deg[i] - opt.torsion_delta_min_deg;
			double s = max(s1, s2) * opt.torsion_score_weight;
			if (s > 0)
				torsion_score[i] = s;
			}

		// Strand hairpin turn:
		// center usually coil, with strand runs on both sides, plus span closure.
		//if (opt.enable_strand_hairpin_turn && i >= 1 && i + 1 < n)
		//	{
		//	uint32_t left_strand_run = 0;
		//	uint32_t right_strand_run = 0;

		//	if (i > 0 && ss[i-1] == sstype_t::SS_Strand)
		//		left_strand_run = run_left[i-1];
		//	if (i + 1 < n && ss[i+1] == sstype_t::SS_Strand)
		//		right_strand_run = run_right[i+1];

		//	bool ok_flanks =
		//		left_strand_run >= opt.strand_run_min &&
		//		right_strand_run >= opt.strand_run_min;

		//	double best_d = 1e30;
		//	double best_curv = 0;

		//	for (uint32_t j = i - 1; j <= i + 1; ++j)
		//		{
		//		if (best_span_k[j] != 0 && best_span_dist[j] < best_d)
		//			best_d = best_span_dist[j];
		//		if (curv_deg[j] > best_curv)
		//			best_curv = curv_deg[j];
		//		}

		//	if (ok_flanks &&
		//		best_d <= opt.strand_hairpin_span_threshold &&
		//		best_curv >= opt.strand_hairpin_curvature_min_deg)
		//		{
		//		double s = (opt.strand_hairpin_span_threshold - best_d) +
		//			0.02*(best_curv - opt.strand_hairpin_curvature_min_deg);
		//		if (s > 0)
		//			strand_hairpin_score[i] = s;
		//		if (s > 0) ++n_emit;
		//		////////////////////
		//		}
		//	////////////////////
		//	if (ok_flanks) ++n_ok_flanks;
		//	if (ok_flanks && best_span_k[i] != 0 &&
		//		best_span_dist[i] <= opt.strand_hairpin_span_threshold) ++n_ok_span;
		//	if (ok_flanks && best_span_k[i] != 0 &&
		//		best_span_dist[i] <= opt.strand_hairpin_span_threshold &&
		//		curv_deg[i] >= opt.strand_hairpin_curvature_min_deg) ++n_ok_curv;
		//	////////////////////
		//	}
		if (opt.enable_strand_hairpin_turn && i >= 2 && i + 2 < n)
			{
			bool ok_left = false;
			bool ok_right = false;

			// accept strand within 2 residues on each side
			for (uint32_t j = i - 2; j <= i - 1; ++j)
				{
				if (ss[j] == sstype_t::SS_Strand &&
					run_left[j] >= opt.strand_run_min)
					{
					ok_left = true;
					break;
					}
				}

			for (uint32_t j = i + 1; j <= i + 2; ++j)
				{
				if (ss[j] == sstype_t::SS_Strand &&
					run_right[j] >= opt.strand_run_min)
					{
					ok_right = true;
					break;
					}
				}

			if (ok_left && ok_right)
				{
				double best_d = 1e30;

				// allow the best turn center to drift by 2 residues
				for (uint32_t j = i - 2; j <= i + 2; ++j)
					{
					if (best_span_k[j] != 0 && best_span_dist[j] < best_d)
						best_d = best_span_dist[j];
					}

				if (best_d <= opt.strand_hairpin_span_threshold)
					{
					double s = opt.strand_hairpin_span_threshold - best_d;
					if (s > 0)
						strand_hairpin_score[i] = s;
					}
				}
			}
		// Helix kink:
		// helix runs on both sides, center bent and/or twist changes.
		if (opt.enable_helix_kink && i >= 1 && i + 1 < n)
			{
			uint32_t left_helix_run = 0;
			uint32_t right_helix_run = 0;

			if (i > 0 && ss[i-1] == sstype_t::SS_Helix)
				left_helix_run = run_left[i-1];
			if (i + 1 < n && ss[i+1] == sstype_t::SS_Helix)
				right_helix_run = run_right[i+1];

			bool ok_flanks =
				left_helix_run >= opt.helix_run_min &&
				right_helix_run >= opt.helix_run_min;

			if (ok_flanks)
				{
				double s = 0;
				if (curv_deg[i] >= opt.helix_kink_curvature_min_deg)
					s += curv_deg[i] - opt.helix_kink_curvature_min_deg;
				if (tors_delta_deg[i] >= opt.helix_kink_torsion_delta_min_deg)
					s += tors_delta_deg[i] - opt.helix_kink_torsion_delta_min_deg;
				helix_kink_score[i] = s;
				}
			}

		// SS transitions.
		if (opt.enable_ss_transitions && i + 1 < n)
			{
			if (opt.suppress_across_chain_breaks && !ok_edge[i])
				{
				}
			else
				{
				sstype_t a = ss[i];
				sstype_t b = ss[i+1];
				if (a != b)
					{
					uint32_t boundary_pos = i;

					double geom = 0;
					if (boundary_pos >= 1 && boundary_pos + 1 < n)
						geom += max(0.0, curv_deg[boundary_pos] - opt.curvature_min_deg);
					if (boundary_pos + 1 >= 2 && boundary_pos + 2 < n)
						geom += 0.5*max(0.0, tors_delta_deg[boundary_pos+1] - opt.torsion_delta_min_deg);

					double s = 1.0 + opt.transition_geom_bonus_weight*geom;

					if (a == sstype_t::SS_Helix && b == sstype_t::SS_Coil)
						tr_hc_score[boundary_pos] = s;
					else if (a == sstype_t::SS_Coil && b == sstype_t::SS_Helix)
						tr_ch_score[boundary_pos] = s;
					else if (a == sstype_t::SS_Strand && b == sstype_t::SS_Coil)
						tr_ec_score[boundary_pos] = s;
					else if (a == sstype_t::SS_Coil && b == sstype_t::SS_Strand)
						tr_ce_score[boundary_pos] = s;
					}
				}
			}

		// Local compactness peak.
		if (opt.enable_local_compactness_peak)
			{
			if (i >= opt.compact_window_radius && i + opt.compact_window_radius < n)
				{
				const uint32_t lo = i - opt.compact_window_radius;
				const uint32_t hi = i + opt.compact_window_radius;

				double cnt = 0;
				for (uint32_t a = lo; a <= hi; ++a)
					{
					for (uint32_t b = a + 1; b <= hi; ++b)
						{
						uint32_t sep = b - a;
						if (sep < opt.compact_min_seq_sep)
							continue;
						if (lm_dist(ca[a], ca[b]) <= opt.compact_contact_dist)
							cnt += 1.0;
						}
					}
				if (cnt >= opt.compactness_min)
					compact_score[i] = cnt * opt.compactness_score_weight;
				}
			}

		// Nonlocal contact peak.
		if (opt.enable_nonlocal_contact_peak)
			{
			double cnt = 0;
			for (uint32_t j = 0; j < n; ++j)
				{
				uint32_t sep = (i > j ? i - j : j - i);
				if (sep < opt.min_seq_sep_for_nonlocal)
					continue;
				if (lm_dist(ca[i], ca[j]) <= opt.nonlocal_contact_dist)
					cnt += 1.0;
				}
			if (cnt >= opt.nonlocal_contacts_min)
				nonlocal_score[i] = cnt * opt.nonlocal_score_weight;
			}
		}

	auto allow_pos = [&](const vector<double> &scorev, uint32_t i) -> bool
		{
		if (scorev[i] < opt.min_score)
			return false;
		if (!opt.require_local_maximum)
			return true;
		return lm_is_local_max(scorev, i, opt.nms_radius);
		};

	auto fill_common = [&](landmark_candidate_t &x, uint32_t i)
		{
		x.pos = i;
		x.span_close_score = span_score[i];
		x.best_span_dist = best_span_dist[i];
		x.best_span_k = best_span_k[i];
		x.curvature_deg = curv_deg[i];
		x.torsion_deg = tors_deg[i];
		x.torsion_delta_deg = tors_delta_deg[i];
		x.compactness = compact_score[i];
		x.nonlocal_contacts = nonlocal_score[i];

		x.ss_center = ss[i];
		x.ss_left = (i > 0 ? ss[i-1] : ss[i]);
		x.ss_right = (i + 1 < n ? ss[i+1] : ss[i]);
		x.ss_run_left = run_left[i];
		x.ss_run_right = run_right[i];
		};

	for (uint32_t i = 0; i < n; ++i)
		{
		if (opt.enable_span_closure && span_score[i] > 0 && allow_pos(span_score, i))
			{
			landmark_candidate_t x;
			fill_common(x, i);
			x.cat = landmark_cat_t::LM_SpanClosure;
			x.score = span_score[i];
			out.push_back(x);
			}

		if (opt.enable_high_curvature && curvature_score[i] > 0 && allow_pos(curvature_score, i))
			{
			landmark_candidate_t x;
			fill_common(x, i);
			x.cat = landmark_cat_t::LM_HighCurvature;
			x.score = curvature_score[i];
			out.push_back(x);
			}

		if (opt.enable_torsion_flip && torsion_score[i] > 0 && allow_pos(torsion_score, i))
			{
			landmark_candidate_t x;
			fill_common(x, i);
			x.cat = landmark_cat_t::LM_TorsionFlip;
			x.score = torsion_score[i];
			out.push_back(x);
			}

		if (opt.enable_strand_hairpin_turn && strand_hairpin_score[i] > 0 && allow_pos(strand_hairpin_score, i))
			{
			landmark_candidate_t x;
			fill_common(x, i);
			x.cat = landmark_cat_t::LM_StrandHairpinTurn;
			x.score = strand_hairpin_score[i];
			out.push_back(x);
			}

		if (opt.enable_helix_kink && helix_kink_score[i] > 0 && allow_pos(helix_kink_score, i))
			{
			landmark_candidate_t x;
			fill_common(x, i);
			x.cat = landmark_cat_t::LM_HelixKink;
			x.score = helix_kink_score[i];
			out.push_back(x);
			}

		if (opt.enable_ss_transitions && tr_hc_score[i] > 0 && allow_pos(tr_hc_score, i))
			{
			landmark_candidate_t x;
			fill_common(x, i);
			x.cat = landmark_cat_t::LM_HelixToCoilTransition;
			x.score = tr_hc_score[i];
			x.transition_score = tr_hc_score[i];
			out.push_back(x);
			}

		if (opt.enable_ss_transitions && tr_ch_score[i] > 0 && allow_pos(tr_ch_score, i))
			{
			landmark_candidate_t x;
			fill_common(x, i);
			x.cat = landmark_cat_t::LM_CoilToHelixTransition;
			x.score = tr_ch_score[i];
			x.transition_score = tr_ch_score[i];
			out.push_back(x);
			}

		if (opt.enable_ss_transitions && tr_ec_score[i] > 0 && allow_pos(tr_ec_score, i))
			{
			landmark_candidate_t x;
			fill_common(x, i);
			x.cat = landmark_cat_t::LM_StrandToCoilTransition;
			x.score = tr_ec_score[i];
			x.transition_score = tr_ec_score[i];
			out.push_back(x);
			}

		if (opt.enable_ss_transitions && tr_ce_score[i] > 0 && allow_pos(tr_ce_score, i))
			{
			landmark_candidate_t x;
			fill_common(x, i);
			x.cat = landmark_cat_t::LM_CoilToStrandTransition;
			x.score = tr_ce_score[i];
			x.transition_score = tr_ce_score[i];
			out.push_back(x);
			}

		if (opt.enable_local_compactness_peak && compact_score[i] > 0 && allow_pos(compact_score, i))
			{
			landmark_candidate_t x;
			fill_common(x, i);
			x.cat = landmark_cat_t::LM_LocalCompactnessPeak;
			x.score = compact_score[i];
			out.push_back(x);
			}

		if (opt.enable_nonlocal_contact_peak && nonlocal_score[i] > 0 && allow_pos(nonlocal_score, i))
			{
			landmark_candidate_t x;
			fill_common(x, i);
			x.cat = landmark_cat_t::LM_NonlocalContactPeak;
			x.score = nonlocal_score[i];
			out.push_back(x);
			}
		}

	if (!opt.keep_multiple_categories_per_pos)
		{
		sort(out.begin(), out.end(),
			[](const landmark_candidate_t &a, const landmark_candidate_t &b)
			{
			if (a.pos != b.pos)
				return a.pos < b.pos;
			return a.score > b.score;
			});

		vector<landmark_candidate_t> tmp;
		tmp.reserve(out.size());
		for (const auto &x : out)
			{
			if (tmp.empty() || tmp.back().pos != x.pos)
				tmp.push_back(x);
			}
		out.swap(tmp);
		}

	lm_apply_nms_same_cat(out, opt.nms_radius);

	if (opt.sort_by_score_desc)
		{
		sort(out.begin(), out.end(),
			[](const landmark_candidate_t &a, const landmark_candidate_t &b)
			{
			if (a.score != b.score)
				return a.score > b.score;
			if (a.pos != b.pos)
				return a.pos < b.pos;
			return (uint32_t) a.cat < (uint32_t) b.cat;
			});
		}

	if (opt.max_candidates > 0 && out.size() > opt.max_candidates)
		out.resize(opt.max_candidates);
	}

static void load_msas(
	const string &msafilesfn,
	vector<SeqDB *> &MSAs,
	vector<string> &msastemnames)
	{
	MSAs.clear();
	msastemnames.clear();
	vector<string> msafns;
	ReadLinesFromFile(msafilesfn, msafns);
	const uint n = uint(msafns.size());
	MSAs.reserve(n);
	for (uint i = 0; i < n; ++i)
		{
		ProgressStep(i, n, "Loading MSAs");
		SeqDB *MSA = new SeqDB;
		MSA->FromFasta(msafns[i], true);
		asserta(MSA->IsAligned());
		MSAs.push_back(MSA);

		vector<string> flds;
		Split(msafns[i], flds, '/');
		msastemnames.push_back(flds[flds.size()-1]);
		}
	}

static const char *cat2str(landmark_cat_t cat)
	{
	switch (cat)
		{
	case landmark_cat_t::LM_SpanClosure:			return "spanclo";
	case landmark_cat_t::LM_HighCurvature:			return "highcur";
	case landmark_cat_t::LM_TorsionFlip:			return "torflip";
	case landmark_cat_t::LM_StrandHairpinTurn:		return "sturn";
	case landmark_cat_t::LM_HelixKink:				return "hkink";
	case landmark_cat_t::LM_LocalCompactnessPeak:	return "lcpeak";
	case landmark_cat_t::LM_NonlocalContactPeak:	return "nlcpeak";
	case landmark_cat_t::LM_HelixToCoilTransition:	return "hlx2coil";
	case landmark_cat_t::LM_CoilToHelixTransition:	return "coil2hlx";
	case landmark_cat_t::LM_StrandToCoilTransition:	return "strnd2coil";
	case landmark_cat_t::LM_CoilToStrandTransition:	return "coil2strnd";
		}
	asserta(false);
	return "?";
	}

static const char cat2char(landmark_cat_t cat)
	{
	return g_LetterToCharMu[uint(cat)];
	}

static double get_maxscore(landmark_cat_t cat)
	{
	switch (cat)
		{
	case landmark_cat_t::LM_SpanClosure:			return 8;
	case landmark_cat_t::LM_HighCurvature:			return 85;
	case landmark_cat_t::LM_TorsionFlip:			return 270;
	case landmark_cat_t::LM_StrandHairpinTurn:		return 10;
	case landmark_cat_t::LM_HelixKink:				return 425;
	case landmark_cat_t::LM_LocalCompactnessPeak:	return 21;
	case landmark_cat_t::LM_NonlocalContactPeak:	return 100;
	case landmark_cat_t::LM_HelixToCoilTransition:	return 6;
	case landmark_cat_t::LM_CoilToHelixTransition:	return 6;
	case landmark_cat_t::LM_StrandToCoilTransition:	return 6;
	case landmark_cat_t::LM_CoilToStrandTransition:	return 6;
		}

	asserta(false);
	return 0;
	}

static const uint nbin = 100;
static vector<vector<uint> > s_cat2counts;

static uint score2bin(landmark_cat_t cat, double score)
	{
	double maxscore = get_maxscore(cat);
	uint bin = uint(round(score*nbin/maxscore));
	if (bin > nbin)
		bin = nbin;
	return bin;
	}

static double bin2score(landmark_cat_t cat, uint bin)
	{
	double maxscore = get_maxscore(cat);
	return bin*maxscore/nbin;
	}

static void alloc()
	{
	const uint N = uint(landmark_cat_t::LM_N);
	s_cat2counts.resize(N);
	for (uint i = 0; i < N; ++i)
		{
		landmark_cat_t cat = landmark_cat_t(i);
		s_cat2counts[i].resize(nbin+1);
		}
	}

static void analyze_landmarks_onecol(
	const vector<landmark_candidate_t> &landmarks)
	{
	if (landmarks.empty())
		return;
	const uint n = uint(landmarks.size());
	for (uint i = 0; i < n; ++i)
		{
		const landmark_candidate_t &lm = landmarks[i];
		double score = lm.score;
		landmark_cat_t cat = lm.cat;
		uint bin = score2bin(cat, score);
		s_cat2counts[uint(cat)][bin] += 1;
		}
	}

static void analyze_col2landmarks(
	vector<vector<landmark_candidate_t> > &col2landmarks)
	{
	const uint ncol = uint(col2landmarks.size());
	for (uint colidx = 0; colidx < ncol; ++colidx)
		analyze_landmarks_onecol(col2landmarks[colidx]);
	}

void cmd_landmarks()
	{
	asserta(optset_input);
	asserta(optset_output);
	asserta(optset_output2);
	const string &chainsfn = g_Arg1;
	const string &msafilesfn = opt(input);
	FILE *ffa = CreateStdioFile(opt(fasta));

	vector<SeqDB *> MSAs;
	vector<string> msastemnames;
	load_msas(msafilesfn, MSAs, msastemnames);
	const uint nmsa = uint(MSAs.size());
	ProgressLog("%u msas\n", nmsa);

	alloc();

	landmark_params_t opt;
	opt.enable_strand_hairpin_turn = true;
	opt.enable_helix_kink = true;
	opt.enable_ss_transitions = true;

	opt.enable_span_closure = true;
	opt.enable_high_curvature = true;
	opt.enable_torsion_flip = true;
	opt.enable_local_compactness_peak = true;
	opt.enable_nonlocal_contact_peak = true;

	vector<flat_chain_t *> chains;
	read_flat_chains(chainsfn, chains);
	const uint nchain = uint(chains.size());
	uint16_t *distmx = myalloc(uint16_t, MAXL*M);
	string ss;
	vector<landmark_candidate_t> landmarks;

	unordered_map<string, uint> label2chainidx;
	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		const string &label = chains[chainidx]->m_label;
		label2chainidx[label] = chainidx;
		}

	uint nfound = 0;
	uint notfound = 0;
	vector<uint> pos2col;
	uint total_residues = 0;
	uint total_landmarks = 0;
	uint total_overlaps = 0;
	uint n_coil = 0;
	uint n_helix = 0;
	uint n_strand = 0;
	for (uint msaidx = 0; msaidx < nmsa; ++msaidx)
		{
		ProgressStep(msaidx, nmsa, "Landmarking (%.1f%% not found)",
			GetPct(notfound, nfound + notfound));

		const string &msastemname = msastemnames[msaidx];
		string output_msafn = opt(output2) + msastemname;
		FILE *foutmsa = CreateStdioFile(output_msafn);

		SeqDB &MSA = *MSAs[msaidx];
		const uint nrow = MSA.GetSeqCount();
		const uint ncol = MSA.GetColCount();

		vector<vector<landmark_candidate_t> > col2landmarks(ncol);

		for (uint rowidx = 0; rowidx < nrow; ++rowidx)
			{
			const string &label = MSA.GetLabel(rowidx);
			const string &row = MSA.GetSeq(rowidx);
			unordered_map<string, uint>::const_iterator iter =
				label2chainidx.find(label);
			if (iter == label2chainidx.end())
				{
				++notfound;
				continue;
				}
			++nfound;

			const uint chainidx = iter->second;
			landmarks.clear();

			flat_chain_t *chain = chains[chainidx];
			const uint L = chain->get_length();
			if (L < MINL || L > MAXL)
				continue;

			total_residues += L;
			const ic_t *xyz = chain->m_xyz->m_data;
			chaq::fill_distmx(xyz, L, M, distmx);

			pos2col.clear();
			vector<xyz_t> coords;
			vector<sstype_t> ss;
			coords.reserve(L);
			ss.reserve(L);
			pos2col.reserve(L);

			uint pos = 0;
			string landmark_row;
			landmark_row.resize(ncol, '!');
			for (uint colidx = 0; colidx < ncol; ++colidx)
				{
				char c = row[colidx];
				if (isgap(c))
					landmark_row[colidx] = '-';
				else
					{
					landmark_row[colidx] = 'A';
					pos2col.push_back(colidx);
					}
				}

			string ss_str;
			chaq::get_ss4_str(distmx, M, L, ss_str);
			asserta(ss_str.size() == L);

			for (uint pos = 0; pos < L; ++pos)
				{
				xyz_t c;
				ic_t xa = xyz[3*pos];
				ic_t ya = xyz[3*pos+1];
				ic_t za = xyz[3*pos+2];
				c.x = ic2coord(xa);
				c.y = ic2coord(ya);
				c.z = ic2coord(za);
				coords.push_back(c);

				char ssc = ss_str[pos];
				if (ssc == 'h')
					{
					++n_helix;
					ss.push_back(sstype_t::SS_Helix);
					}
				else if (ssc == 's')
					{
					++n_strand;
					ss.push_back(sstype_t::SS_Strand);
					}
				else
					{
					++n_coil;
					ss.push_back(sstype_t::SS_Coil);
					}
				}

			find_landmark_candidates(coords, ss, opt, landmarks);

			const uint nlm = uint(landmarks.size());
			total_landmarks += nlm;
			if (nlm > 0)
				total_overlaps += 1;

			string landmark_seq;
			landmark_seq.resize(L, 'A');
			for (uint lmidx = 0; lmidx < nlm; ++lmidx)
				{
				const landmark_candidate_t &lm = landmarks[lmidx];
				uint pos = lm.pos;
				asserta(pos < pos2col.size());
				uint colidx = pos2col[pos];
				asserta(colidx < col2landmarks.size());
				col2landmarks[colidx].push_back(lm);

				landmark_seq[pos] = cat2char(lm.cat);
				landmark_row[colidx] = cat2char(lm.cat);
				}
			SeqToFasta(ffa, label, landmark_seq);
			SeqToFasta(foutmsa, label, landmark_row);
			}
		analyze_col2landmarks(col2landmarks);
		CloseStdioFile(foutmsa);
		}
	ProgressLog("total_residues %u, landmarks %u, overlaps %u\n",
		total_residues, total_landmarks, total_overlaps);
	ProgressLog("helix %u, strand %u, coil %u\n",
		n_helix, n_strand, n_coil);
	ProgressLog("n_ok_flanks %u, n_ok_span %u, n_ok_curv %u, n_emit %u\n",
		n_ok_flanks, n_ok_span, n_ok_curv, n_emit);

	const uint ncat = uint(landmark_cat_t::LM_N);
	FILE *f = CreateStdioFile(opt(output));
	fprintf(f, "bin");
	for (uint cati = 0; cati < ncat; ++cati)
		{
		landmark_cat_t cat = landmark_cat_t(cati);
		if (cat == landmark_cat_t::LM_None)
			continue;
		fprintf(f, "\t%s_sc", cat2str(cat));
		fprintf(f, "\t%s_N", cat2str(cat));
		}
	fprintf(f, "\n");

	for (uint binidx = 0; binidx < nbin; ++binidx)
		{
		fprintf(f, "%u", binidx);
		for (uint cati = 0; cati < ncat; ++cati)
			{
			landmark_cat_t cat = landmark_cat_t(cati);
			if (cat == landmark_cat_t::LM_None)
				continue;
			uint count = s_cat2counts[cati][binidx];
			double score = bin2score(cat, binidx);
			fprintf(f, "\t%.2f\t%u", score, count);
			}
		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	CloseStdioFile(ffa);
	}
