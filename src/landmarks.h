#pragma once

static const double MY_PI = 3.1415926535;

enum class sstype_t : uint8_t
	{
	SS_Coil = 0,
	SS_Helix,
	SS_Strand
	};

enum class landmark_cat_t : uint8_t
	{
	LM_None = 0,

	// Geometry-first generic classes
	LM_SpanClosure,
	LM_HighCurvature,
	LM_TorsionFlip,

	// Secondary-structure-informed classes
	LM_StrandHairpinTurn,
	LM_HelixKink,
	LM_HelixToCoilTransition,
	LM_CoilToHelixTransition,
	LM_StrandToCoilTransition,
	LM_CoilToStrandTransition,

	// Contextual classes
	LM_LocalCompactnessPeak,
	LM_NonlocalContactPeak,

	// for arrays
	LM_N
	};

struct xyz_t
	{
	double x = 0;
	double y = 0;
	double z = 0;
	};

struct landmark_candidate_t
	{
	uint32_t pos = 0;
	landmark_cat_t cat = landmark_cat_t::LM_None;

	double score = 0;

	// Raw features for later calibration / sweeps.
	double span_close_score = 0;
	double best_span_dist = 0;
	uint32_t best_span_k = 0;

	double curvature_deg = 0;
	double torsion_deg = 0;
	double torsion_delta_deg = 0;

	double compactness = 0;
	double nonlocal_contacts = 0;
	double transition_score = 0;

	// Secondary structure context.
	sstype_t ss_center = sstype_t::SS_Coil;
	sstype_t ss_left = sstype_t::SS_Coil;
	sstype_t ss_right = sstype_t::SS_Coil;
	uint32_t ss_run_left = 0;   // same-state run length ending at pos
	uint32_t ss_run_right = 0;  // same-state run length starting at pos

	uint32_t aux = 0;
	};

struct landmark_params_t
	{
	// Enabled detectors.
	bool enable_span_closure = true;
	bool enable_high_curvature = true;
	bool enable_torsion_flip = true;

	bool enable_strand_hairpin_turn = true;
	bool enable_helix_kink = true;
	bool enable_ss_transitions = true;

	bool enable_local_compactness_peak = true;
	bool enable_nonlocal_contact_peak = true;

	// Emission / dedup.
	bool keep_multiple_categories_per_pos = true;
	bool require_local_maximum = true;
	bool sort_by_score_desc = true;

	uint32_t nms_radius = 2;        // same-category suppression radius
	uint32_t max_candidates = 0;    // 0 => no limit
	double min_score = 0;

	// Span closure.
	uint32_t span_k_min = 4;
	uint32_t span_k_max = 10;
	double span_close_threshold = 5.0;   // Angstroms
	double span_close_score_bias = 8.0;  // score = bias - distance

	// Curvature.
	double curvature_min_deg = 95.0;

	// Torsion anomaly.
	double torsion_abs_min_deg = 60.0;
	double torsion_delta_min_deg = 90.0;
	double torsion_score_weight = 1.0;

	// Strand hairpin turn:
	// require strand runs on both sides and a strong span closure.
	uint32_t strand_run_min = 2;
	double strand_hairpin_span_threshold = 6.0;
	double strand_hairpin_curvature_min_deg = 70.0;

	// Helix kink:
	// require helix runs on both sides and a bend / torsion anomaly.
	uint32_t helix_run_min = 2;
	double helix_kink_curvature_min_deg = 70.0;
	double helix_kink_torsion_delta_min_deg = 45.0;

	// SS transitions:
	// score can include geometry anomaly near boundary.
	double transition_geom_bonus_weight = 0.02;

	// Local compactness.
	uint32_t compact_window_radius = 4;  // [i-r, i+r]
	uint32_t compact_min_seq_sep = 3;
	double compact_contact_dist = 8.0;
	double compactness_min = 6.0;
	double compactness_score_weight = 1.0;

	// Nonlocal contacts.
	uint32_t min_seq_sep_for_nonlocal = 8;
	double nonlocal_contact_dist = 10.0;
	double nonlocal_contacts_min = 3.0;
	double nonlocal_score_weight = 1.0;

	// Optional chain break handling:
	// treat adjacent residues farther apart than this as broken.
	bool suppress_across_chain_breaks = true;
	double max_adjacent_ca_dist = 4.5;
	};

static inline double lm_dot(const xyz_t &a, const xyz_t &b)
	{
	return a.x*b.x + a.y*b.y + a.z*b.z;
	}

static inline xyz_t lm_sub(const xyz_t &a, const xyz_t &b)
	{
	xyz_t r;
	r.x = a.x - b.x;
	r.y = a.y - b.y;
	r.z = a.z - b.z;
	return r;
	}

static inline xyz_t lm_cross(const xyz_t &a, const xyz_t &b)
	{
	xyz_t r;
	r.x = a.y*b.z - a.z*b.y;
	r.y = a.z*b.x - a.x*b.z;
	r.z = a.x*b.y - a.y*b.x;
	return r;
	}

static inline double lm_norm(const xyz_t &a)
	{
	return std::sqrt(lm_dot(a, a));
	}

static inline double lm_dist(const xyz_t &a, const xyz_t &b)
	{
	return lm_norm(lm_sub(a, b));
	}

static inline double lm_clamp(double x, double lo, double hi)
	{
	return x < lo ? lo : (x > hi ? hi : x);
	}

static inline double lm_absdiff(double a, double b)
	{
	double d = a - b;
	return d >= 0 ? d : -d;
	}

static inline double lm_angle_deg(const xyz_t &a, const xyz_t &b, const xyz_t &c)
	{
	xyz_t u = lm_sub(a, b);
	xyz_t v = lm_sub(c, b);
	double nu = lm_norm(u);
	double nv = lm_norm(v);
	if (nu <= 0 || nv <= 0)
		return 0;
	double cs = lm_dot(u, v)/(nu*nv);
	cs = lm_clamp(cs, -1.0, 1.0);
	return std::acos(cs)*180.0/MY_PI;
	}

static inline double lm_dihedral_deg(
	const xyz_t &a,
	const xyz_t &b,
	const xyz_t &c,
	const xyz_t &d)
	{
	xyz_t b1 = lm_sub(b, a);
	xyz_t b2 = lm_sub(c, b);
	xyz_t b3 = lm_sub(d, c);

	xyz_t n1 = lm_cross(b1, b2);
	xyz_t n2 = lm_cross(b2, b3);

	double n1n = lm_norm(n1);
	double n2n = lm_norm(n2);
	double b2n = lm_norm(b2);
	if (n1n <= 0 || n2n <= 0 || b2n <= 0)
		return 0;

	xyz_t b2u;
	b2u.x = b2.x / b2n;
	b2u.y = b2.y / b2n;
	b2u.z = b2.z / b2n;

	xyz_t m1 = lm_cross(n1, b2u);

	double x = lm_dot(n1, n2) / (n1n*n2n);
	double y = lm_dot(m1, n2) / (n1n*n2n);

	return std::atan2(y, x)*180.0/MY_PI;
	}

static inline bool lm_is_local_max(const std::vector<double> &v, uint32_t i, uint32_t radius)
	{
	const uint32_t n = (uint32_t) v.size();
	const double x = v[i];
	const uint32_t lo = (i > radius ? i - radius : 0);
	const uint32_t hi = std::min<uint32_t>(n - 1, i + radius);
	for (uint32_t j = lo; j <= hi; ++j)
		{
		if (j == i)
			continue;
		if (v[j] > x)
			return false;
		}
	return true;
	}

static inline void lm_apply_nms_same_cat(
	std::vector<landmark_candidate_t> &cand,
	uint32_t radius)
	{
	if (cand.empty() || radius == 0)
		return;

	std::sort(cand.begin(), cand.end(),
		[](const landmark_candidate_t &a, const landmark_candidate_t &b)
		{
		if (a.cat != b.cat)
			return (uint32_t) a.cat < (uint32_t) b.cat;
		if (a.score != b.score)
			return a.score > b.score;
		return a.pos < b.pos;
		});

	std::vector<landmark_candidate_t> out;
	out.reserve(cand.size());

	for (const auto &x : cand)
		{
		bool keep = true;
		for (const auto &y : out)
			{
			if (x.cat != y.cat)
				continue;
			uint32_t d = (x.pos > y.pos ? x.pos - y.pos : y.pos - x.pos);
			if (d <= radius)
				{
				keep = false;
				break;
				}
			}
		if (keep)
			out.push_back(x);
		}

	cand.swap(out);
	}

static inline uint32_t lm_same_state_run_left(
	const std::vector<sstype_t> &ss,
	uint32_t i)
	{
	const sstype_t s = ss[i];
	uint32_t n = 0;
	for (;;)
		{
		if (ss[i] != s)
			break;
		++n;
		if (i == 0)
			break;
		--i;
		}
	return n;
	}

static inline uint32_t lm_same_state_run_right(
	const std::vector<sstype_t> &ss,
	uint32_t i)
	{
	const uint32_t N = (uint32_t) ss.size();
	const sstype_t s = ss[i];
	uint32_t n = 0;
	while (i < N && ss[i] == s)
		{
		++n;
		++i;
		}
	return n;
	}

static inline bool lm_has_chain_break(
	const std::vector<uint8_t> &ok_edge,
	uint32_t lo,
	uint32_t hi)
	{
	// Tests edges [lo..hi-1]
	for (uint32_t i = lo; i < hi; ++i)
		{
		if (!ok_edge[i])
			return true;
		}
	return false;
	}
