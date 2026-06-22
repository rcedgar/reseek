#include "myutils.h"
#include "hitdata.h"
#include "cigar.h"
#include "flat_chain.h"

string flat_chain_aa_string(const flat_chain_t *chain, uint L)
	{
	asserta(chain != 0);
	string s;
	s.reserve(L);
	for (uint i = 0; i < L; ++i)
		s += chain->get_aa(i);
	return s;
	}

void reseek_hit_fill_aa_cigar(
	reseek_hit &hit,
	const flat_chain_t *query_chain,
	const flat_chain_t *target_chain)
	{
	const uint LQ = hit.query.slice_len;
	const uint LT = hit.target.slice_len;
	asserta(LQ > 0);
	asserta(LT > 0);
	hit.seq_q = flat_chain_aa_string(query_chain, LQ);
	hit.seq_t = flat_chain_aa_string(target_chain, LT);
	if (hit.fwd.ncol > 0 && !hit.fwd.path.empty())
		PathToCIGAR(hit.fwd.path.c_str(), hit.cigar);
	else
		hit.cigar.clear();
	}

void reseek_hit_emit_tsv(FILE *f, bool nu_only, const reseek_hit &hit)
	{
	if (f == 0)
		return;

	string str;
	str = hit.parent_label_q;
	str += "\t" + hit.parent_label_t;
	if (nu_only)
		{
		Psa(str, "\t%.3g", hit.nu_combined_score);
		str += "\n";
		fputs(str.c_str(), f);
		return;
		}

	Psa(str, "\t%.3g", hit.TS);
	Psa(str, "\t%.3g", hit.nu_fwd_score);
	Psa(str, "\t%.3g", hit.nu_combined_score);
	Psa(str, "\t%u", hit.kappa_diag_score);
	str += "\n";
	fputs(str.c_str(), f);
	}

void reseek_hit_emit_aln_aa(FILE *f, const reseek_hit &hit)
	{
	(void) hit;
	if (f == 0)
		return;
	// TODO: stitched aa-only pretty alignment from cigar + seq_q/seq_t + slices
	}
