#include "myutils.h"
#include "reseeker.h"
#include "hitdata.h"
#include "flat_helpers.h"

/***
Called from mulitple threads.
Output functions are responsible for locking.
For single-line output e.g. tsv DO NOT lock, exploit thread-safety of fputs().
For multi-line output e.g. aln use mutex or equivalent.
***/

const char *UFToStr(USERFIELD UF)
	{
	switch (UF)
		{
#define x(name)	case UF_##name : return #name;
#include "userfieldnames.h"
		}
	Die("Invalid USERFIELD=%d", UF);
	return "UF_ERROR";
	}

USERFIELD StrToUF(const char *Str)
	{
	return StrToUF(string(Str));
	}

USERFIELD StrToUF(const string &Str)
	{
#define x(name)	if (Str == #name) return UF_##name;
#include "userfieldnames.h"
	Die("Invalid user field name '%s'", Str.c_str());
	return UF_Undefined;
	}

static const char *PvalueToStr(double P, string &s)
	{
	if (P > 10)
		P = 99;
	if (P > 1)
		Ps(s, "%.1f", P);
	else if (P > 0.001)
		Ps(s, "%.4f", P);
	else
		Ps(s, "%.3g", P);
	return s.c_str();
	}

void reseeker::write_aln(const hitdata &hit)
	{
	if (m_faln == 0)
		return;

	lock_guard<mutex> lock(m_aln_lock);

	fprintf(m_faln, "\n");
	fprintf(m_faln, "_____________________________________________________________________________________________________________\n");

	string labelQ, labelT;
	WriteLocalAln(m_faln,
		hit.query->m_label, (const byte *) hit.query->m_aa->m_data,
		hit.target->m_label, (const byte *) hit.target->m_aa->m_data,
		hit.qlo, hit.tlo,
		hit.path);

	fprintf(m_faln, "%s [%u-%u/%u aa]\n",
		hit.query->m_label.c_str(), hit.qlo + 1, hit.qhi + 1, hit.query->m_L);
	fprintf(m_faln, "%s [%u-%u/%u aa]\n",
		hit.target->m_label.c_str(), hit.tlo + 1, hit.thi + 1, hit.target->m_L);

	string s;
	fprintf(m_faln, "P-value %s, cols %u, gaps %u, ids %u (%.1f%%)\n",
		PvalueToStr(hit.pvalue, s), hit.ncol, hit.gaps, hit.ids,
		GetPct(hit.ids, hit.ncol));
	}

void reseeker::write_tsv(const hitdata &hit)
	{
	if (m_fhit == 0) return;
	string str;
	for (auto uf : m_UFs)
		append_userfield(str, hit, uf);
	fputs(str.c_str(), m_fhit);
	}

void reseeker::init_userfields()
	{
	static const vector<USERFIELD> default_columns = 
		{ UF_query, UF_target, UF_pvalue };

	if (optset_columns)
		{
		vector<string> Fields;
		Split(string(opt(columns)), Fields, '+');
		const uint n = SIZE(Fields);
		if (n == 0)
			Die("Empty -columns option");
		for (uint i = 0; i < n; ++i)
			{
			if (Fields[i] == "std")
				m_UFs = default_columns;
			else
				{
				USERFIELD UF = StrToUF(Fields[i]);
				m_UFs.push_back(UF);
				}
			}
		}
	else
		m_UFs = default_columns;
	}

void reseeker::append_userfield(
	string &s,
	const hitdata &hit,
	USERFIELD UF)
	{
	if (!s.empty()) s += '\t';
	switch (UF)
		{
	case UF_query:
		{
		s += hit.query->m_label;
		break;
		}

	case UF_target:
		{
		s += hit.target->m_label;
		break;
		}

	case UF_pvalue:
		{
		const double P = hit.pvalue;
		if (P >= 1)
			s += "1.000";
		else if (P > 0.001)
			Psa(s, "%.4f", P);
		else
			Psa(s, "%.4g", P);
		break;
		}

	case UF_evalue:
		{
		const double SCOP40X_SIZE = 8291;
		const double E = hit.pvalue*SCOP40X_SIZE;
		if (E >= 100)
			s += "100";
		else if (E > 1)
			Ps(s, "%.2f", E);
		else if (E > 0.001)
			Ps(s, "%.4f", E);
		else
			Ps(s, "%.4g", E);
		break;
		}

	case UF_qlo:
		{
		Psa(s, "%u", hit.qlo+1);
		break;
		}

	case UF_qhi:
		{
		Psa(s, "%u", hit.qhi+1);
		break;
		}

	case UF_tlo:
		{
		Psa(s, "%u", hit.tlo+1);
		break;
		}

	case UF_thi:
		{
		Psa(s, "%u", hit.thi+1);
		break;
		}

	case UF_ql:
		{
		Psa(s, "%u", hit.query->m_L);
		break;
		}

	case UF_tl:
		{
		Psa(s, "%u", hit.target->m_L);
		break;
		}

	case UF_pctid:
		{
		Ps(s, "%.1f", GetPct(hit.ids, hit.ncol));
		break;
		}

	case UF_cigar:
		{
		s.append(hit.cigar_ptr(), hit.cigar_length());
		break;
		}

	case UF_qrow:
		{
		string row;
		row.reserve(hit.ncol);
		uint qpos = hit.qlo;
		const uint LQ = hit.query->m_L;
		const char *qaa = hit.query->m_aa->m_data;
		for (uint i = 0; i < hit.ncol; ++i)
			{
			char c = hit.path[i];
			if (c == 'M' || c == 'I')
				{
				assert(qpos < LQ);
				row += qaa[qpos++];
				}
			else
				row += '-';
			}
		s += row;
		break;
		}

	case UF_trow:
		{
		string row;
		row.reserve(hit.ncol);
		uint tpos = hit.tlo;
		const uint LT = hit.target->m_L;
		const char *taa = hit.target->m_aa->m_data;
		for (uint i = 0; i < hit.ncol; ++i)
			{
			char c = hit.path[i];
			if (c == 'M' || c == 'D')
				{
				assert(tpos < LT);
				row += taa[tpos++];
				}
			else
				row += '-';
			}
		s += row;
		break;
		}

	case UF_raw:
		{
		Psa(s, "%.3g", hit.TS);
		break;
		}

	case UF_mega:
		{
		Psa(s, "%.3g", hit.mega_fwd_score);
		break;
		}

	case UF_lddt:
		{
		Psa(s, "%.3g", hit.lddt);
		break;
		}

	case UF_dali:
		{
		Psa(s, "%.3g", hit.dali);
		break;
		}

	case UF_tm:
		{
		Psa(s, "%.3g", hit.TM);
		break;
		}

	case UF_ids:
		{
		Psa(s, "%u", hit.ids);
		break;
		}

	case UF_gaps:
		{
		Psa(s, "%u", hit.gaps);
		break;
		}

	case UF_aq:
		{
		Die("AQ TODO");
		break;
		}

	case UF_cols:
		{
		Psa(s, "%u", hit.ncol);
		break;
		}

	case UF_qcovpct:
		{
		uint QL = hit.query->m_L;
		uint qsegl = hit.qhi - hit.qlo + 1;
		double pct = GetPct(qsegl, QL);
		Psa(s, "%.1f", pct);
		break;
		}

	case UF_tcovpct:
		{
		uint TL = hit.target->m_L;
		uint tsegl = hit.thi - hit.tlo + 1;
		double pct = GetPct(tsegl, TL);
		Psa(s, "%.1f", pct);
		break;
		}

	default: Die("Bad column index %u (%s)", uint(UF), UFToStr(UF));
		}
	}
