#include "myutils.h"
#include "reseeker.h"
#include "hitdata.h"

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
	if (m_faln == 0) return;
	Die("reseeker::write_aln() TODO");
	}

void reseeker::write_tsv(const hitdata &hit)
	{
	if (m_fhit == 0) return;
	string str;
	str = hit.query->m_label;
	str += "\t" + hit.target->m_label;
	Psa(str, "\t%.3g", hit.TS);
	str += "\n";
	// fprintf & fputs are thread-safe
	fputs(str.c_str(), m_fhit);
	}

void reseeker::init_userfields()
	{
	static const vector<USERFIELD> default_columns = 
		{ UF_query, UF_target, UF_raw };

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
			Ps(s, "%.4f", P);
		else
			Ps(s, "%.4g", P);
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
		asserta(false);
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
		Ps(s, "%.1f", GetPct(hit.ncol, hit.ids));
		break;
		}

	case UF_cigar:
		{
		s += hit.cigar;
		break;
		}

	case UF_qrow:
		{
		string row;
		row.reserve(hit.ncol);
		uint qpos = hit.qlo;
		const char *qaa = hit.query->m_aa->m_data;
		for (uint i = 0; i < hit.ncol; ++i)
			{
			char c = hit.path[i];
			if (c == 'M' || c == 'D')
				row += qaa[qpos++];
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
		const char *taa = hit.target->m_aa->m_data;
		for (uint i = 0; i < hit.ncol; ++i)
			{
			char c = hit.path[i];
			if (c == 'M' || c == 'I')
				row += taa[tpos++];
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
