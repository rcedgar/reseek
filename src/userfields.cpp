#include "myutils.h"
#include "userfields.h"
#include "dssaligner.h"
#include "cigar.h"

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

static const char *EvalueToStr(double E, string &s)
	{
	if (E > 10)
		E = 99;
	if (E > 1)
		Ps(s, "%.1f", E);
	else if (E > 0.001)
		Ps(s, "%.4f", E);
	else
		Ps(s, "%.3g", E);
	return s.c_str();
	}

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

// Up is true  if alignment is Query=A, Target=B
// Up is false if alignment is Query=B, Target=A
void DSSAligner::WriteUserField(FILE *f, USERFIELD UF, bool aUp)
	{
	if (f == 0)
		return;

	const bool Up = (optset_fast ? aUp : !aUp);

	string TmpStr;
	switch (UF)
		{
	case UF_query:	fputs(GetLabel(Up), f); break;
	case UF_target: fputs(GetLabel(!Up), f); break;
	case UF_evalue:	fprintf(f, "%s", EvalueToStr(GetEvalue(Up), TmpStr)); break;
	case UF_pvalue:	fprintf(f, "%.3g", GetPvalue(Up)); break;
	case UF_ql:		fprintf(f, "%u", GetQL(Up)); break;
	case UF_tl:		fprintf(f, "%u", GetTL(Up)); break;
	case UF_qlo:	fprintf(f, "%u", GetLo(Up) + 1); break;
	case UF_qhi:	fprintf(f, "%u", GetHi(Up) + 1); break;
	case UF_tlo:	fprintf(f, "%u", GetLo(!Up) + 1); break;
	case UF_thi:	fprintf(f, "%u", GetHi(!Up) + 1); break;
	case UF_pctid:	fprintf(f, "%.1f", GetPctId()); break;
	case UF_ts:		fprintf(f, "%.3g", GetTestStatistic(Up)); break;
	case UF_newts:	fprintf(f, "%.3g", GetNewTestStatistic(Up)); break;
	case UF_raw:	fprintf(f, "%.3g", m_AlnFwdScore); break;
	case UF_ids:	fprintf(f, "%u", m_Ids); break;
	case UF_gaps:	fprintf(f, "%u", m_Gaps); break;
	case UF_muscore:	fprintf(f, "%.3g", GetMuScore()); break;

	case UF_cigar:
		{
		string CIGAR;
		PathToCIGAR(m_Path.c_str(), CIGAR, Up);
		fputs(CIGAR.c_str(), f);
		break;
		}

	case UF_qrow:
		{
		string Row;
		GetRow(Up, true, false, Row);
		fputs(Row.c_str(), f);
		break;
		}

	case UF_trow:
		{
		string Row;
		GetRow(Up, false, false, Row);
		fputs(Row.c_str(), f);
		break;
		}

	case UF_qrowg:
		{
		string Row;
		GetRow(Up, true, true, Row);
		fputs(Row.c_str(), f);
		break;
		}

	case UF_trowg:
		{
		string Row;
		GetRow(Up, false, true, Row);
		fputs(Row.c_str(), f);
		break;
		}

	case UF_dpscore:
		{
		fprintf(f, "%.4g", m_AlnFwdScore);
		break;
		}

	case UF_lddt:
		{
		fprintf(f, "%.4g", GetLDDT());
		break;
		}

	case UF_aq:
		{
		fprintf(f, "%.4f", GetAQ(Up));
		break;
		}

	case UF_muhsp:
		{
		fprintf(f, "%d", m_MKF.m_BestHSPScore);
		break;
		}

	case UF_muchain:
		{
		fprintf(f, "%d", m_MKF.m_BestChainScore);
		break;
		}

	case UF_gscore:
		{
		fprintf(f, "%.1f", m_GlobalScore);
		break;
		}

	default:
		Die("Unsupported user field %d='%s'", UF, UFToStr(UF));
		}
	}
