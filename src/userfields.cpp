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
void DSSAligner::WriteUserField(FILE *f, USERFIELD UF, bool Up) const
	{
	if (f == 0)
		return;

	switch (UF)
		{
	case UF_query:	fputs(GetLabel(Up), f); break;
	case UF_target: fputs(GetLabel(!Up), f); break;
	case UF_evalue:	fprintf(f, "%.3g", GetEvalue(Up)); break;
	case UF_qlo:	fprintf(f, "%u", GetLo(Up) + 1); break;
	case UF_qhi:	fprintf(f, "%u", GetHi(Up) + 1); break;
	case UF_tlo:	fprintf(f, "%u", GetLo(!Up) + 1); break;
	case UF_thi:	fprintf(f, "%u", GetHi(!Up) + 1); break;
	case UF_pctid:	fprintf(f, "%.1f", GetPctId()); break;
	case UF_ts:		fprintf(f, "%.3g", GetTestStatistic(Up)); break;

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

	case UF_rigid:
		{
		vector<double> t;
		vector<vector<double> > R;
		float RMS = GetKabsch(t, R, Up);
		fprintf(f, "%.1f", RMS);
		asserta(SIZE(t) == 3);
		asserta(SIZE(R) == 3);
		asserta(SIZE(R[0]) == 3);
		fprintf(f, "\t%.4g\t%.4g\t%.4g", t[0], t[1], t[2]);
		for (uint i = 0; i < 3; ++i)
			fprintf(f, "\t%.4g\t%.4g\t%.4g", R[i][0], R[i][1], R[i][2]);
		break;
		}

	default:
		Die("Unsupported user field %d='%s'", UF, UFToStr(UF));
		}
	}
