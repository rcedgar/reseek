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

void DSSAligner::AppendUserField(string &s, USERFIELD UF, bool IsBA)
	{
	if (!s.empty())
		s += '\t';

	switch (UF)
		{
	case UF_query:
		{
		s += m_ChainA->m_Label;
		break;
		}

	case UF_target:
		{
		s += m_ChainB->m_Label;
		break;
		}

	case UF_evalue:
		{
		Psa(s, "%.3g", IsBA ? m_EvalueBA : m_EvalueAB);
		break;
		}

	case UF_qlo:
		{
		Psa(s, "%.3g", IsBA ? m_LoA + 1 : m_LoB + 1);
		break;
		}

	case UF_qhi:
		{
		//Psa(s, "%.3g", IsBA ? m_HiA + 1 : m_HiB + 1);
		break;
		}

	case UF_tlo:
		{
		Psa(s, "%.3g", IsBA ? m_LoB + 1 : m_LoA + 1);
		break;
		}

	case UF_thi:
		{
		//Psa(s, "%.3g", IsBA ? m_HiB + 1 : m_HiA + 1);
		break;
		}

	case UF_pctid:
		{
		Psa(s, "%.1f", GetPctId());
		break;
		}

	case UF_cigar:
		{
		string CIGAR;
		PathToCIGAR(m_PathAB.c_str(), CIGAR, IsBA);
		s += CIGAR;
		break;
		}

	case UF_qrow:
		{
		string QRow;
		break;
		}

	case UF_trow:
		{
		break;
		}

	case UF_qrowg:
		{
		break;
		}

	case UF_trowg:
		{
		break;
		}

	case UF_ts:
		{
		break;
		}

	default:
		Die("Unsupported user field %d='%s'", UF, UFToStr(UF));
		}
	}
