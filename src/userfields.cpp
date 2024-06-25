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

void DSSAligner::AppendUserField(string &s, USERFIELD UF, bool IsA)
	{
	if (!s.empty())
		s += '\t';

	switch (UF)
		{
	case UF_query:
		{
		s += GetLabelA(IsA);
		break;
		}

	case UF_target:
		{
		s += GetLabelB(IsA);
		break;
		}

	case UF_evalue:
		{
		Psa(s, "%.3g", GetEvalueA(IsA));
		break;
		}

	case UF_qlo:
		{
		Psa(s, "%.3g", GetLoA(IsA));
		break;
		}

	case UF_qhi:
		{
		Psa(s, "%.3g", GetHiA(IsA));
		break;
		}

	case UF_tlo:
		{
		Psa(s, "%.3g", GetLoB(IsA));
		break;
		}

	case UF_thi:
		{
		Psa(s, "%.3g", GetHiB(IsA));
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
		PathToCIGAR(m_PathA.c_str(), CIGAR, IsA);
		s += CIGAR;
		break;
		}

	case UF_qrow:
		{
		string Row;
		if (IsA) GetRowB(Row); else GetRowA(Row);
		s += Row;
		break;
		}

	case UF_trow:
		{
		string Row;
		if (IsA) GetRowA(Row); else GetRowB(Row);
		s += Row;
		break;
		}

	case UF_qrowg:
		{
		string Row;
		if (IsA) GetRowBg(Row); else GetRowAg(Row);
		s += Row;
		break;
		}

	case UF_trowg:
		{
		string Row;
		if (IsA) GetRowAg(Row); else GetRowBg(Row);
		s += Row;
		break;
		}

	case UF_ts:
		{
		Psa(s, "%.3g", IsA ? m_TestStatisticB : m_TestStatisticA);
		break;
		}

	case UF_rigid:
		{
		double Kabsch(const PDBChain &ChainA, const PDBChain &ChainB,
		  uint LoA, uint LoB, const string &Path,
		  vector<double> &t, vector<vector<double> > &u);

		vector<double> t;
		vector<vector<double> > R;
		double RMS = DBL_MAX;
		if (IsA)
			RMS = Kabsch(*m_ChainA, *m_ChainB, m_LoA, m_LoB, m_PathA, t, R);
		else
			{
			void InvertPath(const string &Path, string &InvPath);
			string PathBA;
			InvertPath(m_PathA, PathBA);
			RMS = Kabsch(*m_ChainB, *m_ChainA, m_LoB, m_LoA, PathBA, t, R);
			}
		Ps(s, "%.1f", RMS);
		asserta(SIZE(t) == 3);
		asserta(SIZE(R) == 3);
		asserta(SIZE(R[0]) == 3);
		Psa(s, "\t%.4g\t%.4g\t%.4g", t[0], t[1], t[2]);
		for (uint i = 0; i < 3; ++i)
			Psa(s, "\t%.4g\t%.4g\t%.4g", R[i][0], R[i][1], R[i][2]);
		}

	default:
		Die("Unsupported user field %d='%s'", UF, UFToStr(UF));
		}
	}
