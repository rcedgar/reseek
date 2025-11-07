#include "myutils.h"
#include "peaker.h"

double Peaker::GetIncreaseRateFactor()
	{
	double rf = GetGlobalRateFactor();
	asserta(rf > 1);
	asserta(rf < 3);
	return rf;
	}

double Peaker::GetDecreaseRateFactor()
	{
	double rf = GetGlobalRateFactor();
	asserta(rf > 1);
	asserta(rf < 3);
	return 1.0/rf;
	}

double Peaker::GetGlobalRateFactor()
	{
	if (m_GlobalVarRateFactorIdx == UINT_MAX)
		m_GlobalVarRateFactorIdx = 0;
	string s;
	GetGlobalStr("rates", s, "1.1");
	vector<string> Fields;
	Split(s, Fields, ',');
	const uint n = SIZE(Fields);
	asserta(m_GlobalVarRateFactorIdx < n);
	double Rate = StrToFloat(Fields[m_GlobalVarRateFactorIdx]);
	return Rate;
	}

bool Peaker::ReduceGlobalRateFactor()
	{
	string s;
	GetGlobalStr("rates", s, "1.1");
	vector<string> Fields;
	Split(s, Fields, ',');
	const uint n = SIZE(Fields);
	asserta(m_GlobalVarRateFactorIdx < n);
	if (m_GlobalVarRateFactorIdx + 1 == n)
		{
		ProgressLog("ReduceGlobalRateFactor() no more\n");
		return false;
		}

	++m_GlobalVarRateFactorIdx;
	ProgressLog("\n");
	ProgressLog("ReduceGlobalRateFactor() => %.3f\n", GetGlobalRateFactor());
	ProgressLog("\n");
	return true;
	}

double Peaker::GetRateFactor(bool Plus)
	{
	if (Plus)
		return GetIncreaseRateFactor();
	else
		return GetDecreaseRateFactor();
	}

void Peaker::HJ_Explore()
	{
	const uint VarCount = GetVarCount();

	uint BestNewDirection = UINT_MAX;
	m_HJ_ExtendPlus = false;
	double Best_dy = 0;
	double Start_Best_y = m_Best_y;
	vector<string> strs_plus(VarCount);
	vector<string> strs_minus(VarCount);
	vector<double> dys_plus(VarCount);
	vector<double> dys_minus(VarCount);
	vector<string> Try_xv;
	uint ImprovementCount = 0;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		if (VarIsConstant(VarIdx))
			continue;
		double Saved_Best_y = m_Best_y;
		double y_plus = HJ_TryDelta("explore", m_Best_xv, VarIdx, true, Try_xv);
		double dy_plus = y_plus - Saved_Best_y;
		dys_plus[VarIdx] = dy_plus;
		strs_plus[VarIdx] = Try_xv[VarIdx];
		if (dy_plus > 0)
			++ImprovementCount;
		if (dy_plus > Best_dy)
			{
			Best_dy = dy_plus;
			m_HJ_ExtendPlus = true;
			BestNewDirection = VarIdx;
			continue;
			}

		Saved_Best_y = m_Best_y;
		double y_minus = HJ_TryDelta("explore", m_Best_xv, VarIdx, false, Try_xv);
		double dy_minus = y_minus - Saved_Best_y;
		dys_minus[VarIdx] = dy_minus;
		strs_minus[VarIdx] = Try_xv[VarIdx];
		if (dy_minus > 0)
			++ImprovementCount;
		if (dy_minus > Best_dy)
			{
			Best_dy = dy_minus;
			m_HJ_ExtendPlus = false;
			BestNewDirection = VarIdx;
			}
		}
	m_HJ_Direction = BestNewDirection;
	if (m_HJ_Direction == UINT_MAX)
		{
		Log("HJ_Explore(), no improvement found\n");
		return;
		}
	Log("HJ_Explore(), new direction %s%c\n",
		GetVarName(m_HJ_Direction), pom(m_HJ_ExtendPlus));

	ProgressLog("\n");
	ProgressLog("HJ_Explore %u improves\n", ImprovementCount);
	double Track_Best_y = Start_Best_y;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		ProgressLog(">%-10.10s", strs_plus[VarIdx].c_str());
		ProgressLog("  %10.2g", dys_plus[VarIdx]);
		ProgressLog("  <%-10.10s", strs_minus[VarIdx].c_str());
		ProgressLog("  %10.2g", dys_minus[VarIdx]);
		ProgressLog("  %s", GetVarName(VarIdx));
		if (dys_plus[VarIdx] > 0 || dys_minus[VarIdx] > 0)
			ProgressLog(" +++");
		ProgressLog("\n");
		}
	ProgressLog("\n");

	const uint N = SIZE(m_Best_ys);
	asserta(SIZE(m_Best_descs) == N);
	for (uint k = 0; k < 5; ++k)
		{
		if (k > N -1)
			break;
		uint i = N-k-1;
		double dy = (i > 0 ? m_Best_ys[i] - m_Best_ys[i-1] : 0);
		ProgressLog("%10.5g", m_Best_ys[i]);
		ProgressLog("  %+10.2g", dy);
		ProgressLog("  %s", m_Best_descs[i].c_str());
		ProgressLog("\n");
		}

	double rf = GetGlobalRateFactor();
	ProgressLog("Global rate factor %.2f\n", GetGlobalRateFactor());
	ProgressLog("\n");
	}

void Peaker::HJ_Extend()
	{
	if (m_HJ_Direction == UINT_MAX)
		return;
	const uint VarCount = GetVarCount();
	asserta(m_HJ_Direction < VarCount);
	asserta(!VarIsConstant(m_HJ_Direction));
	const uint VarIdx = m_HJ_Direction;
	const char *Name = GetVarName(VarIdx);
	vector<string> Try_xv;
	double Start_Best_y = m_Best_y;
	for (uint Iter = 0; Iter < m_HJ_MaxExtendIters; ++Iter)
		{
		string reason;
		Ps(reason, "extend%u", Iter+1);

		double Saved_Best_y = m_Best_y;
		HJ_TryDelta(reason, m_Best_xv, VarIdx, m_HJ_ExtendPlus, Try_xv);
		if (m_Best_y <= Saved_Best_y)
			break;
		}
	}

double Peaker::HJ_TryDelta(const string &reason,
	const vector<string> &Start_xv, uint VarIdx, bool Plus,
	vector<string> &Try_xv)
	{
	const char *VarName = GetVarName(VarIdx);
	uint Idx = Find_xv(Start_xv);
	asserta(Idx != UINT_MAX);
	asserta(Idx < SIZE(m_ys));
	const double Start_y = m_ys[Idx];
	vector<string> Saved_Best_xv;
	vector<string> Saved_Start_xv;
	const uint VarCount = GetVarCount();
	for (uint i = 0; i < VarCount; ++i)
		{
		Saved_Start_xv.push_back(Start_xv[i]);
		Saved_Best_xv.push_back(m_Best_xv[i]);
		}

	string NewStr;
	const string OldStr = Start_xv[VarIdx];
	double OldValue = VarStrToFloat(VarIdx, OldStr);
	DeltaVar(VarIdx, Plus, OldStr, NewStr);
	if (NewStr == OldStr)
		{
		ProgressLog("HJ_TryDelta(%s%c) DeltaVar %s=%s no change\n",
			reason.c_str(), pom(Plus), VarName, OldStr.c_str());
		return Start_y;
		}

	Try_xv.clear();
	for (uint i = 0; i < VarCount; ++i)
		Try_xv.push_back(Start_xv[i]);
	Try_xv[VarIdx] = NewStr;

	string why;
	Ps(why, "%s%c%s", reason.c_str(), pom(Plus), VarName);
	double y = Evaluate(Try_xv, why);
	if (y == DBL_MAX)
		return DBL_MAX;
	if (Start_y == DBL_MAX)
		return y;

	double absdy = fabs(Start_y - y);
	double targetdy = GetGlobalFloat("targetdy", DBL_MAX);
	asserta(targetdy != DBL_MAX);

	Log("HJ_TryDelta(%s%c)", reason.c_str(), pom(Plus));
	Log(" %s", VarName);
	Log(" %s,", OldStr.c_str());
	Log(" %s", NewStr.c_str());
	Log(" y %.4g,", Start_y);
	Log("%.4g", y);
	Log(" dy %.3g (target %.3g)", absdy, targetdy);
	Log("\n");

	return y;
	}

bool Peaker::HJ_Iter()
	{
	double Saved_Best_y = m_Best_y;
	HJ_Explore();
	double Height = m_Best_y - Saved_Best_y;
	if (Height > 0)
		{
		if (GetGlobalBool("extend", false))
			HJ_Extend();
		return true;
		}
	bool ok = ReduceGlobalRateFactor();
	return ok;
	}

void Peaker::HJ_RunHookeJeeves()
	{
	InitRates();
	for (uint Iters = 0; ; ++Iters)
		{
		if (Iters >= m_HJ_MaxIters)
			{
			Warning("HJ max iters, not converged");
			break;
			}
		bool ok = HJ_Iter();
		if (!ok)
			return;
		}
	}

void Peaker::NormalizeVarStr(uint VarIdx, const string &Str,
	string &NormalizedStr) const
	{
	double Value = StrToFloat(Str);
	double z = VarSpecGetFloat(VarIdx, "zero", 0);
	if (Value < z)
		Value = 0;
	uint SigFig = VarSpecGetInt(VarIdx, "sigfig", 2);
	GetRoundedStr(Value, SigFig, NormalizedStr);
	}

void Peaker::DeltaVar(uint VarIdx, bool Plus,
	const string &OldStr, string &NewStr)
	{
	NewStr.clear();
	uint SigFig = VarSpecGetInt(VarIdx, "sigfig", 2);
	double OldValue = VarStrToFloat(VarIdx, OldStr);
	if (OldValue == 0)
		{
		if (!Plus)
			{
			NewStr = OldStr;
			return;
			}
		double IfZero = VarSpecGetFloat(VarIdx, "ifzero", 0);
		string TmpStr;
		VarFloatToStr(VarIdx, IfZero, TmpStr);
		NormalizeVarStr(VarIdx, TmpStr, NewStr);
		ProgressLog("%s ifzero %s => %s\n", GetVarName(VarIdx),
			OldStr.c_str(), NewStr.c_str());
		return;
		}
	string TmpStr;
	double Factor = GetRateFactor(Plus);
	double NewValue = OldValue*Factor;
	VarFloatToStr(VarIdx, NewValue, TmpStr);
	NormalizeVarStr(VarIdx, TmpStr, NewStr);
	if (NewStr == OldStr)
		{
		IncFloat(OldStr, Plus, TmpStr);
		NormalizeVarStr(VarIdx, TmpStr, NewStr);
		}
	}
