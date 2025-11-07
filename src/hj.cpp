#include "myutils.h"
#include "peaker.h"

double Peaker::GetIncreaseRateFactor(uint Rate) const
	{
	double rf = GetGlobalRateFactor();
	if (rf != DBL_MAX)
		{
		asserta(rf > 1);
		asserta(rf < 3);
		return rf;
		}
	if (optset_rate_factor)
		return opt(rate_factor);

	asserta(Rate >= MIN_RATE && Rate <= MAX_RATE);
	switch (Rate)
		{
	case 1:	return 1.02 + (opt(ratenoise) ? randf(0.02) : 0);
	case 2: return 1.05 + (opt(ratenoise) ? randf(0.05) : 0);
	case 3: return 1.10 + (opt(ratenoise) ? randf(0.1) : 0);
	case 4:	return 1.30 + (opt(ratenoise) ? randf(0.2) : 0);
	case 5:	return 1.40 + (opt(ratenoise) ? randf(0.2) : 0);
		}
	asserta(false);
	return DBL_MAX;
	}

double Peaker::GetDecreaseRateFactor(uint Rate) const
	{
	double rf = GetGlobalRateFactor();
	if (rf != DBL_MAX)
		{
		asserta(rf > 1);
		asserta(rf < 3);
		return 1.0/rf;
		}
	if (optset_rate_factor)
		return 1.0/opt(rate_factor);

	asserta(Rate >= MIN_RATE && Rate <= MAX_RATE);
	switch (Rate)
		{
	case 1:	return 0.98 - randf(0.02);
	case 2: return 0.95 - randf(0.05);
	case 3: return 0.90 - randf(0.1);
	case 4:	return 0.70 - randf(0.2);
	case 5:	return 0.60 - randf(0.2);
		}
	asserta(false);
	return DBL_MAX;
	}

double Peaker::GetGlobalRateFactor() const
	{
	string s;
	GetGlobalStr("rates", s, "");
	if (s == "")
		return DBL_MAX;
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
	GetGlobalStr("rates", s, "");
	if (s == "")
		return false;
	vector<string> Fields;
	Split(s, Fields, ',');
	const uint n = SIZE(Fields);
	asserta(m_GlobalVarRateFactorIdx < n);
	if (m_GlobalVarRateFactorIdx == n)
		return false;

	++m_GlobalVarRateFactorIdx;
	ProgressLog("\n");
	ProgressLog("Reduce global rate factor => %.3f\n", GetGlobalRateFactor());
	ProgressLog("\n");
	return true;
	}

double Peaker::GetRateFactor(uint Rate, bool Plus) const
	{
	if (Plus)
		return GetIncreaseRateFactor(Rate);
	else
		return GetDecreaseRateFactor(Rate);
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
	vector<uint> Start_Rates = m_VarRates;
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
		int dRate = int(m_VarRates[VarIdx]) - int(Start_Rates[VarIdx]);
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
	if (rf == DBL_MAX)
		{
		ProgressLog("Var rates ");
		for (uint i = 0; i < VarCount; ++i)
			ProgressLog("%u", m_VarRates[i]);
		ProgressLog("\n");
		}
	else
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

void Peaker::DecreaseRate(uint VarIdx)
	{
	asserta(VarIdx < SIZE(m_VarRates));
	uint OldRate = m_VarRates[VarIdx];
	if (OldRate <= MIN_RATE)
		return;
	m_VarRates[VarIdx] = OldRate - 1;
	}

void Peaker::IncreaseRate(uint VarIdx)
	{
	asserta(VarIdx < SIZE(m_VarRates));
	uint OldRate = m_VarRates[VarIdx];
	if (OldRate >= MAX_RATE)
		return;
	m_VarRates[VarIdx] = OldRate + 1;
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
		ProgressLog("HJ_TryDelta(%s%c) DeltaVar %s=%s no change (rate %u)\n",
			reason.c_str(), pom(Plus), VarName, OldStr.c_str(), m_VarRates[VarIdx]);
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
	Log(" [%u]", m_VarRates[VarIdx]);
	Log(" %s", NewStr.c_str());
	Log(" y %.4g,", Start_y);
	Log("%.4g", y);
	Log(" dy %.3g (target %.3g)", absdy, targetdy);

	if (OldValue != 0)
		{
		if (absdy < targetdy)
			{
			IncreaseRate(VarIdx);
			Log(" ++rate %u", m_VarRates[VarIdx]);
			}
		else if (absdy > targetdy)
			{
			DecreaseRate(VarIdx);
			Log(" --rate %u", m_VarRates[VarIdx]);
			}
		}
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
	asserta(VarIdx < SIZE(m_VarRates));
	string TmpStr;

	// Defensive to avoid infinite loop
	for (uint CheckCounter = MIN_RATE; CheckCounter < MAX_RATE;
		++CheckCounter)
		{
		uint Rate = m_VarRates[VarIdx];
		double Factor = GetRateFactor(Rate, Plus);
		double NewValue = OldValue*Factor;
		VarFloatToStr(VarIdx, NewValue, TmpStr);
		if (OldStr == TmpStr)
			{
			if (m_VarRates[VarIdx] > MIN_RATE)
				{
				--m_VarRates[VarIdx];
				continue;
				}
			else
				{
				IncFloat(OldStr, Plus, TmpStr);
				break;
				}
			}
		}
	NormalizeVarStr(VarIdx, TmpStr, NewStr);
	}
