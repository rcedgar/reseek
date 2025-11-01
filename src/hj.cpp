#include "myutils.h"
#include "peaker.h"

double Peaker::GetIncreaseRateFactor(uint Rate)
	{
	asserta(Rate >= MIN_RATE && Rate <= MAX_RATE);
	switch (Rate)
		{
	case 1:	return 1.05 + randf(0.05);
	case 2: return 1.10 + randf(0.1);
	case 3: return 1.30 + randf(0.2);
	case 4:	return 1.50 + randf(0.3);
	case 5:	return 2.00 + randf(0.4);
		}
	asserta(false);
	return DBL_MAX;
	}

double Peaker::GetDecreaseRateFactor(uint Rate)
	{
	asserta(Rate >= MIN_RATE && Rate <= MAX_RATE);
	switch (Rate)
		{
	case 1:	return 0.95 - randf(0.05);
	case 2: return 0.90 - randf(0.05);
	case 3: return 0.80 - randf(0.1);
	case 4:	return 0.60 - randf(0.1);
	case 5:	return 0.50 - randf(0.1);
		}
	asserta(false);
	return DBL_MAX;
	}

double Peaker::GetRateFactor(uint Rate, bool Plus)
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
	double Saved_Best_y = m_Best_y;
	double BestNew_y = m_Best_y;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		if (VarIsConstant(VarIdx))
			continue;
		for (int iPlus = 0; iPlus <= 1; ++iPlus)
			{
			bool Plus = (iPlus == 1);
			if (VarIdx == m_HJ_Direction)
				continue;
			double y = HJ_TryDelta("explore", m_Best_xv, VarIdx, Plus);
			if (y == DBL_MAX)
				continue;
			if (y > BestNew_y)
				{
				BestNew_y = y;
				m_HJ_ExtendPlus = Plus;
				BestNewDirection = VarIdx;
				}
			}
		}
	m_HJ_Direction = BestNewDirection;
	if (m_HJ_Direction == UINT_MAX)
		Log("HJ_Expore(), no improvement found\n");
	else
		Log("HJ_Expore(), new direction %s%c\n",
		  GetVarName(m_HJ_Direction), pom(m_HJ_ExtendPlus));
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
	for (uint Iter = 0; Iter < m_HJ_MaxExtendIters; ++Iter)
		{
		string reason;
		Ps(reason, "extend%u", Iter+1);

		double Saved_Best_y = m_Best_y;
		HJ_TryDelta(reason, m_Best_xv, VarIdx, m_HJ_ExtendPlus);
		if (m_Best_y <= Saved_Best_y)
			return;
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
	const vector<string> &Start_xv, uint VarIdx, bool Plus)
	{
	const char *VarName = GetVarName(VarIdx);
	uint Idx = Find_xv(Start_xv);
	asserta(Idx != UINT_MAX);
	asserta(Idx < SIZE(m_ys));
	const double Start_y = m_ys[Idx];

	string NewStr;
	const string &OldStr = Start_xv[VarIdx];
	DeltaVar(VarIdx, Plus, OldStr, NewStr);
	if (NewStr == Start_xv[VarIdx])
		{
		ProgressLog("HJ_TryDelta(%s) DeltaVar %s=%s no change (rate %u)\n",
			reason.c_str(), VarName, OldStr.c_str(), m_VarRates[VarIdx]);
		return Start_y;
		}

	vector<string> Try_xv = Start_xv;
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

	Log("HJ_TryDelta(%s)", reason.c_str());
	Log(" %s", VarName);
	Log(" %s,", OldStr.c_str());
	Log(" %s", NewStr.c_str());
	Log(" y %.4g,", Start_y);
	Log("%.4g", y);
	Log(" dy %.3g (target %.3g)", absdy, targetdy);

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
		HJ_Extend();
		return true;
		}

	asserta(Height <= 0);
	const uint VarCount = GetVarCount();
	bool Any = false;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		uint Rate = m_VarRates[VarIdx];
		if (Rate > MIN_RATE)
			{
			m_VarRates[VarIdx] = MIN_RATE;
			Any = true;
			}
		}
	ProgressLog("HJ converged, no improvement found\n");
	return false;
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
		double MinValue = VarSpecGetFloat(VarIdx, "min", DBL_MAX);
		string TmpStr;
		VarFloatToStr(VarIdx, MinValue, TmpStr);
		NormalizeVarStr(VarIdx, TmpStr, NewStr);
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
