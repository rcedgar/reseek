#include "myutils.h"
#include "peaker.h"
#include "sort.h"

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

bool Peaker::CanReduceGlobalRateFactor()
	{
	string s;
	GetGlobalStr("rates", s, "1.1");
	vector<string> Fields;
	Split(s, Fields, ',');
	const uint n = SIZE(Fields);
	asserta(m_GlobalVarRateFactorIdx < n);
	if (m_GlobalVarRateFactorIdx + 1 == n)
		return false;
	return true;
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
		ProgressLogNoPrefix("%s: ReduceGlobalRateFactor() no more, elapsed %s [%.5g]\n",
			m_Name.c_str(), GetElapsedTimeStr(s), m_Best_y);
		return false;
		}

	++m_GlobalVarRateFactorIdx;
	ProgressLogNoPrefix("\n");
	ProgressLogNoPrefix("%s: ReduceGlobalRateFactor() => %.3f, elapsed %s [%.5g]\n",
		m_Name.c_str(), GetGlobalRateFactor(), GetElapsedTimeStr(s), m_Best_y);
	ProgressLogNoPrefix("\n");
	return true;
	}

double Peaker::GetRateFactor(bool Plus)
	{
	if (Plus)
		return GetIncreaseRateFactor();
	else
		return GetDecreaseRateFactor();
	}

bool Peaker::VarIsStalled(uint VarIdx) const
	{
	asserta(VarIdx < SIZE(varidx2last_improved_hjiter));
	uint n = varidx2last_improved_hjiter[VarIdx];
	return m_HJIter - n >= 3;
	}

bool Peaker::AnyStalledVars() const
	{
	const uint VarCount = GetVarCount();
	for (uint i = 0; i < VarCount; ++i)
		if (VarIsStalled(i)) return true;
	return false;
	}

void Peaker::HJ_Explore(bool try_stalled, bool stalled_only)
	{
	const uint VarCount = GetVarCount();
	if (stalled_only) try_stalled = true;

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
	vector<uint> Order;
	Range(Order, VarCount);
	if (GetGlobalBool("shuffle", true))
		{
		ProgressLog("-- shuffle vars --\n");
		Shuffle(Order);
		}
	uint LastVarIdx = UINT_MAX;
	for (uint k = 0; k < VarCount; ++k)
		{
		uint VarIdx = Order[k];
		if (VarIdx == LastVarIdx)
			continue;
		if (VarIsConstant(VarIdx))
			continue;
		bool stalled = VarIsStalled(VarIdx);
		if (stalled)
			{
			ProgressLog("{{ %s stalled", GetVarName(VarIdx));
			if (try_stalled)
				ProgressLog(" ... try }}\n");
			else
				ProgressLog(" ... skip }}\n");
			continue;
			}
		if (!stalled && stalled_only)
			continue;
		LastVarIdx = VarIdx;
		double Saved_Best_y = m_Best_y;
		string reason;
		Ps(reason, "ex%u/%u", k+1, VarCount);
		double y_plus = HJ_TryDelta(reason, m_Best_xv, VarIdx, true, Try_xv);
		double dy_plus = y_plus - Saved_Best_y;
		dys_plus[VarIdx] = dy_plus;
		strs_plus[VarIdx] = Try_xv[VarIdx];
		if (dy_plus > 0)
			{
			varidx2last_improved_hjiter[VarIdx] = m_HJIter;
			++ImprovementCount;
			}
		if (dy_plus > Best_dy)
			{
			Best_dy = dy_plus;
			m_HJ_ExtendPlus = true;
			BestNewDirection = VarIdx;
			continue;
			}

		Saved_Best_y = m_Best_y;
		Ps(reason, "ex%u/%u", k+1, VarCount);
		double y_minus = HJ_TryDelta(reason, m_Best_xv, VarIdx, false, Try_xv);
		double dy_minus = y_minus - Saved_Best_y;
		dys_minus[VarIdx] = dy_minus;
		strs_minus[VarIdx] = Try_xv[VarIdx];
		if (dy_minus > 0)
			++ImprovementCount;
		if (dy_minus > Best_dy)
			{
			varidx2last_improved_hjiter[VarIdx] = m_HJIter;
			Best_dy = dy_minus;
			m_HJ_ExtendPlus = false;
			BestNewDirection = VarIdx;
			}
		}
	m_HJ_Direction = BestNewDirection;
	if (m_HJ_Direction == UINT_MAX)
		{
		Log("%s: HJ_Explore(), no improvement found\n", m_Name.c_str());
		return;
		}
	Log("%s: HJ_Explore(), new direction %s%c\n", m_Name.c_str(),
		GetVarName(m_HJ_Direction), pom(m_HJ_ExtendPlus));

	ProgressLogNoPrefix("\n");
	ProgressLogNoPrefix("%s: HJ_Explore(iter=%u) %u improves\n",
		m_Name.c_str(), m_HJIter+1, ImprovementCount);
	double Track_Best_y = Start_Best_y;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		double dy_plus = dys_plus[VarIdx];
		double dy_minus = dys_minus[VarIdx];
		double pct_plus = GetPct(abs(dy_plus), m_Best_y);
		double pct_minus = GetPct(abs(dy_minus), m_Best_y);
		char sign_plus = pom(dy_plus >= 0);
		char sign_minus = pom(dy_minus >= 0);
		bool VarIsStalled(VarIdx);
		uint di = varidx2last_improved_hjiter[VarIdx] - m_HJIter;

		ProgressLogNoPrefix("%c%5.2f%%", sign_plus, pct_plus);
		ProgressLogNoPrefix("  <%-10.10s", strs_minus[VarIdx].c_str());
		ProgressLogNoPrefix("  %c%5.2f%%", sign_minus, pct_minus);
		ProgressLogNoPrefix("  %s", GetVarName(VarIdx));
		if (dys_plus[VarIdx] > 0 || dys_minus[VarIdx] > 0)
			ProgressLogNoPrefix(" +++ improved");
		else
			{
			ProgressLogNoPrefix("  di=%u", di);
			if (VarIsStalled)
				ProgressLogNoPrefix(" STALLED");
			}
		ProgressLogNoPrefix("\n");
		}
	ProgressLogNoPrefix("\n");

	const uint N = SIZE(m_Best_ys);
	asserta(SIZE(m_Best_descs) == N);
	for (uint k = 0; k < 5; ++k)
		{
		if (k > N -1)
			break;
		uint i = N-k-1;
		double dy = (i > 0 ? m_Best_ys[i] - m_Best_ys[i-1] : 0);
		ProgressLogNoPrefix("%10.5g", m_Best_ys[i]);
		ProgressLogNoPrefix("  %+10.2g", dy);
		ProgressLogNoPrefix("  %s", m_Best_descs[i].c_str());
		ProgressLogNoPrefix("\n");
		}
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
		ProgressLogNoPrefix("%s: HJ_TryDelta(%s%c) DeltaVar %s=%s no change\n",
			m_Name.c_str(), reason.c_str(), pom(Plus), VarName, OldStr.c_str());
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

	Log("HJ_TryDelta(%s%c)", reason.c_str(), pom(Plus));
	Log(" %s", VarName);
	Log(" %s,", OldStr.c_str());
	Log(" %s", NewStr.c_str());
	Log(" y %.4g,", Start_y);
	Log("%.4g", y);
	Log(" dy %.3g", absdy);
	Log("\n");

	return y;
	}

bool Peaker::HJ_Iter()
	{
	double Saved_Best_y = m_Best_y;
	bool try_stalled = (m_HJIter%4 == 0);
	HJ_Explore(true, false);
	double Height = m_Best_y - Saved_Best_y;
	if (Height == 0 && !try_stalled)
		{
		HJ_Explore(false, true);
		Height = m_Best_y - Saved_Best_y;
		}

	double Pct = 0;
	if (m_Best_y != 0 && Height > 0 && Height < m_Best_y)
		Pct = GetPct(Height, m_Best_y);
	ProgressLogNoPrefix("%s: [%.5g] HJ_Iter() height %.3g (%.2f%%) /%.2f/\n\n",
		m_Name.c_str(), m_Best_y, Height, Pct, GetGlobalRateFactor());
	if (Height > 0)
		{
		bool CanReduce = CanReduceGlobalRateFactor();
		if (Pct < m_ConvergePct && !CanReduce)
			return false;
		if (Pct < m_ConvergeReducePct && CanReduce)
			{
			ProgressLogNoPrefix("%s: HJ_Iter() small (<0.1) height pct %.3g%%\n",
				m_Name.c_str(), Pct);
			bool ok = ReduceGlobalRateFactor();
			asserta(ok);
			return true;
			}
		if (GetGlobalBool("extend", false))
			HJ_Extend();
		return true;
		}
	bool ok = ReduceGlobalRateFactor();
	ProgressLogNoPrefix("%s: HJ_Iter() ReduceGlobalRateFactor %c\n", m_Name.c_str(), tof(ok));
	return ok;
	}

void Peaker::HJ_RunHookeJeeves()
	{
	InitRates();
	varidx2last_improved_hjiter.clear();
	varidx2last_improved_hjiter.resize(GetVarCount());
	for (m_HJIter = 0; ; ++m_HJIter)
		{
		ProgressLog("\n\n============ HJ iter %u ============\n\n", m_HJIter+1);
		if (m_HJIter >= m_HJ_MaxIters)
			{
			Warning("HJ max iters, not converged");
			break;
			}
		bool ok = HJ_Iter();
		if (!ok)
			return;
		}
	}

void Peaker::GetInitialVarStr(string &Str) const
	{
	Str.clear();
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

void Peaker::DeltaVarInt(uint VarIdx, bool Plus, const string &OldStr,
	string &NewStr)
	{
	double OldValue = VarStrToFloat(VarIdx, OldStr);
	int iOldValue = int(round(OldValue));
	asserta(float(iOldValue) == OldValue);
	int iNewValue = (Plus ? iOldValue + 1 : iOldValue - 1);
	if (iNewValue < 0)
		iNewValue = 0;
	VarFloatToStr(VarIdx, float(iNewValue), NewStr);
	}

void Peaker::DeltaVar(uint VarIdx, bool Plus,
	const string &OldStr, string &NewStr)
	{
	if (VarIsInt(VarIdx))
		{
		DeltaVarInt(VarIdx, Plus, OldStr, NewStr);
		return;
		}

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
		ProgressLogNoPrefix("%s ifzero %s => %s\n", GetVarName(VarIdx),
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
		IncFloat(OldStr, Plus, TmpStr, SigFig);
		NormalizeVarStr(VarIdx, TmpStr, NewStr);
		}
	}
