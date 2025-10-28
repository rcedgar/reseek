#include "myutils.h"
#include "peaker.h"

double Peaker::IncreaseWithNoise() const
	{
	double Noise = rr(0, m_Noise);
	double Change = DBL_MAX;
	if (randu32()%2 == 0)
		Change = m_Change + Noise;
	else
		Change = m_Change - Noise;
	asserta(Change > 1 && Change < 10);
	return Change;
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
	LogDeltas();
	}

void Peaker::HJ_Extend()
	{
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
			{
			// Decelerate overshoot
			if (Iter > 1)
				AdjustDeltaDown(VarIdx);
			return;
			}
		// Accelerate if 2+ iters
		if (Iter > 0)
			AdjustDeltaUp(VarIdx);
		}
	}

double Peaker::DeltaVar(uint VarIdx, bool Plus, double OldValue)
	{
	uint SigFig = GetSigFig(VarIdx);
	double OldValue_Rounded = GetRounded(OldValue, SigFig);
	if (OldValue_Rounded <= 0)
		return 0;
	double Delta = m_Deltas[VarIdx];
	asserta(Delta > 0 && Delta != DBL_MAX);
	if (!feq(OldValue, OldValue_Rounded))
		Die("DeltaVar(%s) old %.4g rounded %.4g",
			OldValue, OldValue_Rounded);
	double NewValue = DBL_MAX;
	for (uint Iter = 0; Iter < 10; ++Iter)
		{
		NewValue = (Plus ? OldValue + Delta : OldValue - Delta);
		if (NewValue <= 0)
			{
			NewValue = 0;
			break;
			}
		NewValue = GetRounded(NewValue, SigFig);
		if (!feq(OldValue_Rounded, NewValue))
			break;

		Delta *= IncreaseWithNoise();
		}

	if (Delta != m_Deltas[VarIdx])
		{
		ProgressLog("Rounding increase delta %s %.3g => %3g\n",
			GetVarName(VarIdx), m_Deltas[VarIdx], Delta);
		m_Deltas[VarIdx] = Delta;
		}
	return NewValue;
	}

bool Peaker::AdjustDeltaDown(uint VarIdx)
	{
	asserta(VarIdx < SIZE(m_Deltas));
	double OldDelta = m_Deltas[VarIdx];
	double NewDelta = OldDelta/IncreaseWithNoise();
	double MinDelta = GetMinDelta(VarIdx);
	if (NewDelta < MinDelta)
		{
		Log("Adjust delta %s down %.3g undeflow\n",
			GetVarName(VarIdx), OldDelta);
		return false;
		}

	m_Deltas[VarIdx] = NewDelta;
	Log("Adjust delta %s down %.3g => %.3g\n",
		GetVarName(VarIdx), OldDelta, NewDelta);
	return true;
	}

bool Peaker::AdjustDeltaUp(uint VarIdx)
	{
	asserta(VarIdx < SIZE(m_Deltas));
	double OldDelta = m_Deltas[VarIdx];
	double NewDelta = OldDelta*IncreaseWithNoise();
	Log("Adjust delta %s up %.3g => %.3g\n",
		GetVarName(VarIdx), OldDelta, NewDelta);
	m_Deltas[VarIdx] = NewDelta;
	return true;
	}

double Peaker::HJ_TryDelta(const string &reason,
	const vector<double> &Start_xv, uint VarIdx, bool Plus)
	{
	const char *VarName = GetVarName(VarIdx);
	uint Idx = Find_xs(Start_xv);
	asserta(Idx != UINT_MAX);
	asserta(Idx < SIZE(m_ys));
	const double Start_y = m_ys[Idx];

	double Delta = m_Deltas[VarIdx];
	asserta(Delta != DBL_MAX);
	asserta(!isnan(Delta));

	vector<double> Try_xv = Start_xv;
	asserta(VarIdx < SIZE(Try_xv));

	double NewValue = DeltaVar(VarIdx, Plus, Start_xv[VarIdx]);
	if (feq(NewValue, Try_xv[VarIdx]))
		{
		ProgressLog("HJ_TryDelta(%s) DeltaVar %s no change\n",
			reason.c_str(), VarName);
		return Start_y;
		}

	Try_xv[VarIdx] = NewValue;

	string why;
	Ps(why, "%s%c%s", reason.c_str(), pom(Plus), VarName);
	double y = Evaluate(Try_xv, why);
	if (y == DBL_MAX)
		return DBL_MAX;
	if (Start_y == DBL_MAX)
		return y;

	double dy = fabs(y - Start_y);
	if (dy < m_Min_dy)
		AdjustDeltaUp(VarIdx);
	else if (dy > m_Max_dy)
		AdjustDeltaDown(VarIdx);
	return y;
	}

void Peaker::HJ_RunHookeJeeves()
	{
	InitDeltas();
	asserta(m_Min_Height != DBL_MAX);
	for (uint Iters = 0; ; ++Iters)
		{
		if (Iters >= m_HJ_MaxIters)
			{
			Warning("HJ max iters, not converged");
			break;
			}
		double Saved_Best_y = m_Best_y;
		HJ_Explore();
		double Height = m_Best_y - Saved_Best_y;
		asserta(Height >= 0);
		ProgressLog("HJ height=%.3g (minh %.3g)\n",
			Height, m_Min_Height);
		LogDeltas();
		if (Height <= m_Min_Height)
			{
			ProgressLog("HJ converged by height\n");
			break;
			}
		if (m_HJ_Direction == UINT_MAX)
			{
			ProgressLog("HJ converged by no improvement found\n");
			break;
			}
		HJ_Extend();
		}
	}
