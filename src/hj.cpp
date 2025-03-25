#include "myutils.h"
#include "peaker.h"

double Peaker::ChangeWithNoise() const
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

// Returns height = largest |dy|
double Peaker::HJ_Explore()
	{
	const uint VarCount = GetVarCount();

	uint BestNewDirection = UINT_MAX;
	m_HJ_Plus = false;
	double Saved_Best_y = m_Best_y;
	double BestNew_y = m_Best_y;
	uint Saved_Best_xIdx = m_Best_xIdx;
	double maxdy = 0;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const char *Name = GetVarName(VarIdx);
		for (int iPlus = 0; iPlus <= 1; ++iPlus)
			{
			bool Plus = (iPlus == 1);
			if (VarIdx == m_HJ_Direction)
				continue;
			Ps(m_Msg, "hj_%s%c/%u/%u",
			  Name, pom(Plus), VarIdx+1, VarCount);
			double y = HJ_EvalDelta(Saved_Best_xIdx, VarIdx, Plus);
			maxdy = max(maxdy, fabs(y - Saved_Best_y));
			if (y > BestNew_y)
				{
				ProgressLog("HJ_Explore() improved %s%c dy=%.3g\n",
				  Name, pom(Plus), y - Saved_Best_y);
				BestNew_y = y;
				m_HJ_Plus = Plus;
				BestNewDirection = VarIdx;
				}
			}
		}
	m_HJ_Direction = BestNewDirection;
	if (m_HJ_Direction == UINT_MAX)
		ProgressLog("HJ_Expore(), no improvement found\n");
	else
		ProgressLog("HJ_Expore(), new direction %s%c\n",
		  GetVarName(m_HJ_Direction), pom(m_HJ_Plus));
	LogDeltas();
	return maxdy;
	}

void Peaker::HJ_Extend()
	{
	const uint VarCount = GetVarCount();
	asserta(m_HJ_Direction < VarCount);
	const uint VarIdx = m_HJ_Direction;
	const char *Name = GetVarName(VarIdx);
	for (uint Iter = 0; ; ++Iter)
		{
		if (Iter > m_HJ_MaxExtendIters)
			return;
		
		double Delta = m_Deltas[VarIdx];
		double Saved_Best_y = m_Best_y;
		Ps(m_Msg, "HJ_%s%c", GetVarName(VarIdx), pom(m_HJ_Plus));
		double y = HJ_EvalDelta(m_Best_xIdx, VarIdx, Delta);
		ProgressLog("HJ_extend() %s%c %u/%u dy=%+.3g\n",
		  Name, pom(m_HJ_Plus), Iter+1, m_HJ_MaxExtendIters, y - Saved_Best_y);
		if (y <= Saved_Best_y)
			{
			double dy = Saved_Best_y - y;
			asserta(dy >= 0);
			if (Iter == 0 && dy > m_Min_dy)
				{
				ProgressLog("HJ_Extend(%s) delta down\n", Name);
				m_Deltas[VarIdx] /= ChangeWithNoise();
				LogDeltas();
				}
			return;
			}
		if (Iter > 0)
			{
			ProgressLog("HJ_Extend(%s) delta up\n", Name, Iter+1);
			m_Deltas[VarIdx] *= ChangeWithNoise();
			LogDeltas();
			}
		}
	}

double Peaker::HJ_EvalDelta(uint xIdx, uint VarIdx, bool Plus)
	{
	double Delta = m_Deltas[VarIdx];
	vector<double> xv = Get_xv(xIdx);

	const char *Name = GetVarName(VarIdx);
	bool Vanished = false;
	double Oldx = xv[VarIdx];
	double Newx = DBL_MAX;
	uint SigFig = GetSigFig(VarIdx);
	double Oldx_Rounded = GetRounded(Oldx, SigFig);
	for (uint Iter = 0; ; ++Iter)
		{
		if (Iter > 100)
			Die("HJ_EvalDelta rounding vanished");

		Newx = Plus ? Oldx + Delta : Oldx - Delta;

		double Newx_Rounded = GetRounded(Newx, SigFig);
		if (Oldx_Rounded != Newx_Rounded)
			{
			xv[VarIdx] = Newx;
			break;
			}

		Warning("Delta %s vanished, increasing", Name);
		m_Deltas[VarIdx] *= ChangeWithNoise();
		Vanished = true;
		}

	uint xIdx2 = Add_xv(xv);
	double y = Evaluate(xIdx);
	double y2 = Evaluate(xIdx2);
	double dy = fabs(y2 - y);
	Log("HJ_EvalDelta %s%c %.3g => %.3g, y %.6g => %.6g (%+.3g)\n",
	  Name, pom(Plus), Oldx, Newx, y, y2, y2 - y);
	if (dy < m_Min_dy)
		{
		double NewDelta = Delta*ChangeWithNoise();
		m_Deltas[VarIdx] = NewDelta;
		ProgressLog("Delta %s adjust up %.3g => %.3g\n",
		  Name, Delta, NewDelta);
		}
	else if (dy > m_Max_dy)
		{
		double NewDelta = Delta/ChangeWithNoise();
		m_Deltas[VarIdx] = NewDelta;
		ProgressLog("Delta %s adjust down %.3g => %.3g\n",
		  Name, Delta, NewDelta);
		}
	return y2;
	}

void Peaker::HJ_RunHookeJeeves()
	{
	InitDeltas();
	asserta(m_Best_xIdx != UINT_MAX);
	asserta(m_Min_Height != DBL_MAX);
	for (uint Iters = 0; ; ++Iters)
		{
		if (Iters >= m_HJ_MaxIters)
			{
			Warning("HJ max iters, not converged");
			break;
			}
		double Height = HJ_Explore();
		ProgressLog("HJ height=%.3g\n", Height);
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
	ProgressLogSummary();
	}
