#include "myutils.h"
#include "squeezer_tor.h"

void Squeezer_tor::InitDecode(SqueezeState *State) const
	{
	SqueezeState_tor *State_tor = (SqueezeState_tor *) State;
	State_tor->m_A.invalidate();
	State_tor->m_B.invalidate();
	State_tor->m_C.invalidate();
	State_tor->m_ResetCount = 0;
	}

void Squeezer_tor::DecodePos(const SqueezeState &State,
						const SqueezeDelta &Delta,
						uint i,
						coords &D) const
	{
	SqueezeState_tor &State_tor = (SqueezeState_tor &) State;
	const SqueezeDelta_tor &Delta_tor =
		(const SqueezeDelta_tor &) Delta;
	DecodePos_tor(State_tor, Delta_tor, i, D);
	}

void Squeezer_tor::DecodePos_tor(const SqueezeState_tor &State_tor,
						const SqueezeDelta_tor &Delta_tor,
						uint i,
						coords &D) const
	{
	if (Delta_tor.m_tord == TORD_tor)
		{
		State_tor.m_A.assert_valid();
		State_tor.m_B.assert_valid();
		State_tor.m_C.assert_valid();
		calculate_D(State_tor.m_A,
					State_tor.m_B,
					State_tor. m_C,
					Delta_tor.m_theta_rad,
					Delta_tor.m_phi_rad,
					D);
		}

	else if (Delta_tor.m_tord == TORD_set_coords)
		{
		D.x = Delta_tor.m_dx;
		D.y = Delta_tor.m_dy;
		D.z = Delta_tor.m_dz;
		}

	else if (Delta_tor.m_tord == TORD_delta_coords)
		{
		asserta(false);
		}
	else
		asserta(false);
	}

SqueezeDelta *Squeezer_tor::EncodePos(const SqueezeState &State, uint i) const
	{
	const SqueezeState_tor &State_tor = (const SqueezeState_tor &) State;
	SqueezeState_tor BeforeState;
	BeforeState = State_tor;

	SqueezeDelta_tor *Delta_tor = EncodePos_tor(State_tor, i);

#if DEBUG
	coords D;
	DecodePos_tor(State_tor, *Delta_tor, i, D);
	coords TrueD;
	GetTrueD(i, TrueD);
	float e = GetErr(D, TrueD);
	Log("[%3u]", i);
	Log("  %.1f (%.1f)", D.x, TrueD.x);
	Log(",  %.1f (%.1f)", D.y, TrueD.y);
	Log(",  %.1f (%.1f)", D.z, TrueD.z);
	Log(" e=%.3g", e);
	Log(" resets=%u", State_tor.m_ResetCount);
	Log("\n");
#endif
	return Delta_tor;
	}

SqueezeDelta_tor *Squeezer_tor::EncodePos_tor(
	const SqueezeState_tor &State_tor, uint i) const
	{
	SqueezeDelta_tor *Delta_tor = new SqueezeDelta_tor;
	coords TrueD;
	m_Chain->GetCoords(i, TrueD);
	if (i < 3)
		{
		Delta_tor->m_tord = TORD_set_coords;
		Delta_tor->m_dx = TrueD.x;
		Delta_tor->m_dy = TrueD.y;
		Delta_tor->m_dz = TrueD.z;
		return Delta_tor;
		}

	float theta_rad, phi_rad;
	State_tor.m_A.assert_valid();
	State_tor.m_B.assert_valid();
	State_tor.m_C.assert_valid();
	calculate_theta_phi(
		State_tor.m_A,
		State_tor.m_B,
		State_tor.m_C,
		TrueD,
		theta_rad,
		phi_rad);

	coords D2;
	calculate_D(
		State_tor.m_A,
		State_tor.m_B,
		State_tor.m_C,
		theta_rad,
		phi_rad,
		D2);

	float e = coords::dist(TrueD, D2);
	if (e <= m_MaxErr)
		{
		Delta_tor->m_tord = TORD_tor;
		Delta_tor->m_theta_rad = theta_rad;
		Delta_tor->m_phi_rad = phi_rad;
		}
	else
		{
		Delta_tor->m_tord = TORD_set_coords;
		Delta_tor->m_dx = TrueD.x;
		Delta_tor->m_dy = TrueD.y;
		Delta_tor->m_dz = TrueD.z;
		}

	return Delta_tor;
	}

void Squeezer_tor::UpdateState(SqueezeState &State,
					   const SqueezeDelta &Delta) const
	{
	SqueezeState_tor &State_tor = (SqueezeState_tor &) State;
	const SqueezeDelta_tor &Delta_tor = (SqueezeDelta_tor &) Delta;
	UpdateState_tor(State_tor, Delta_tor);
	}

void Squeezer_tor::UpdateState_tor(SqueezeState_tor &State_tor,
					   const SqueezeDelta_tor &Delta_tor) const
	{

	if (Delta_tor.m_tord == TORD_tor)
		{
		coords D;
		calculate_D(
			State_tor.m_A,
			State_tor.m_B,
			State_tor.m_C,
			Delta_tor.m_theta_rad,
			Delta_tor.m_phi_rad,
			D);

		State_tor.m_A = State_tor.m_B;
		State_tor.m_B = State_tor.m_C;
		State_tor.m_C = D;
		}
	else if (Delta_tor.m_tord == TORD_set_coords)
		{
		State_tor.m_A = State_tor.m_B;
		State_tor.m_B = State_tor.m_C;

		State_tor.m_C.x = Delta_tor.m_dx;
		State_tor.m_C.y = Delta_tor.m_dy;
		State_tor.m_C.z = Delta_tor.m_dz;

		State_tor.m_ResetCount++;
		}
	}
