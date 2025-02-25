#include "myutils.h"
#include "squeezer_tor.h"

// theta 25 .. 100 deg.
// phi -180 .. 180 deg.

void Squeezer_tor::DecodePos(const SqueezeState &State,
						const SqueezeDelta &Delta,
						coords &D) const
	{
	const SqueezeDelta_tor &Delta_tor =
		(const SqueezeDelta_tor &) Delta;
	DecodePos_tor(State, Delta_tor, D);
	}

void Squeezer_tor::DecodePos_tor(const SqueezeState &State_tor,
						const SqueezeDelta_tor &Delta_tor,
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
		D.x = Delta_tor.m_x;
		D.y = Delta_tor.m_y;
		D.z = Delta_tor.m_z;
		}

	else
		asserta(false);
	}

SqueezeDelta *Squeezer_tor::EncodePos(
	const SqueezeState &State, uint i) const
	{
	SqueezeDelta_tor *Delta_tor = new SqueezeDelta_tor;
	coords TrueD;
	m_Chain->GetCoords(i, TrueD);
	if (i < 3)
		{
		Delta_tor->m_tord = TORD_set_coords;
		Delta_tor->m_x = TrueD.x;
		Delta_tor->m_y = TrueD.y;
		Delta_tor->m_z = TrueD.z;
		return Delta_tor;
		}

	float theta_rad, phi_rad;
	State.m_A.assert_valid();
	State.m_B.assert_valid();
	State.m_C.assert_valid();
	calculate_theta_phi(
		State.m_A,
		State.m_B,
		State.m_C,
		TrueD,
		theta_rad,
		phi_rad);

	coords D2;
	calculate_D(
		State.m_A,
		State.m_B,
		State.m_C,
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
		Delta_tor->m_x = TrueD.x;
		Delta_tor->m_y = TrueD.y;
		Delta_tor->m_z = TrueD.z;
		}

	return Delta_tor;
	}
