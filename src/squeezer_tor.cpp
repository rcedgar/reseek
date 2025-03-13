#include "myutils.h"
#include "squeezer.h"

SqueezeDelta *Squeezer::EncodePos_algo_tor(
	const SqueezeState &State, const coords &TrueD) const
	{
	SqueezeDelta *Delta = new SqueezeDelta;
	if (State.m_i < 3)
		{
		Delta->m_op = SDOP_set_coords;
		Delta->m_x = TrueD.x;
		Delta->m_y = TrueD.y;
		Delta->m_z = TrueD.z;
		return Delta;
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
		Delta->m_op = SDOP_tor;
		Delta->m_theta_rad = theta_rad;
		Delta->m_phi_rad = phi_rad;
		}
	else
		{
		Delta->m_op = SDOP_set_coords;
		Delta->m_x = TrueD.x;
		Delta->m_y = TrueD.y;
		Delta->m_z = TrueD.z;
		}

	return Delta;
	}
