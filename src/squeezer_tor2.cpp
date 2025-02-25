#include "myutils.h"
#include "squeezer.h"

SqueezeDelta *Squeezer::EncodePos_algo_tor2(
	const SqueezeState &State, const coords &TrueD) const
	{
	SqueezeDelta *Delta = new SqueezeDelta;
	if (State.m_i < 3)
		{
		Delta->m_op = SDOP_IC;
		Delta->m_ICx = PDBChain::CoordToIC(TrueD.x);
		Delta->m_ICy = PDBChain::CoordToIC(TrueD.y);
		Delta->m_ICz = PDBChain::CoordToIC(TrueD.z);
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

	int8_t itheta = theta2i(theta_rad);
	int8_t iphi = phi2i(phi_rad);
	float theta_rad2 = i2theta(itheta);
	float phi_rad2 = i2phi(iphi);

	coords D2;
	calculate_D(
		State.m_A,
		State.m_B,
		State.m_C,
		theta_rad2,
		phi_rad2,
		D2);

	float e = coords::dist(TrueD, D2);
	if (e <= m_MaxErr)
		{
		Delta->m_op = SDOP_tor2;
		Delta->m_itheta = itheta;
		Delta->m_iphi = iphi;
#if DEBUG
		{
		coords TestD2;
		DecodePos_tor2(State, *Delta, TestD2);
		float e2 = coords::dist(TrueD, TestD2);
		asserta(e2 <= m_MaxErr);
		}
#endif
		}
	else
		{
		int16_t icx = PDBChain::CoordToIC(State.m_C.x);
		int16_t icy = PDBChain::CoordToIC(State.m_C.y);
		int16_t icz = PDBChain::CoordToIC(State.m_C.z);

		int16_t itruedx = PDBChain::CoordToIC(TrueD.x);
		int16_t itruedy = PDBChain::CoordToIC(TrueD.y);
		int16_t itruedz = PDBChain::CoordToIC(TrueD.z);

		int8_t fix_x = itruedx - icx;
		int8_t fix_y = itruedy - icy;
		int8_t fix_z = itruedz - icz;

		int16_t idx = icx + fix_x;
		int16_t idy = icy + fix_y;
		int16_t idz = icz + fix_z;

		coords TestD;
		TestD.x = PDBChain::ICToCoord(idx);
		TestD.y = PDBChain::ICToCoord(idy);
		TestD.z = PDBChain::ICToCoord(idz);

		float teste = coords::dist(TrueD, TestD);
		if (teste <= m_MaxErr)
			{
			Delta->m_ICdx = fix_x;
			Delta->m_ICdy = fix_y;
			Delta->m_ICdz = fix_z;
			Delta->m_op = SDOP_ICd;
#if DEBUG
			{
			coords TestD2;
			DecodePos_ICd(State, *Delta, TestD2);
			float e2 = coords::dist(TrueD, TestD2);
			asserta(e2 <= m_MaxErr);
			}
#endif
			}
		else
			{
			Delta->m_ICx = PDBChain::CoordToIC(TrueD.x);
			Delta->m_ICy = PDBChain::CoordToIC(TrueD.y);
			Delta->m_ICz = PDBChain::CoordToIC(TrueD.z);
			Delta->m_op = SDOP_IC;
			}
		}

#if DEBUG
	coords D3;
	DecodePos(State, *Delta, D3);
	float e3 = coords::dist(TrueD, D3);
	asserta(e3 <= m_MaxErr);
#endif
	return Delta;
	}
