#include "myutils.h"
#include "squeezer_tor2.h"
#include "mypi.h"

// theta	 25 .. 100 deg.		0.4363323 .. 1.745329 radians.
// phi	   -180 .. 180 deg.		-PI .. +PI

static const float MIN_THETA = 0.436f;
static const float MAX_THETA = 1.746f;

static const float MIN_PHI = -PI;
static const float MAX_PHI = PI;

void Squeezer_tor2::test_cvt()
	{
	for (int itheta = 25; itheta <= 100; ++itheta)
		{
		float theta_deg = float(itheta);
		float theta_rad = radians(theta_deg);
		int8_t i2 = theta2i(theta_rad);
		float theta_rad2 = i2theta(i2);
		float theta_deg2 = degrees(theta_rad2);
		Log("theta %4.1f  %4.1f\n", theta_deg, theta_deg2);
		}

	for (int iphi = -180; iphi <= 180; ++iphi)
		{
		float phi_deg = float(iphi);
		float phi_rad = radians(phi_deg);
		int8_t i2 = phi2i(phi_rad);
		float phi_rad2 = i2phi(i2);
		float phi_deg2 = degrees(phi_rad2);
		Log("phi %4.1f  %4.1f\n", phi_deg, phi_deg2);
		}
	}

int8_t Squeezer_tor2::theta2i(float theta)
	{
	float zero_to_one = (theta - MIN_THETA)/(MAX_THETA - MIN_THETA);
	float minus_one_to_plus_one = 2.0f*(zero_to_one - 0.5f);
	float minus_127_to_plus_127 = 127*minus_one_to_plus_one;
	int8_t i = int8_t(minus_127_to_plus_127);
	return i;
	}

int8_t Squeezer_tor2::phi2i(float phi)
	{
	float zero_to_one = (phi - MAX_PHI)/(MAX_PHI - MIN_PHI);
	float minus_one_to_plus_one = 2.0f*(zero_to_one - 0.5f);
	float minus_127_to_plus_127 = 127*minus_one_to_plus_one;
	int8_t i = int8_t(minus_127_to_plus_127);
	return i;
	}

float Squeezer_tor2::i2theta(int8_t itheta)
	{
	float minus_127_to_plus_127 = float(itheta);
	float minus_one_to_plus_one = minus_127_to_plus_127/127.0f;
	float zero_to_one = 0.5f + minus_one_to_plus_one/2.0f;
	float theta = MIN_THETA + zero_to_one*(MAX_THETA - MIN_THETA);
	return theta;
	}

float Squeezer_tor2::i2phi(int8_t iphi)
	{
	float minus_127_to_plus_127 = float(iphi);
	float minus_one_to_plus_one = minus_127_to_plus_127/127.0f;
	float zero_to_one = 0.5f + minus_one_to_plus_one/2.0f;
	float phi = MIN_PHI + zero_to_one*(MAX_PHI - MIN_PHI);
	return phi;
	}

void Squeezer_tor2::DecodePos(const SqueezeState &State,
						const SqueezeDelta &Delta,
						coords &D) const
	{
	const SqueezeDelta_tor2 &Delta_tor2 = (const SqueezeDelta_tor2 &) Delta;
	DecodePos_tor2(State, Delta_tor2, D);
	}

void Squeezer_tor2::DecodePos_tor2(const SqueezeState &State_tor2,
						const SqueezeDelta_tor2 &Delta_tor2,
						coords &D) const
	{
	if (Delta_tor2.m_tor2d == TOR2D_tor)
		{
		State_tor2.m_A.assert_valid();
		State_tor2.m_B.assert_valid();
		State_tor2.m_C.assert_valid();

		float theta_rad = i2theta(Delta_tor2.m_itheta);
		float phi_rad = i2phi(Delta_tor2.m_iphi);
		calculate_D(State_tor2.m_A,
					State_tor2.m_B,
					State_tor2. m_C,
					theta_rad,
					phi_rad,
					D);
		}

	else if (Delta_tor2.m_tor2d == TOR2D_fix_coords)
		{
		int16_t icx = PDBChain::CoordToIC(State_tor2.m_C.x) + Delta_tor2.m_fix_x;
		int16_t icy = PDBChain::CoordToIC(State_tor2.m_C.y) + Delta_tor2.m_fix_y;
		int16_t icz = PDBChain::CoordToIC(State_tor2.m_C.z) + Delta_tor2.m_fix_z;

		D.x = PDBChain::ICToCoord(icx);
		D.y = PDBChain::ICToCoord(icy);
		D.z = PDBChain::ICToCoord(icz);
		}

	else if (Delta_tor2.m_tor2d == TOR2D_set_coords)
		{
		D.x = PDBChain::ICToCoord(Delta_tor2.m_ix);
		D.y = PDBChain::ICToCoord(Delta_tor2.m_iy);
		D.z = PDBChain::ICToCoord(Delta_tor2.m_iz);
		}

	else
		asserta(false);
	}

SqueezeDelta *Squeezer_tor2::EncodePos(
	const SqueezeState &State, uint i) const
	{
	SqueezeDelta_tor2 *Delta_tor2 = new SqueezeDelta_tor2;
	coords TrueD;
	m_Chain->GetCoords(i, TrueD);
	if (i < 3)
		{
		Delta_tor2->m_tor2d = TOR2D_set_coords;
		Delta_tor2->m_ix = PDBChain::CoordToIC(TrueD.x);
		Delta_tor2->m_iy = PDBChain::CoordToIC(TrueD.y);
		Delta_tor2->m_iz = PDBChain::CoordToIC(TrueD.z);
		return Delta_tor2;
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
		Delta_tor2->m_tor2d = TOR2D_tor;
		Delta_tor2->m_itheta = theta2i(theta_rad);
		Delta_tor2->m_iphi = phi2i(phi_rad);
		}
	else
		{
		int16_t icx = PDBChain::CoordToIC(State.m_C.x);
		int16_t icy = PDBChain::CoordToIC(State.m_C.y);
		int16_t icz = PDBChain::CoordToIC(State.m_C.z);

		int16_t idx = PDBChain::CoordToIC(TrueD.x);
		int16_t idy = PDBChain::CoordToIC(TrueD.y);
		int16_t idz = PDBChain::CoordToIC(TrueD.z);

		int fix_x = idx - icx;
		int fix_y = idy - icy;
		int fix_z = idz - icz;

		bool FixOk = true;
		FixOk = FixOk && (fix_x >= INT8_MIN && fix_x <= INT8_MAX);
		FixOk = FixOk && (fix_y >= INT8_MIN && fix_y <= INT8_MAX);
		FixOk = FixOk && (fix_z >= INT8_MIN && fix_z <= INT8_MAX);
		if (FixOk)
			{
			Delta_tor2->m_fix_x = int8_t(fix_x);
			Delta_tor2->m_fix_y = int8_t(fix_y);
			Delta_tor2->m_fix_z = int8_t(fix_z);
			Delta_tor2->m_tor2d = TOR2D_fix_coords;
			}
		else
			{
			Delta_tor2->m_ix = PDBChain::CoordToIC(TrueD.x);
			Delta_tor2->m_iy = PDBChain::CoordToIC(TrueD.y);
			Delta_tor2->m_iz = PDBChain::CoordToIC(TrueD.z);
			Delta_tor2->m_tor2d = TOR2D_set_coords;
			}
		}

	return Delta_tor2;
	}
