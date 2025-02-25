#include "myutils.h"
#include "squeezer.h"
#include "chainreader2.h"
#include "mypi.h"

// theta	 25 .. 100 deg.		0.4363323 .. 1.745329 radians.
// phi	   -180 .. 180 deg.		-PI .. +PI

static const float MIN_THETA = 0.436f;
static const float MAX_THETA = 1.746f;

static const float MIN_PHI = -PI;
static const float MAX_PHI = PI;

void Squeezer::test_cvt()
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

int8_t Squeezer::theta2i(float theta)
	{
	float zero_to_one = (theta - MIN_THETA)/(MAX_THETA - MIN_THETA);
	float minus_one_to_plus_one = 2.0f*(zero_to_one - 0.5f);
	float minus_127_to_plus_127 = 127*minus_one_to_plus_one;
	int8_t i = int8_t(minus_127_to_plus_127);
	return i;
	}

int8_t Squeezer::phi2i(float phi)
	{
	float zero_to_one = (phi - MAX_PHI)/(MAX_PHI - MIN_PHI);
	float minus_one_to_plus_one = 2.0f*(zero_to_one - 0.5f);
	float minus_127_to_plus_127 = 127*minus_one_to_plus_one;
	int8_t i = int8_t(minus_127_to_plus_127);
	return i;
	}

float Squeezer::i2theta(int8_t itheta)
	{
	float minus_127_to_plus_127 = float(itheta);
	float minus_one_to_plus_one = minus_127_to_plus_127/127.0f;
	float zero_to_one = 0.5f + minus_one_to_plus_one/2.0f;
	float theta = MIN_THETA + zero_to_one*(MAX_THETA - MIN_THETA);
	return theta;
	}

float Squeezer::i2phi(int8_t iphi)
	{
	float minus_127_to_plus_127 = float(iphi);
	float minus_one_to_plus_one = minus_127_to_plus_127/127.0f;
	float zero_to_one = 0.5f + minus_one_to_plus_one/2.0f;
	float phi = MIN_PHI + zero_to_one*(MAX_PHI - MIN_PHI);
	return phi;
	}

SqueezeDelta *Squeezer::EncodePos_algo_undefined(
	const SqueezeState &State, const coords &TrueD) const
	{
	asserta(false);
	return 0;
	}

SqueezeDelta *Squeezer::EncodePos_algo_ICd(
	const SqueezeState &State, const coords &TrueD) const
	{
	SqueezeDelta *Delta = new SqueezeDelta;

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
	return Delta;
	}

SqueezeDelta *Squeezer::EncodePos_algo_ICv(
	const SqueezeState &State, const coords &TrueD) const
	{
	SqueezeDelta *Delta = new SqueezeDelta;

	int16_t icx = PDBChain::CoordToIC(State.m_C.x);
	int16_t icy = PDBChain::CoordToIC(State.m_C.y);

	int16_t itruedx = PDBChain::CoordToIC(TrueD.x);
	int16_t itruedy = PDBChain::CoordToIC(TrueD.y);

	int8_t fix_x = itruedx - icx;
	int8_t fix_y = itruedy - icy;

	int16_t idx = icx + fix_x;
	int16_t idy = icy + fix_y;

	coords TestD;
	TestD.x = PDBChain::ICToCoord(idx);
	TestD.y = PDBChain::ICToCoord(idy);

// 3.81 = sqrt(dx*dx + dy*dy + dz*dz)
// 3.81^2 = dx*dx + dy*dy + dz*dz
// dz^2 = 3.81^2 - dx*dx - dy*dy
	float dx = TestD.x - State.m_C.x;
	float dy = TestD.y - State.m_C.y;
	float dz2 = 3.81f*3.81f - dx*dx - dy*dy;
	if (dz2 < 0)
		dz2 = 0;
	TestD.z = State.m_C.z + sqrt(dz2);
	float d = coords::dist(TestD, State.m_C);
	float teste = coords::dist(TrueD, TestD);
	if (teste <= m_MaxErr)
		{
		Delta->m_ICvx = fix_x;
		Delta->m_ICvy = fix_y;
		Delta->m_ICvsignz = +1;
		Delta->m_op = SDOP_ICv;
#if DEBUG
		{
		coords TestD2;
		DecodePos_ICv(State, *Delta, TestD2);
		float e2 = coords::dist(TrueD, TestD2);
		asserta(e2 <= m_MaxErr);
		}
#endif
		return Delta;
		}

	TestD.z = State.m_C.z - sqrt(dz2);
	d = coords::dist(TestD, State.m_C);
	float teste_minus = coords::dist(TrueD, TestD);
	if (teste_minus <= m_MaxErr)
		{
		Delta->m_ICvx = fix_x;
		Delta->m_ICvy = fix_y;
		Delta->m_ICvsignz = -1;
		Delta->m_op = SDOP_ICv;
#if DEBUG
		{
		coords TestD2;
		DecodePos_ICv(State, *Delta, TestD2);
		float e2 = coords::dist(TrueD, TestD2);
		asserta(e2 <= m_MaxErr);
		}
#endif
		return Delta;
		}

	Delta->m_ICx = PDBChain::CoordToIC(TrueD.x);
	Delta->m_ICy = PDBChain::CoordToIC(TrueD.y);
	Delta->m_ICz = PDBChain::CoordToIC(TrueD.z);
	Delta->m_op = SDOP_IC;
	return Delta;
	}

SqueezeDelta *Squeezer::EncodePos_algo_IC(
	const SqueezeState &State, const coords &TrueD) const
	{
	SqueezeDelta *Delta = new SqueezeDelta;
	Delta->m_op = SDOP_IC;
	Delta->m_ICx = PDBChain::CoordToIC(TrueD.x);
	Delta->m_ICy = PDBChain::CoordToIC(TrueD.y);
	Delta->m_ICz = PDBChain::CoordToIC(TrueD.z);
	return Delta;
	}

SqueezeDelta *Squeezer::EncodePos_algo_null(
	const SqueezeState &State, const coords &TrueD) const
	{
	SqueezeDelta *Delta = new SqueezeDelta;
	Delta->m_op = SDOP_set_coords;
	Delta->m_x = TrueD.x;
	Delta->m_y = TrueD.y;
	Delta->m_z = TrueD.z;
	return Delta;
	}

void Squeezer::DecodePos_undefined(const SqueezeState &State,
						const SqueezeDelta &Delta,
						coords &D) const
	{
	asserta(false);
	}

void Squeezer::DecodePos_IC(const SqueezeState &State,
						const SqueezeDelta &Delta,
						coords &D) const
	{
	asserta(Delta.m_op == SDOP_IC);
	D.x = PDBChain::ICToCoord(Delta.m_ICx);
	D.y = PDBChain::ICToCoord(Delta.m_ICy);
	D.z = PDBChain::ICToCoord(Delta.m_ICz);
	}

void Squeezer::DecodePos_ICd(const SqueezeState &State,
						const SqueezeDelta &Delta,
						coords &D) const
	{
	asserta(Delta.m_op == SDOP_ICd);

	int16_t icx = PDBChain::CoordToIC(State.m_C.x);
	int16_t icy = PDBChain::CoordToIC(State.m_C.y);
	int16_t icz = PDBChain::CoordToIC(State.m_C.z);

	D.x = PDBChain::ICToCoord(icx + Delta.m_ICdx);
	D.y = PDBChain::ICToCoord(icy + Delta.m_ICdy);
	D.z = PDBChain::ICToCoord(icz + Delta.m_ICdz);
	}

void Squeezer::DecodePos_ICv(const SqueezeState &State,
						const SqueezeDelta &Delta,
						coords &D) const
	{
	asserta(Delta.m_op == SDOP_ICv);

	int16_t icx = PDBChain::CoordToIC(State.m_C.x);
	int16_t icy = PDBChain::CoordToIC(State.m_C.y);

	D.x = PDBChain::ICToCoord(icx + Delta.m_ICvx);
	D.y = PDBChain::ICToCoord(icy + Delta.m_ICvy);

	float dx = D.x - State.m_C.x;
	float dy = D.y - State.m_C.y;
	float dz2 = 3.81f*3.81f - dx*dx - dy*dy;
	if (dz2 < 0)
		dz2 = 0;

	D.z = State.m_C.z + sqrt(dz2)*Delta.m_ICvsignz;
	}

void Squeezer::DecodePos_set_coords(const SqueezeState &State,
						const SqueezeDelta &Delta,
						coords &D) const
	{
	asserta(Delta.m_op == SDOP_set_coords);
	D.x = Delta.m_x;
	D.y = Delta.m_y;
	D.z = Delta.m_z;
	}

void Squeezer::DecodePos_tor(const SqueezeState &State,
						const SqueezeDelta &Delta,
						coords &D) const
	{
	asserta(Delta.m_op == SDOP_tor);
	State.m_A.assert_valid();
	State.m_B.assert_valid();
	State.m_C.assert_valid();
	calculate_D(State.m_A,
				State.m_B,
				State. m_C,
				Delta.m_theta_rad,
				Delta.m_phi_rad,
				D);
	}

void Squeezer::DecodePos_tor2(const SqueezeState &State,
						const SqueezeDelta &Delta,
						coords &D) const
	{
	asserta(Delta.m_op == SDOP_tor2);

	float theta_rad = i2theta(Delta.m_itheta);
	float phi_rad = i2phi(Delta.m_iphi);

	State.m_A.assert_valid();
	State.m_B.assert_valid();
	State.m_C.assert_valid();
	calculate_D(State.m_A,
				State.m_B,
				State.m_C,
				theta_rad,
				phi_rad,
				D);
	}

uint Squeezer::GetBits() const
	{
	uint Bits = 0;
	for (uint i = 0; i < SIZE(m_Deltas); ++i)
		Bits += m_Deltas[i]->GetBits();
	return Bits;
	}

void Squeezer::Clear()
	{
	const uint L = SIZE(m_Deltas);
	for (uint i = 0; i < L; ++i)
		delete m_Deltas[i];
	m_Deltas.clear();
	}

void Squeezer::EncodeChain(const PDBChain &Chain)
	{
	Clear();
	SqueezeState *State = new SqueezeState;
	m_Chain = &Chain;
	const uint L = Chain.GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		coords TrueD;
		Chain.GetCoords(i, TrueD);
		SqueezeDelta *Delta = EncodePos(*State, TrueD);
		m_Deltas.push_back(Delta);
		UpdateState(*State, *Delta);
		State->m_i = i;
		}
	delete State;
	}

void Squeezer::DecodeChain(PDBChain &Chain) const
	{
	SqueezeState *State = new SqueezeState;
	const uint L = SIZE(m_Deltas);
	Chain.m_Xs.reserve(L);
	Chain.m_Ys.reserve(L);
	Chain.m_Zs.reserve(L);

	Chain.m_Xs.clear();
	Chain.m_Ys.clear();
	Chain.m_Zs.clear();

	coords D;
	for (uint i = 0; i < L; ++i)
		{
		const SqueezeDelta *Delta = m_Deltas[i];
		DecodePos(*State, *Delta, D);
		Chain.m_Xs.push_back(D.x);
		Chain.m_Ys.push_back(D.y);
		Chain.m_Zs.push_back(D.z);

		UpdateState(*State, *Delta);
		State->m_i = i;
		}
	delete State;
	}

void Squeezer::UpdateState(SqueezeState &State, const SqueezeDelta &Delta) const
	{
	coords D;
	DecodePos(State, Delta, D);

	State.m_A = State.m_B;
	State.m_B = State.m_C;
	State.m_C = D;
	}

void Squeezer::GetErr(float &avge, float &maxe) const
	{
	const PDBChain &Chain = *m_Chain;
	PDBChain Chain2;
	Chain2.m_Label = Chain.m_Label;
	Chain2.m_Seq = Chain.m_Seq;
	DecodeChain(Chain2);
	const uint L = Chain.GetSeqLength();
	const uint L2 = Chain2.GetSeqLength();
	asserta(L2 == L);
	float sume = 0;
	for (uint i = 0; i < L; ++i)
		{
		float x = Chain.m_Xs[i];
		float y = Chain.m_Ys[i];
		float z = Chain.m_Zs[i];

		float x2 = Chain2.m_Xs[i];
		float y2 = Chain2.m_Ys[i];
		float z2 = Chain2.m_Zs[i];

		float ex = x2 - x;
		float ey = y2 - y;
		float ez = z2 - z;
		float e = sqrt(ex*ex + ey*ey + ez*ez);
		sume += e;
		maxe = max(e, maxe);
		}
	avge = maxe/L;
	}

void Squeezer::LogMe() const
	{
	const PDBChain &Chain = *m_Chain;
	PDBChain Chain2;
	Chain2.m_Label = Chain.m_Label;
	Chain2.m_Seq = Chain.m_Seq;
	DecodeChain(Chain2);
	const uint L = Chain.GetSeqLength();
	const uint L2 = Chain2.GetSeqLength();
	asserta(L2 == L);
	float maxe = 0;
	float sume = 0;
	for (uint i = 0; i < L; ++i)
		{
		float x = Chain.m_Xs[i];
		float y = Chain.m_Ys[i];
		float z = Chain.m_Zs[i];

		float x2 = Chain2.m_Xs[i];
		float y2 = Chain2.m_Ys[i];
		float z2 = Chain2.m_Zs[i];

		float ex = x2 - x;
		float ey = y2 - y;
		float ez = z2 - z;
		float e = sqrt(ex*ex + ey*ey + ez*ez);
		sume += e;
		maxe = max(e, maxe);

		Log("%8.1f", x);
		Log("  %8.1f", x2);
		Log(" |  %8.1f", y);
		Log("  %8.1f", y2);
		Log("  |  %8.1f", z);
		Log("  %8.1f", z2);
		Log("  > %3.1f", e);
		Log("\n");
		}
	}

void Squeezer::GetOpCounts(vector<uint> &Counts) const
	{
	Counts.clear();
	Counts.resize(SDOP_count);
	const uint L = SIZE(m_Deltas);
	for (uint i = 0; i < L; ++i)
		Counts[m_Deltas[i]->m_op] += 1;
	}

void Squeezer::LogOpCounts() const
	{
	vector<uint> Counts;
	GetOpCounts(Counts);
	for (uint op = 0; op < SDOP_count; ++op)
		{
		uint n = Counts[op];
		if (n == 0)
			continue;
		Log(" %s/%u", op2str(SDOP(op)), n);
		}
	}

static void LogCADists(const PDBChain &Chain)
	{
	Log("\n");
	const uint L = Chain.GetSeqLength();
	Log(">%s(%u)\n", Chain.m_Label.c_str(), L);
	for (uint i = 0; i < L; ++i)
		{
		float x = Chain.m_Xs[i];
		float y = Chain.m_Ys[i];
		float z = Chain.m_Zs[i];

		Log("[%5u]  %6.1f  %8.1f  %8.1f", i, x, y, z);
		if (i > 0)
			{
			float x1 = Chain.m_Xs[i-1];
			float y1 = Chain.m_Ys[i-1];
			float z1 = Chain.m_Zs[i-1];
			float d = Chain.GetDist(i, i-1);
			Log("  %3.1f", d);
			Log(" (%.2f, %.2f, %.2f)",
				x - x1, 
				y - y1, 
				z - z1);
			}
		Log("\n");
		}
	}

void cmd_squeeze()
	{
	ChainReader2 CR;
	CR.Open(g_Arg1);

	vector<string> Algos;
#define a(x)	Algos.push_back(#x);
#include "squeezealgos.h"
	const uint AlgoCount = SIZE(Algos);

	Squeezer S;
	for (;;)
		{
		PDBChain *ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		const PDBChain &Chain = *ptrChain;
		LogCADists(Chain);
		for (uint ia = 0; ia < AlgoCount; ++ia)
			{
			S.m_Algo = Algos[ia];
			if (S.m_Algo == "undefined")
				continue;
			S.EncodeChain(Chain);
			uint Bytes = S.GetBits()/8;
			float avge, maxe;
			S.GetErr(avge, maxe);
			Log("%8.8s", S.m_Algo.c_str());
			Log("   avge %9.3g", avge);
			Log("   maxe %9.3g", maxe);
			Log("   bytes %5u", Bytes);
			S.LogOpCounts();
			Log("\n");
			}

		delete ptrChain;
		}
	}
