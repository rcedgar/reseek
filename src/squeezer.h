#pragma once

#include "pdbchain.h"

enum SDOP
	{
#define s(x)	SDOP_##x,
#include "squeezeops.h"
	};

static const uint SDOP_count = 0 +
#define s(x)	+ 1
#include "squeezeops.h"
	;

static const char *op2str(SDOP s)
	{
	switch (s)
		{
#define s(x)	case SDOP_##x: return #x;
#include "squeezeops.h"
		}
	asserta(false);
	return "SDOP_**ERROR";
	}

static SDOP str2op(const char *s)
	{
	if (0) void(0) ;
#define s(x)	else if (!strcmp(s, #x)) return SDOP_##x;
#include "squeezeops.h"
	asserta(false);
	return SDOP_undefined;
	}

class SqueezeDelta
	{
public:
	SDOP m_op = SDOP_undefined;

	float m_theta_rad = FLT_MAX;
	float m_phi_rad = FLT_MAX;

	int8_t m_itheta = INT8_MAX;
	int8_t m_iphi = INT8_MAX;

	float m_x = FLT_MAX;
	float m_y = FLT_MAX;
	float m_z = FLT_MAX;

	int16_t m_ICx = INT16_MAX;
	int16_t m_ICy = INT16_MAX;
	int16_t m_ICz = INT16_MAX;

	int8_t m_ICdx = INT8_MAX;
	int8_t m_ICdy = INT8_MAX;
	int8_t m_ICdz = INT8_MAX;

	int8_t m_ICvx = INT8_MAX;
	int8_t m_ICvy = INT8_MAX;
	int8_t m_ICvsignz = INT8_MAX;

public:
	uint GetBits() const
		{
		switch (m_op)
			{
		case SDOP_set_coords:	return 8*3*sizeof(float);
		case SDOP_ICv:	return 8*2*sizeof(int8_t) + 1;
		case SDOP_ICd:	return 8*3*sizeof(int16_t);
		case SDOP_IC:	return 8*3*sizeof(int16_t);
		case SDOP_tor:	return 8*2*sizeof(float);
		case SDOP_tor2:	return 8*2*sizeof(int8_t);
			}
		asserta(false);
		return 0;
		}

	void LogMe() const
		{
		Die("Not implemented");
		}
	};

class SqueezeState
	{
public:
	uint m_i;
	coords m_A;
	coords m_B;
	coords m_C;

public:
	SqueezeState()
		{
		m_i = 0;
		m_A.invalidate();
		m_B.invalidate();
		m_C.invalidate();
		}

	void LogMe() const
		{
		Log("State[%u]", m_i);
		m_A.logme("A", false);
		m_B.logme("B", false);
		m_C.logme("C", false);
		Log("\n");
		}
	};

class Squeezer
	{
public:
	string m_Algo = "invalid";
	const PDBChain *m_Chain = 0;
	vector<SqueezeDelta *> m_Deltas;
	float m_MaxErr = 0.1f;

public:
	SqueezeDelta *EncodePos(const SqueezeState &State, const coords &TrueD) const
		{
		if (State.m_i < 3)
			{
			SqueezeDelta *Delta = new SqueezeDelta;
			Delta->m_ICx = PDBChain::CoordToIC(TrueD.x);
			Delta->m_ICy = PDBChain::CoordToIC(TrueD.y);
			Delta->m_ICz = PDBChain::CoordToIC(TrueD.z);
			Delta->m_op = SDOP_IC;
			return Delta;
			}

#define	a(x)	if (m_Algo == #x) return EncodePos_algo_##x(State, TrueD);
#include "squeezealgos.h"
		asserta(false);
		return 0;
		}

#define a(x)	SqueezeDelta *EncodePos_algo_##x(const SqueezeState &State, const coords &TrueD) const;
#include "squeezealgos.h"

#define s(x)	void DecodePos_##x(const SqueezeState &State, const SqueezeDelta &Delta, coords &D) const;
#include "squeezeops.h"

	void DecodePos(const SqueezeState &State, const SqueezeDelta &Delta, coords &D) const
		{
		switch (Delta.m_op)
			{
#define s(x)	case SDOP_##x: DecodePos_##x(State, Delta, D); break;
#include "squeezeops.h"
		default: asserta(false);
			}
		}

public:
	void Clear();
	void EncodeChain(const PDBChain &Chain);
	void DecodeChain(PDBChain &Chain) const;
	void UpdateState(SqueezeState &State, const SqueezeDelta &Delta) const;
	uint GetBits() const;
	void LogMe() const;
	void GetErr(float &avge, float &maxe) const;
	void GetOpCounts(vector<uint> &Counts) const;
	void LogOpCounts() const;

public:
	static float i2phi(int8_t iphi);
	static float i2theta(int8_t itheta);
	static int8_t phi2i(float phi);
	static int8_t theta2i(float theta);
	static void test_cvt();
	};

void calculate_D(const coords &A, const coords &B, const coords &C,
				 float theta_radians, float phi_radians,
				 coords &D);
void calculate_theta_phi(const coords &A, const coords &B,
						 const coords &C, const coords &D,
						 float &theta_radians, float &phi_radians);
