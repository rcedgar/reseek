#pragma once

#include "squeezer.h"

void calculate_theta_phi(const coords &A, const coords &B,
						 const coords &C, const coords &D,
						 float &theta_radians, float &phi_radians);

void calculate_D(const coords &A, const coords &B, const coords &C,
				 float theta_radians, float phi_radians,
				 coords &D);

enum TOR2D
	{
	TOR2D_undefined,
	TOR2D_tor,
	TOR2D_set_coords,
	};

class SqueezeDelta_tor2 : public SqueezeDelta
	{
public:
	SqueezeDelta_tor2()
		{
		m_Name = "tor2";
		}

public:
	TOR2D m_tor2d = TOR2D_undefined;
	int8_t m_itheta = INT8_MAX;
	int8_t m_iphi = INT8_MAX;
	int16_t m_ix = INT16_MAX;
	int16_t m_iy = INT16_MAX;
	int16_t m_iz = INT16_MAX;

public:
	virtual uint GetBytes() const
		{
		switch (m_tor2d)
			{
		case TOR2D_tor:			return 2*sizeof(int8_t);
		case TOR2D_set_coords: return 3*sizeof(int16_t);
			}
		Die("SqueezeDelta_tor2::GetBytes m_tor2d=%d", m_tor2d);
		return 0;
		}

	virtual void LogMe() const
		{
		Die("TODO");
		}
	};

class SqueezeState_tor2 : public SqueezeState
	{
public:
	coords m_A;
	coords m_B;
	coords m_C;
	uint m_ResetCount;

	SqueezeState_tor2()
		{
		m_A.invalidate();
		m_B.invalidate();
		m_C.invalidate();
		m_ResetCount = 0;
		}

	void LogStats() const
		{
		Log("Resets %u\n", m_ResetCount);
		}

	virtual void LogMe() const
		{
		Log("S_tor2 ");
		m_A.logme("A", false);
		m_B.logme("B", false);
		m_C.logme("C", false);
		}
	};

class Squeezer_tor2 : public Squeezer
	{
public:
	Squeezer_tor2()
		{
		m_Name = "tor2";
		}

public:
	virtual SqueezeState *NewState() const
		{
		return new SqueezeState_tor2;
		}

	virtual void InitDecode(SqueezeState *State) const;

	virtual void DecodePos(const SqueezeState &State,
						   const SqueezeDelta &Delta,
						   uint i,
						   coords &D) const;


	virtual SqueezeDelta *EncodePos(const SqueezeState &State,
									uint i) const;

	virtual void UpdateState(SqueezeState &State,
						   const SqueezeDelta &Delta) const;

public:
	void DecodePos_tor2(const SqueezeState_tor2 &State,
						   const SqueezeDelta_tor2 &Delta,
						   uint i,
						   coords &D) const;
	SqueezeDelta_tor2 *EncodePos_tor2(const SqueezeState_tor2 &State_tor2,
									uint i) const;

	void UpdateState_tor2(SqueezeState_tor2 &State_tor2,
						   const SqueezeDelta_tor2 &Delta_tor2) const;

public:
	static float i2phi(int8_t iphi);
	static float i2theta(int8_t itheta);
	static int8_t phi2i(float phi);
	static int8_t theta2i(float theta);
	static void test_cvt();
	};
