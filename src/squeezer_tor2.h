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
	TOR2D_fix_coords,
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
	int8_t m_fix_x = INT8_MAX;
	int8_t m_fix_y = INT8_MAX;
	int8_t m_fix_z = INT8_MAX;

public:
	virtual uint GetBytes() const
		{
		switch (m_tor2d)
			{
		case TOR2D_tor:			return 2*sizeof(int8_t);
		case TOR2D_set_coords:	return 3*sizeof(int16_t);
		case TOR2D_fix_coords:	return 3*sizeof(int8_t);
			}
		Die("SqueezeDelta_tor2::GetBytes m_tor2d=%d", m_tor2d);
		return 0;
		}

	virtual void LogMe() const
		{
		Die("TODO");
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
	virtual void DecodePos(const SqueezeState &State,
						   const SqueezeDelta &Delta,
						   coords &D) const;


	virtual SqueezeDelta *EncodePos(const SqueezeState &State,
									uint i) const;

public:
	void DecodePos_tor2(const SqueezeState &State,
						   const SqueezeDelta_tor2 &Delta,
						   coords &D) const;

public:
	static float i2phi(int8_t iphi);
	static float i2theta(int8_t itheta);
	static int8_t phi2i(float phi);
	static int8_t theta2i(float theta);
	static void test_cvt();
	};
