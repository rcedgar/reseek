#pragma once

#include "squeezer.h"

void calculate_theta_phi(const coords &A, const coords &B,
						 const coords &C, const coords &D,
						 float &theta_radians, float &phi_radians);

void calculate_D(const coords &A, const coords &B, const coords &C,
				 float theta_radians, float phi_radians,
				 coords &D);

enum TORD
	{
	TORD_undefined,
	TORD_tor,
	TORD_set_coords,
	};

class SqueezeDelta_tor : public SqueezeDelta
	{
public:
	SqueezeDelta_tor()
		{
		m_Name = "tor";
		}

public:
	TORD m_tord = TORD_undefined;
	float m_theta_rad = FLT_MAX;
	float m_phi_rad = FLT_MAX;
	float m_x = FLT_MAX;
	float m_y = FLT_MAX;
	float m_z = FLT_MAX;

public:
	virtual uint GetBytes() const
		{
		switch (m_tord)
			{
		case TORD_tor:	return 2*sizeof(float);
		case TORD_set_coords: return 3*sizeof(float);
			}
		Die("SqueezeDelta_tor::GetBytes m_tord=%d", m_tord);
		return 0;
		}

	virtual void LogMe() const
		{
		Log("D_tor(%d)", int(m_tord));
		switch (m_tord)
			{
		case TORD_tor: Log(" theta=%.1f, phi=%.1f\n", degrees(m_theta_rad), degrees(m_phi_rad)); return;
		case TORD_set_coords: Log(" SET x=%.1f y=%.1f z=%.1f\n", m_x, m_y, m_z); return;
			}
		Die("bad");
		}
	};

class Squeezer_tor : public Squeezer
	{
public:
	Squeezer_tor()
		{
		m_Name = "tor";
		}

public:
	virtual void DecodePos(const SqueezeState &State,
						   const SqueezeDelta &Delta,
						   coords &D) const;


	virtual SqueezeDelta *EncodePos(const SqueezeState &State,
									uint i) const;

public:
	void DecodePos_tor(const SqueezeState &State,
						   const SqueezeDelta_tor &Delta,
						   coords &D) const;
	};
