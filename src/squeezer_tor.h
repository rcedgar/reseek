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
	TORD_delta_coords,
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
	float m_dx = FLT_MAX;
	float m_dy = FLT_MAX;
	float m_dz = FLT_MAX;

public:
	virtual uint GetBytes() const
		{
		switch (m_tord)
			{
		case TORD_tor:	return 2*sizeof(float);
		case TORD_delta_coords: return 3*sizeof(float);
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
		case TORD_delta_coords: Log(" ERROR\n", degrees(m_theta_rad), degrees(m_phi_rad)); return;
		case TORD_set_coords: Log(" SET x=%.1f y=%.1f z=%.1f\n", m_dx, m_dy, m_dz); return;
			}
		Die("bad");
		}
	};

class SqueezeState_tor : public SqueezeState
	{
public:
	coords m_A;
	coords m_B;
	coords m_C;
	uint m_ResetCount;

	SqueezeState_tor()
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
		Log("S_tor ");
		m_A.logme("A", false);
		m_B.logme("B", false);
		m_C.logme("C", false);
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
	virtual SqueezeState *NewState() const
		{
		return new SqueezeState_tor;
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
	void DecodePos_tor(const SqueezeState_tor &State,
						   const SqueezeDelta_tor &Delta,
						   uint i,
						   coords &D) const;
	SqueezeDelta_tor *EncodePos_tor(const SqueezeState_tor &State_tor,
									uint i) const;

	void UpdateState_tor(SqueezeState_tor &State_tor,
						   const SqueezeDelta_tor &Delta_tor) const;
	};
