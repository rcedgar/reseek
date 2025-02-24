#pragma once

#include "squeezer.h"

class SqueezeDelta_int16 : public SqueezeDelta
	{
public:
	SqueezeDelta_int16()
		{
		m_Name = "null";
		}

public:
	int16_t m_x16 = INT16_MAX;
	int16_t m_y16 = INT16_MAX;
	int16_t m_z16 = INT16_MAX;

public:
	virtual uint GetBytes() const { return 3*sizeof(int16_t); }
	};

class Squeezer_int16 : public Squeezer
	{
public:
	Squeezer_int16()
		{
		m_Name = "int16";
		}

public:
	virtual void InitDecode() {}

	SqueezeState *NewState() const
		{
		return new SqueezeState;
		}

	virtual void DecodePos(const SqueezeState &State,
							const SqueezeDelta &Delta,
							uint i,
							coords &D) const
		{
		const SqueezeDelta_int16 &Delta_int16 =
			(const SqueezeDelta_int16 &) Delta;
		D.x = PDBChain::ICToCoord(Delta_int16.m_x16);
		D.y = PDBChain::ICToCoord(Delta_int16.m_y16);
		D.z = PDBChain::ICToCoord(Delta_int16.m_z16);
		}

	virtual SqueezeDelta *EncodePos(const SqueezeState &State, uint i) const
		{
		SqueezeDelta_int16 *Delta_int16 = new SqueezeDelta_int16;
		coords TrueD;
		m_Chain->GetCoords(i, TrueD);
		Delta_int16->m_x16 = PDBChain::CoordToIC(TrueD.x);
		Delta_int16->m_y16 = PDBChain::CoordToIC(TrueD.y);
		Delta_int16->m_z16 = PDBChain::CoordToIC(TrueD.z);
		return Delta_int16;
		}

	virtual void UpdateState(
		SqueezeState &State,
		const SqueezeDelta &Delta) const {}
	};
