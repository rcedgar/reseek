#pragma once

#include "squeezer.h"

class SqueezeDelta_null : public SqueezeDelta
	{
public:
	SqueezeDelta_null()
		{
		m_Name = "null";
		}

public:
	coords m_TrueD;

public:
	virtual uint GetBytes() const { return 3*sizeof(float); }
	};

class Squeezer_null : public Squeezer
	{
public:
	Squeezer_null()
		{
		m_Name = "null";
		}

public:
	virtual SqueezeState *NewState() const
		{
		return new SqueezeState;
		}

	virtual void DecodePos(const SqueezeState &State,
							const SqueezeDelta &Delta,
							uint i,
							coords &D) const
		{
		const SqueezeDelta_null &Delta_null =
			(const SqueezeDelta_null &) Delta;
		D = Delta_null.m_TrueD;
		}

	virtual SqueezeDelta *EncodePos(const SqueezeState &, uint i) const
		{
		SqueezeDelta_null *Delta = new SqueezeDelta_null;
		m_Chain->GetCoords(i, Delta->m_TrueD);
		return Delta;
		}

	virtual void UpdateState(
		SqueezeState &State,
		const SqueezeDelta &Delta) const {}
	};
