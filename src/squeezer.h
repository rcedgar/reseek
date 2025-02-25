#pragma once

#include "pdbchain.h"

class SqueezeDelta
	{
public:
	string m_Name;

public:
	virtual uint GetBytes() const = 0;
	virtual void LogMe() const { Die("Not implemented"); }
	};

class SqueezeState
	{
public:
	coords m_A;
	coords m_B;
	coords m_C;

public:
	SqueezeState()
		{
		m_A.invalidate();
		m_B.invalidate();
		m_C.invalidate();
		}

	void LogMe() const
		{
		Log("State: ");
		m_A.logme("A", false);
		m_B.logme("B", false);
		m_C.logme("C", false);
		Log("\n");
		}
	};

class Squeezer
	{
public:
	string m_Name;
	const PDBChain *m_Chain = 0;
	vector<SqueezeDelta *> m_Deltas;
	float m_MaxErr = 0.1f;

public:
	virtual SqueezeDelta *EncodePos(
		const SqueezeState &State,
		uint i) const = 0;

	virtual void DecodePos(
		const SqueezeState &State,
		const SqueezeDelta &Delta,
		coords &D) const = 0;

public:
	SqueezeState *NewState() const
		{
		return new SqueezeState;
		}

	void GetTrueD(uint i, coords &TrueD) const
		{
		m_Chain->GetCoords(i, TrueD);
		}

	float GetErr(const coords &D, const coords &TrueD) const
		{
		return coords::dist(D, TrueD);
		}

public:
	void Clear();
	void EncodeChain(const PDBChain &Chain);
	void DecodeChain(PDBChain &Chain);
	uint GetBytes() const;
	void UpdateState(SqueezeState &State, const SqueezeDelta &Delta) const;

public:
	static Squeezer *NewSqueezer(const string &Name);
	};
