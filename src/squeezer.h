#pragma once

#include "pdbchain.h"

class SqueezeDelta
	{
public:
	string m_Name;
	};

class SqueezeDelta_null : public SqueezeDelta
	{
public:
	SqueezeDelta_null()
		{
		m_Name = "null";
		}

public:
	coords m_TrueD;
	};

class SqueezeDelta_int16 : public SqueezeDelta
	{
public:
	SqueezeDelta_int16()
		{
		m_Name = "null";
		}

public:
	int16_t m_x16;
	int16_t m_y16;
	int16_t m_z16;
	};

class Squeezer
	{
public:
	string m_Name;
	const PDBChain *m_Chain = 0;
	vector<SqueezeDelta *> m_Deltas;
	uint m_NextDecodePos = UINT_MAX;

public:
	virtual uint GetBits() const = 0;
	virtual SqueezeDelta *EncodePos(uint i) = 0;
	virtual void InitDecode() const = 0;
	virtual void DecodeNext(coords &D) const = 0;

public:
	void Clear();
	void EncodeChain(const PDBChain &Chain);
	void DecodeChain(PDBChain &Chain);
	};

class Squeezer_null : public Squeezer
	{
public:
	Squeezer_null()
		{
		m_Name = "null";
		}

public:
	virtual uint GetBits() const { return 3*sizeof(float); }

	virtual SqueezeDelta *EncodePos(uint i)
		{
		SqueezeDelta_null *Delta = new SqueezeDelta_null;
		m_Chain->GetCoords(i, Delta->m_TrueD);
		return Delta;
		}

	virtual void DecodePos(const SqueezeDelta *Delta, coords &D) const
		{
		const SqueezeDelta_null *Delta_null = (const SqueezeDelta_null *) Delta;
		D = Delta_null->m_TrueD;
		}

	virtual void InitDecode() const
		{
		}

	virtual void DecodeNext(coords &D) const
		{
		assert(m_NextDecodePos < SIZE(m_Deltas));
		const SqueezeDelta_null *Delta =
			(const SqueezeDelta_null *) m_Deltas[m_NextDecodePos];
		D = Delta->m_TrueD;
		}
	};

class Squeezer_int16 : public Squeezer
	{
public:
	Squeezer_int16()
		{
		m_Name = "int16";
		}

public:
	virtual uint GetBits() const { return 3*sizeof(int16_t); }

	virtual void DecodePos(const SqueezeDelta *Delta, coords &D) const
		{
		const SqueezeDelta_int16 *Delta_int16 =
			(const SqueezeDelta_int16 *) Delta;
		D.x = PDBChain::ICToCoord(Delta_int16->m_x16);
		D.y = PDBChain::ICToCoord(Delta_int16->m_y16);
		D.z = PDBChain::ICToCoord(Delta_int16->m_z16);
		}

	virtual SqueezeDelta *EncodePos(uint i)
		{
		SqueezeDelta_int16 *Delta_int16 = new SqueezeDelta_int16;
		coords TrueD;
		m_Chain->GetCoords(i, TrueD);
		Delta_int16->m_x16 = PDBChain::CoordToIC(TrueD.x);
		Delta_int16->m_y16 = PDBChain::CoordToIC(TrueD.y);
		Delta_int16->m_z16 = PDBChain::CoordToIC(TrueD.z);
		return Delta_int16;
		}

	virtual void InitDecode() const {}

	virtual void DecodeNext(coords &D) const
		{
		assert(m_NextDecodePos < SIZE(m_Deltas));
		const SqueezeDelta_int16 *Delta =
			(const SqueezeDelta_int16 *) m_Deltas[m_NextDecodePos];
		D.x = PDBChain::ICToCoord(Delta->m_x16);
		D.y = PDBChain::ICToCoord(Delta->m_y16);
		D.z = PDBChain::ICToCoord(Delta->m_z16);
		}
	};
