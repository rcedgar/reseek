#include "myutils.h"
#include "squeezer_null.h"
#include "squeezer_int16.h"
#include "squeezer_tor.h"
#include "squeezer_tor2.h"
#include "chainreader2.h"

Squeezer *Squeezer::NewSqueezer(const string &Name)
	{
	if (Name == "null")
		return new Squeezer_null;
	else if (Name == "int16")
		return new Squeezer_int16;
	else if (Name == "tor")
		return new Squeezer_tor;
	else if (Name == "tor2")
		return new Squeezer_tor2;
	Die("NewSqueezer(%s)\n", Name.c_str());
	return 0;
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
	SqueezeState *State = NewState();
	m_Chain = &Chain;
	const uint L = Chain.GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		//Log("\n");
		//Log(" ========== %u =========\n", i);
		SqueezeDelta *Delta = EncodePos(*State, i);
		//Delta->LogMe();
		m_Deltas.push_back(Delta);
		UpdateState(*State, *Delta);
		//State->LogMe();
		}
	delete State;
	}

void Squeezer::DecodeChain(PDBChain &Chain)
	{
	SqueezeState *State = NewState();
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
		}
	delete State;
	}

uint Squeezer::GetBytes() const
	{
	const uint L = SIZE(m_Deltas);
	uint Bytes = 0;
	for (uint i = 0; i < L; ++i)
		Bytes += m_Deltas[i]->GetBytes();
	return Bytes;
	}

void cmd_squeeze()
	{
	ChainReader2 CR;
	CR.Open(g_Arg1);
	Squeezer &S = *Squeezer::NewSqueezer("tor2");
	for (;;)
		{
		PDBChain *ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		const PDBChain &Chain = *ptrChain;
		S.EncodeChain(Chain);

		PDBChain Chain2;
		Chain2.m_Label = Chain.m_Label;
		Chain2.m_Seq = Chain.m_Seq;
		S.DecodeChain(Chain2);
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
		uint Bytes = S.GetBytes();
		Log("%s >%s(%u) avge %.3g maxe %.3g bytes %u\n",
			S.m_Name.c_str(),
			Chain.m_Label.c_str(),
			L,
			sume/L, maxe,
			Bytes);

		delete ptrChain;
		}
	}

void Squeezer::UpdateState(SqueezeState &State, const SqueezeDelta &Delta) const
	{
	coords D;
	DecodePos(State, Delta, D);

	State.m_A = State.m_B;
	State.m_B = State.m_C;
	State.m_C = D;
	}
