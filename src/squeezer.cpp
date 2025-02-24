#include "myutils.h"
#include "squeezer.h"
#include "chainreader2.h"

/////////////////////////////////////////////
// Base class
/////////////////////////////////////////////
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
	m_Chain = &Chain;
	const uint L = Chain.GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		SqueezeDelta *Delta = EncodePos(i);
		m_Deltas.push_back(Delta);
		}
	}

void Squeezer::DecodeChain(PDBChain &Chain)
	{
	InitDecode();
	m_NextDecodePos = 0;
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
		DecodeNext(D);
		Chain.m_Xs.push_back(D.x);
		Chain.m_Ys.push_back(D.y);
		Chain.m_Zs.push_back(D.z);
		++m_NextDecodePos;
		}
	}

void cmd_squeeze()
	{
	ChainReader2 CR;
	CR.Open(g_Arg1);
	Squeezer_int16 S;
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
		Log("%s >%s avge %.3g maxe %.3g\n",
			S.m_Name.c_str(), Chain.m_Label.c_str(), sume/L, maxe);

		delete ptrChain;
		}
	}