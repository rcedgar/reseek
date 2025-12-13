#include "myutils.h"
#include "hexintseq.h"
#include "alpha.h"

void ByteSeqToHexIntStr(uint AlphaSize, const vector<byte> &Letters, string &Str)
	{
	string Tmp;
	const uint L = SIZE(Letters);
	Str.clear();
	Str.reserve(L);
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		byte Letter = Letters[Pos];
		if (AlphaSize < 36)
			Str += g_LetterToCharMu[Letter];
		else
			{
			Ps(Tmp, "%02x", Letter);
			Str += Tmp;
			}
		}
	}
