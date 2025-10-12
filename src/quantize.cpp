#include "myutils.h"
#include "dss.h"

/***
Alphabet has AlphaSize (AS) letters.
Quantize by AS-1 thresholds
Undefined value can be FLT_MAX or in range.
***/

uint DSS::ValueToInt(float Value, uint AlphaSize, const vector<float> &Ts,
	uint DefaultLetter)
	{
	if (Value == FLT_MAX)
		{
		asserta(DefaultLetter < AlphaSize);
		return DefaultLetter;
		}

	asserta(SIZE(Ts) + 1 == AlphaSize);
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		if (Value < Ts[i])
			return i;
	return AlphaSize - 1;
	}
