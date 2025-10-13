#pragma once

template <bool AllowUndefLetter> 
	uint ValueToIntTpl(float Value, uint AlphaSize,
		const vector<float> &Ts, uint DefaultLetter)
	{
	if (Value == FLT_MAX)
		{
		if (!AllowUndefLetter)
			asserta(DefaultLetter < AlphaSize);
		return DefaultLetter;
		}

	asserta(SIZE(Ts) + 1 == AlphaSize);
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		if (Value < Ts[i])
			return i;
	return AlphaSize - 1;
	}
