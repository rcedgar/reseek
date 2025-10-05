#pragma once

enum UNDEF_BINNING
	{
	UB_Invalid,
	UB_NeverUndefined,
	UB_IgnoreUndefined,
	UB_IncludeUndefined,
	UB_UndefinedIsDefaultLetter,
	UB_UndefinedIsOnlyZero,
	UB_UndefinedIsZeroOverload,
	UB_IntFeatureNoBinning,
	};

UNDEF_BINNING StrToUB(const string &s);
const char *UBToStr(UNDEF_BINNING UB);