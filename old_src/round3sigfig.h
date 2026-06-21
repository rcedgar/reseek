#pragma once

static inline float round3sigfig(float x)
	{
	char s[16];
	sprintf(s, "%.3g", x);
	float rounded = (float) atof(s);
	return rounded;
	}
