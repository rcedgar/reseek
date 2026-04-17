#pragma once

#include "flat_aligner.h"

class flat_alignx
	{
public:

public:
	static float alignx(
		const flat_aligner &fa,
		const uint8_t *profQ,
		const uint8_t *profT,
		const sid_t *distmxQ,
		const sid_t *distmxT,
		float selfT,
		float selfQ,
		uint M);
	};