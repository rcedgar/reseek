#pragma once

#pragma once

#include "flat_base.h"

class museq_t : public flat_vec<char, FE_museq>
	{
protected:
    museq_t(uint32_t L) : flat_vec<char, FE_museq>(L) {}

public:
	static museq_t *newflat(uint32_t L)
		{
		museq_t *p = new museq_t(L);
		return p;
		}
	};
