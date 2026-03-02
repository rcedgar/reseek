#include "myutils.h"
#include "flat_base.h"
#include "chaq.h"

const chaindistmx_t *chaq::get_distmx()
	{
	if (m_distmx) return m_distmx;
	const uint32_t L = get_length();
	m_distmx = create_chaindistmx(L);
	return m_distmx;
	}
