#include "myutils.h"
#include "flat_base.h"
#include "chaq.h"
#include "fast_dist_mx2.h"

static const BandIndexLite s_bi(M);

const chaindistmx_t *chaq::get_distmx()
	{
	if (m_distmx) return m_distmx;
	const uint32_t L = get_length();
	m_distmx = create_chaindistmx(L);
	const uint16_t *chain_ics = m_chain->m_xyz->m_data;
	uint16_t *distmx = m_distmx->m_data;
	banded_distances_avx2_u16_xyz(chain_ics, L, s_bi, distmx);
	return m_distmx;
	}
