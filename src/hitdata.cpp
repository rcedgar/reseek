#include "myutils.h"
#include "flat_chain.h"
#include "hitdata.h"
#include "cigar.h"

void hitdata::fill()
	{
	ids = 0;
	diffs = 0;
	gaps = 0;
	const char *qaa = query->m_aa->m_data;
	const char *taa = target->m_aa->m_data;
	const uint LQ = query->m_L;
	const uint LT = target->m_L;
	cigar.clear();
	cigar.reserve(2*ncol);
	uint qpos = qlo;
	uint tpos= tlo;
	char lastc = path[0];
	assert(lastc == 'M');
	uint cigaropn = 0;
	if (qlo > 0)
		Psa(cigar, "%uS", qlo);
	for (uint i = 0; i < ncol; ++i)
		{
		char c = path[i];
		if (c == lastc)
			++cigaropn;
		else
			{
			Psa(cigar, "%u%c", cigaropn, lastc);
			cigaropn = 1;
			}

		switch (c)
			{
		case 'M':
			{
			assert(qpos < LQ);
			assert(tpos < LT);
			char q = qaa[qpos];
			char t = taa[tpos];
			if (toupper(q) == toupper(t))
				++ids;
			else
				++diffs;
			++qhi;
			++thi;
			break;
			}

		case 'D':
			{
			++gaps;
			++qhi;
			break;
			}

		case 'I':
			{
			++thi;
			++gaps;
			break;
			}
		default: asserta(false);
			}
		lastc = c;
		}
	Psa(cigar, "%u%c", cigaropn, lastc);
	qhi = qpos - 1;
	thi = tpos - 1;
	}
