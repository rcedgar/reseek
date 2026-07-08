#include "myutils.h"
#include <algorithm>
#include <string>
#include "parasail.h"

static void path_append_run(string &s, uint32_t count, char op)
	{
	if (count > 0)
		s.append(count, op);
	}

extern uint32_t parasail_cigar_encode(uint32_t length, char op_letter);

static const int PARASAIL_RCE_REQUIRED_FLAGS =
PARASAIL_FLAG_SW | PARASAIL_FLAG_STRIPED | PARASAIL_FLAG_TRACE |
PARASAIL_FLAG_BITS_16 | PARASAIL_FLAG_LANES_16;

static void parasail_result_verify_rce_source(const parasail_result_t *result)
	{
	asserta(result != NULL);
	asserta(parasail_result_is_trace(result));
	asserta(!(result->flag & PARASAIL_FLAG_SATURATED));
	asserta((result->flag & PARASAIL_RCE_REQUIRED_FLAGS) == PARASAIL_RCE_REQUIRED_FLAGS);
	asserta((result->flag & ~PARASAIL_RCE_REQUIRED_FLAGS) == 0);
	asserta(result->trace != NULL);
	asserta(result->trace->trace_table != NULL);
	}

/* Fold leading reseek D/I prefix into lo when parasail reports beg=0,0
   but encodes the offset as a leading run in the path (cigar.h quirk). */
static void path_fold_leading_prefix_reseek(uint &lo_i, uint &lo_j, string &path)
	{
	if (lo_i != 0 || lo_j != 0 || path.empty())
		return;
	if (path[0] == 'D')
		{
		size_t n = 0;
		while (n < path.size() && path[n] == 'D')
			++n;
		lo_i += (uint) n;
		path.erase(0, n);
		}
	else if (path[0] == 'I')
		{
		size_t n = 0;
		while (n < path.size() && path[n] == 'I')
			++n;
		lo_j += (uint) n;
		path.erase(0, n);
		}
	}

void parasail_result_to_path(
	parasail_result_t *result,
	int lena,
	int lenb,
	uint &lo_i,
	uint &lo_j,
	string &path)
	{
	string path_rev;
	uint32_t c_mat = 0;
	uint32_t c_del = 0;
	uint32_t c_ins = 0;
	int64_t i;
	int64_t j;
	int where;
	int16_t *HT;
	const int64_t segWidth = 16;
	int64_t segLen;
	bool ok = true;

	path.clear();
	lo_i = 0;
	lo_j = 0;

	parasail_result_verify_rce_source(result);

	asserta(lena > 0);
	asserta(lenb > 0);
	asserta(result->end_query >= 0 && result->end_query < lena);
	asserta(result->end_ref >= 0 && result->end_ref < lenb);

	segLen = (lena + segWidth - 1) / segWidth;
	HT = (int16_t *) result->trace->trace_table;

	path_rev.reserve((size_t) lena + (size_t) lenb);

	i = result->end_query;
	j = result->end_ref;
	where = PARASAIL_DIAG;

#define RESET_PATH \
do { \
	c_mat = 0; \
	c_del = 0; \
	c_ins = 0; \
} while (0)

#define WRITE_ANY_PATH \
do { \
	if (c_mat) \
		path_append_run(path_rev, c_mat, 'M'); \
	else if (c_del) \
		path_append_run(path_rev, c_del, 'D'); \
	else if (c_ins) \
		path_append_run(path_rev, c_ins, 'I'); \
	RESET_PATH; \
} while (0)

	while (ok && (i >= 0 || j >= 0)) {
		const int64_t loc = j * segLen * segWidth
			+ (i % segLen) * segWidth + (i / segLen);

		if (i < 0) {
			if (c_ins == 0)
				WRITE_ANY_PATH;
			while (j >= 0) {
				++c_ins;
				--j;
				}
			break;
			}
		if (j < 0) {
			if (c_del == 0)
				WRITE_ANY_PATH;
			while (i >= 0) {
				++c_del;
				--i;
				}
			break;
			}

		if (where == PARASAIL_DIAG) {
			if (HT[loc] & PARASAIL_DIAG) {
				if (c_mat == 0)
					WRITE_ANY_PATH;
				c_mat += 1;
				--i;
				--j;
				}
			else if (HT[loc] & PARASAIL_INS) {
				where = PARASAIL_INS;
				}
			else if (HT[loc] & PARASAIL_DEL) {
				where = PARASAIL_DEL;
				}
			else {
				break;
				}
			}
		else if (where == PARASAIL_INS) {
			if (c_ins == 0)
				WRITE_ANY_PATH;
			c_ins += 1;
			--j;
			if (HT[loc] & PARASAIL_DIAG_E)
				where = PARASAIL_DIAG;
			else if (HT[loc] & PARASAIL_INS_E)
				where = PARASAIL_INS;
			else
				ok = false;
			}
		else if (where == PARASAIL_DEL) {
			if (c_del == 0)
				WRITE_ANY_PATH;
			c_del += 1;
			--i;
			if (HT[loc] & PARASAIL_DIAG_F)
				where = PARASAIL_DIAG;
			else if (HT[loc] & PARASAIL_DEL_F)
				where = PARASAIL_DEL;
			else
				ok = false;
			}
		else if (where == PARASAIL_ZERO) {
			break;
			}
		else {
			ok = false;
			}
		}

#undef WRITE_ANY_PATH
#undef RESET_PATH

	if (!ok) {
		path.clear();
		lo_i = 0;
		lo_j = 0;
		return;
		}

	if (c_mat || c_del || c_ins) {
		if (c_mat)
			path_append_run(path_rev, c_mat, 'M');
		else if (c_del)
			path_append_run(path_rev, c_del, 'D');
		else if (c_ins)
			path_append_run(path_rev, c_ins, 'I');
		}

	path.assign(path_rev.rbegin(), path_rev.rend());

	lo_i = (uint) (i + 1);
	lo_j = (uint) (j + 1);

	path_fold_leading_prefix_reseek(lo_i, lo_j, path);
	}