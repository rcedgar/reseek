#include "myutils.h"
#include <stdint.h>
#include <stdlib.h>
#include "parasail.h"

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

parasail_cigar_t* parasail_result_get_cigar_rce(
	parasail_result_t *result,
	int lena,
	int lenb)
	{
	size_t size;
	parasail_cigar_t *cigar;
	uint32_t *cigar_reverse;
	uint32_t c_mat = 0;
	uint32_t c_del = 0;
	uint32_t c_ins = 0;
	int64_t i;
	int64_t j;
	int where;
	int16_t *HT;
	const int64_t segWidth = 16;
	int64_t segLen;

	parasail_result_verify_rce_source(result);

	asserta(lena > 0);
	asserta(lenb > 0);
	asserta(result->end_query >= 0 && result->end_query < lena);
	asserta(result->end_ref >= 0 && result->end_ref < lenb);

	segLen = (lena + segWidth - 1) / segWidth;
	HT = (int16_t *)result->trace->trace_table;

	size = (size_t)lena + (size_t)lenb;
	cigar = (parasail_cigar_t *)malloc(sizeof(parasail_cigar_t));
	if (!cigar) return NULL;
	cigar->seq = (uint32_t *)malloc(sizeof(uint32_t) * size);
	if (!cigar->seq) {
		free(cigar);
		return NULL;
		}
	cigar->len = 0;
	cigar->beg_query = 0;
	cigar->beg_ref = 0;

	i = result->end_query;
	j = result->end_ref;
	where = PARASAIL_DIAG;

#define INC                                                       \
do {                                                              \
	cigar->len += 1;                                              \
	if ((size_t)cigar->len >= size) {                             \
		size *= 2;                                                  \
		cigar->seq = (uint32_t *) realloc(cigar->seq, sizeof(uint32_t) * size); \
		if (!cigar->seq) {                                          \
			free(cigar);                                            \
			return NULL;                                            \
		}                                                           \
	}                                                               \
} while (0)

#define RESET  \
do {           \
	c_mat = 0; \
	c_del = 0; \
	c_ins = 0; \
} while (0)

#define WRITE(VAL, CHAR)                                          \
do {                                                              \
	INC;                                                          \
	cigar->seq[cigar->len - 1] = parasail_cigar_encode(VAL, CHAR); \
} while (0)

#define WRITE_ANY_RCE      \
do {                         \
	if (c_mat) {             \
		WRITE(c_mat, 'M');   \
	}                        \
	else if (c_del) {        \
		WRITE(c_del, 'D');   \
	}                        \
	else if (c_ins) {        \
		WRITE(c_ins, 'I');   \
	}                        \
	RESET;                   \
} while (0)

	while (i >= 0 || j >= 0) {
		const int64_t loc = j * segLen * segWidth + (i % segLen) * segWidth + (i / segLen);

		if (i < 0) {
			if (c_ins == 0) WRITE_ANY_RCE;
			while (j >= 0) {
				++c_ins;
				--j;
				}
			break;
			}
		if (j < 0) {
			if (c_del == 0) WRITE_ANY_RCE;
			while (i >= 0) {
				++c_del;
				--i;
				}
			break;
			}

		if (where == PARASAIL_DIAG) {
			if (HT[loc] & PARASAIL_DIAG) {
				if (c_mat == 0) WRITE_ANY_RCE;
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
			if (c_ins == 0) WRITE_ANY_RCE;
			c_ins += 1;
			--j;
			if (HT[loc] & PARASAIL_DIAG_E) {
				where = PARASAIL_DIAG;
				}
			else if (HT[loc] & PARASAIL_INS_E) {
				where = PARASAIL_INS;
				}
			else {
				parasail_cigar_free(cigar);
				return NULL;
				}
			}
		else if (where == PARASAIL_DEL) {
			if (c_del == 0) WRITE_ANY_RCE;
			c_del += 1;
			--i;
			if (HT[loc] & PARASAIL_DIAG_F) {
				where = PARASAIL_DIAG;
				}
			else if (HT[loc] & PARASAIL_DEL_F) {
				where = PARASAIL_DEL;
				}
			else {
				parasail_cigar_free(cigar);
				return NULL;
				}
			}
		else if (where == PARASAIL_ZERO) {
			break;
			}
		else {
			parasail_cigar_free(cigar);
			return NULL;
			}
		}

	WRITE_ANY_RCE;

#undef WRITE_ANY_RCE
#undef WRITE
#undef RESET
#undef INC

	cigar_reverse = parasail_reverse_uint32_t(cigar->seq, (size_t)cigar->len);
	free(cigar->seq);
	cigar->seq = cigar_reverse;
	cigar->beg_query = i + 1;
	cigar->beg_ref = j + 1;

	return cigar;
	}