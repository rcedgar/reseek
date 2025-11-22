#include "myutils.h"
#include "parasail.h"

extern parasail_matrix_t parasail_mu_matrix;

void cmd_test_parasail()
	{
	uint LA = 10;
	char *A = myalloc(char, LA);
	for (uint i = 0; i < 10; ++i)
		A[i] = 0;

	const parasail_profile_t *profile =
		parasail_profile_create_avx_256_16(A, LA, &parasail_mu_matrix);

	parasail_result_t* result =
	  parasail_sw_striped_profile_avx2_256_16(profile, A, LA, 2, 1);

	ProgressLog("score=%d\n", result->score);
	}
