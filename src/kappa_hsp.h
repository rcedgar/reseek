#pragma once

#include "myutils.h"

extern int16_t kappa32_flat_logodds[32 * 32];

void kappa_get_hsp_limits(int LQ, int LT, int Diag,
	int &mini, int &minj, int &n);

int kappa_find_hsp(const byte *QSeq, const byte *TSeq,
	int LQ, int LT, int Diag);

int kappa_max_pos_logodds();
