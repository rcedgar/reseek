#pragma once

// Search-only helpers for BCAData. Include this (not bcadata.h alone)
// when building struct_data / parasail profiles for search.
// Keeps STRUCTS I/O free of parasail.

#include "bcadata.h"
#include "struct_data.h"
#include "flat_params.h"

struct_data *bca_get_struct_data(
	const BCAData &bca,
	const flat_params &params,
	uint idx,
	chaq_vecs2 *cv,
	uint8_t *scratch_buffer,
	uint scratch_buffer_bytes);

struct_data **bca_get_struct_data_vec(
	BCAData &bca,
	const flat_params &params,
	vector<string> &kept_labels);
