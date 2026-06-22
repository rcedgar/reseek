#pragma once

#include "hitdata.h"

void reseek_hit_sink_begin(uint thread_count);
void reseek_hit_sink_submit(const reseek_hit &hit, bool nu_only);
void reseek_hit_sink_thread_end(uint threadidx);
void reseek_hit_sink_flush(FILE *fhit);
