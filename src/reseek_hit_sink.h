#pragma once

#include "hitdata.h"

void reseek_hit_sink_begin();
void reseek_hit_sink_submit(const reseek_hit &hit, bool nu_only);
void reseek_hit_sink_thread_end();
void reseek_hit_sink_flush(FILE *fhit);
