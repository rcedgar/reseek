#include "myutils.h"
#include "reseeker.h"
#include "reseeker_thread_body_impl.h"

void reseeker::static_thread_body_nusort(uint threadidx)
	{
	reseeker_thread_body_impl(threadidx, true);
	}
