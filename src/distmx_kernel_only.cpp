#if 0
#include "distmx_kernel.h"
extern uint16_t *gx;
extern uint16_t *gy;
extern uint16_t *gz;
extern uint32_t *gdist;

void distmx_kernel_only()
{
    dist_fill_avx2(gx, gy, gz, 100, gdist);
}
#endif // 0