#pragma once
#include <cstdint>
#include "flat_vec.h"

class dist_list_i8 : public flat_vec<std::int8_t>
{
private:
    dist_list_i8(std::size_t n, const char* file, int line)
        : flat_vec<std::int8_t>(ft_dist_i8, n, "dist_list_i8", file, line)
    {}

public:
    static dist_list_i8* create(std::size_t n, const char* file, int line)
    {
        return new dist_list_i8(n, file, line);
    }

    static dist_list_i8* reference_copy(dist_list_i8* p)
    {
        if (p) p->add_ref();
        return p;
    }

    static void destroy(dist_list_i8* p)
    {
        if (p) p->release();
    }
};
