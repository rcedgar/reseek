#pragma once
#include <cstdint>
#include "flat_vec.h"

class aa_seq_i8 : public flat_vec<std::int8_t>
{
private:
    aa_seq_i8(std::size_t n, const char* file, int line)
        : flat_vec<std::int8_t>(ft_aa_seq_i8, n, "aa_seq_i8", file, line)
    {}

public:
    static aa_seq_i8* create(std::size_t n, const char* file, int line)
    {
        return new aa_seq_i8(n, file, line);
    }

    static aa_seq_i8* reference_copy(aa_seq_i8* p)
    {
        if (p) p->add_ref();
        return p;
    }

    static void destroy(aa_seq_i8* p)
    {
        if (p) p->release();
    }
};