#pragma once
#include "flat_base.h"

template <typename T>
class flat_vec : public flat_base<T>
{
protected:
    flat_vec(flat_type ft, std::size_t n,
             const char* type_name,
             const char* file, int line)
        : flat_base<T>(ft, n, type_name, file, line)
    {}

public:
    // Convenience indexing (same as get/set)
    T operator[](std::size_t i) const { return this->get(i); }
    T& at(std::size_t i) { assert(i < this->m_size); return this->m_data[i]; }
};
