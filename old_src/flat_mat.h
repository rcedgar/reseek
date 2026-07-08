#pragma once
#include "flat_base.h"

template <typename T>
class flat_mat : public flat_base<T>
{
protected:
    std::size_t m_rows;
    std::size_t m_cols;

    flat_mat(flat_type ft, std::size_t rows, std::size_t cols,
             const char* type_name,
             const char* file, int line)
        : flat_base<T>(ft, rows * cols, type_name, file, line),
          m_rows(rows), m_cols(cols)
    {}

public:
    std::size_t rows() const { return m_rows; }
    std::size_t cols() const { return m_cols; }

    T get(std::size_t r, std::size_t c) const
    {
        assert(r < m_rows && c < m_cols);
        return this->m_data[r * m_cols + c];
    }

    void set(std::size_t r, std::size_t c, T v)
    {
        assert(r < m_rows && c < m_cols);
        this->m_data[r * m_cols + c] = v;
    }

    T& at(std::size_t r, std::size_t c)
    {
        assert(r < m_rows && c < m_cols);
        return this->m_data[r * m_cols + c];
    }
};