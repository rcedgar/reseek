#pragma once
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cassert>
#include <type_traits>

#include "flat_alloc.h"
#include "flat_stats.h"
#include "flat_round.h"

// Alignment policy: adjust as desired (e.g. 32 or 64 for AVX512 friendliness).
static const std::size_t FLAT_ALIGNMENT = 32;

template <typename T>
class flat_base
{
public:
    // Trivial types only (as you want).
    // If you need pre-C++11, remove static_assert and rely on discipline.
    static_assert(std::is_trivial<T>::value, "flat_base requires trivial T");

    std::atomic<uint32_t> m_RefCount;  // intrusive
    flat_type m_ft;

    T* m_data;
    std::size_t m_size;   // requested elements
    std::size_t m_cap;    // allocated capacity (rounded)

protected:
    // Protected: leaf types control construction via create().
    flat_base(flat_type ft, std::size_t n,
              const char* type_name,
              const char* file, int line)
        : m_RefCount(1), m_ft(ft), m_data(0), m_size(n), m_cap(0)
    {
        m_cap = round_capacity(m_size);

        const std::size_t bytes = m_cap * sizeof(T);
        void* p = flatmem::aligned_malloc(bytes, FLAT_ALIGNMENT);
        assert(p && "aligned_malloc failed");

        m_data = (T*)p;

        flatstats::on_create(m_ft, bytes, type_name, file, line, this);
    }

    virtual ~flat_base()
    {
        // Called only when refcount hits zero.
        const std::size_t bytes = m_cap * sizeof(T);

        flatstats::on_delete(m_ft, bytes, this);
        flatmem::aligned_free(m_data);

        m_data = 0;
        m_size = m_cap = 0;
    }

public:
    // No copy (prevents double-free, counter corruption)
    flat_base(const flat_base&);
    flat_base& operator=(const flat_base&);

    void add_ref()
    {
        m_RefCount.fetch_add(1, std::memory_order_relaxed);
    }

    void release()
    {
        // Underflow check: if it was already 0, that's a bug.
        uint32_t prev = m_RefCount.fetch_sub(1, std::memory_order_acq_rel);
        assert(prev != 0 && "refcount underflow (destroy called too many times)");

        if (prev == 1)
            delete this;
    }

    std::size_t size() const { return m_size; }
    std::size_t capacity() const { return m_cap; }

    // Fast unchecked accessors (you can add checked versions if desired).
    T get(std::size_t i) const
    {
        assert(i < m_size);
        return m_data[i];
    }

    void set(std::size_t i, T v)
    {
        assert(i < m_size);
        m_data[i] = v;
    }

    T* data() { return m_data; }
    const T* data() const { return m_data; }
};