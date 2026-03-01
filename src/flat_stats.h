#pragma once
#include <cstddef>
#include <cstdint>

enum flat_type
{
    ft_aa_seq_i8,
    ft_dist_i8,
    ft_matrix_f32,
    ft_n
};

namespace flatstats
{
    void on_create(flat_type ft, std::size_t bytes,
                   const char* type_name,
                   const char* file, int line,
                   void* obj_ptr);

    void on_delete(flat_type ft, std::size_t bytes, void* obj_ptr);

    // Snapshot counters (thread-safe reads).
    uint64_t creates(flat_type ft);
    uint64_t deletes(flat_type ft);
    uint64_t bytes(flat_type ft);

    // Leak report (call single-threaded after workers complete).
    // Returns number of currently-live objects tracked.
    std::size_t report_leaks();
}