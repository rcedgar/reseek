#include <vector>
#include <algorithm>
#include <cstdint>
#include <cassert>

double roc_auc_binary(const std::vector<float>& scores,
                      const std::vector<bool>& is_tp)
{
    const size_t N = scores.size();
    assert(is_tp.size() == N);
    assert(N > 0);

    size_t n_tp = 0;
    for (size_t i = 0; i < N; ++i)
        n_tp += (is_tp[i] ? 1 : 0);

    const size_t n_fp = N - n_tp;

    assert(n_tp > 0);
    assert(n_fp > 0);

    std::vector<size_t> order(N);
    for (size_t i = 0; i < N; ++i)
        order[i] = i;

    std::sort(order.begin(), order.end(),
        [&](size_t a, size_t b)
        {
            return scores[a] < scores[b];
        });

    double tp_rank_sum = 0.0;

    size_t i = 0;
    while (i < N)
    {
        size_t j = i + 1;
        const float s = scores[order[i]];
        while (j < N && scores[order[j]] == s)
            ++j;

        // Tied block [i, j), 0-based positions in sorted order.
        // Average 1-based rank = ((i+1) + j) / 2.
        const double avg_rank = 0.5 * double(i + 1 + j);

        size_t tp_in_tie = 0;
        for (size_t k = i; k < j; ++k)
            tp_in_tie += (is_tp[order[k]] ? 1 : 0);

        tp_rank_sum += avg_rank * double(tp_in_tie);
        i = j;
    }

    const double npos = double(n_tp);
    const double nneg = double(n_fp);

    const double auc =
        (tp_rank_sum - npos * (npos + 1.0) * 0.5) / (npos * nneg);

    return auc;
}