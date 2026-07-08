#include "myutils.h"
#include <numeric>

/**
 * @brief Calculates the arithmetic mean (average) of a vector.
 *
 * @tparam T The numeric type of the vector elements (e.g., int, float, double).
 * @param v The input vector.
 * @return T The calculated mean.
 */
template <typename T>
T calculate_mean(const vector<T>& v) {
    if (v.empty()) {
        throw invalid_argument("Vector cannot be empty for mean calculation.");
    }
    // accumulate calculates the sum, starting with an initial value of 0.0
    // We explicitly cast the initial value to double to ensure floating-point division
    // and avoid potential overflow for large integer sums, before casting back to T if needed.
    double sum = accumulate(v.begin(), v.end(), 0.0);
    return static_cast<T>(sum / v.size());
}

/**
 * @brief Calculates Pearson's correlation coefficient (r) between two vectors.
 *
 * Pearson's r measures the linear correlation between two sets of data. The result is a value
 * between -1 and +1 inclusive, where +1 is total positive linear correlation, 0 is no linear
 * correlation, and -1 is total negative linear correlation.
 *
 * Formula: r = Sum[(xi - mean_x) * (yi - mean_y)] / (sqrt[Sum(xi - mean_x)^2] * sqrt[Sum(yi - mean_y)^2])
 *
 * @tparam T The numeric type of the vector elements (e.g., int, float, double).
 * @param x The first input vector.
 * @param y The second input vector.
 * @return double The Pearson's correlation coefficient.
 * @throws invalid_argument If the vectors are of unequal size or are empty.
 */
template <typename T>
double pearson_correlation(const vector<T>& x, const vector<T>& y) {
    // 1. Input Validation
    if (x.size() != y.size()) {
        throw invalid_argument("Vectors must be of the same size.");
    }
    size_t n = x.size();
    if (n == 0) {
        throw invalid_argument("Vectors cannot be empty.");
    }

    // 2. Calculate Means
    double mean_x = calculate_mean(x);
    double mean_y = calculate_mean(y);

    // 3. Calculate Components
    // Initialize components for the numerator and the two parts of the denominator
    double numerator_sum = 0.0;           // Sum[(xi - mean_x) * (yi - mean_y)] (Covariance component)
    double x_std_dev_sum_sq = 0.0;        // Sum[(xi - mean_x)^2]
    double y_std_dev_sum_sq = 0.0;        // Sum[(yi - mean_y)^2]

    for (size_t i = 0; i < n; ++i) {
        // Calculate differences from the mean
        double diff_x = static_cast<double>(x[i]) - mean_x;
        double diff_y = static_cast<double>(y[i]) - mean_y;

        // Numerator component (cross-product)
        numerator_sum += (diff_x * diff_y);

        // Denominator components (squared differences)
        x_std_dev_sum_sq += (diff_x * diff_x);
        y_std_dev_sum_sq += (diff_y * diff_y);
    }

    // 4. Calculate Denominator (product of standard deviations)
    double denominator = sqrt(x_std_dev_sum_sq) * sqrt(y_std_dev_sum_sq);

    // 5. Calculate Correlation Coefficient
    if (denominator == 0.0) {
        // This occurs if one or both vectors have zero variance (all elements are the same).
        // In this edge case, correlation is mathematically undefined, but statistically
        // often treated as 0 or 1 depending on context. For a robust result, we return 0.
        return 0.0;
    }

    return numerator_sum / denominator;
}