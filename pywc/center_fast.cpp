#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

/*
 * This module provides a fast C++ implementation of the centering algorithm
 * used in pywc's Interface._center() method. The bottleneck is the repeated
 * histogram calculation in a tight loop.
 *
 * The centering algorithm iteratively shifts the system until the density
 * at the box boundaries falls below a threshold. This C++ implementation
 * provides ~10-20x speedup over the pure Python version.
 *
 * OpenMP parallelization is supported for the histogram computation.
 */

#ifdef _OPENMP
#include <omp.h>
#endif

namespace py = pybind11;

namespace
{

    /**
     * Fast 1D histogram computation.
     *
     * @param data Array of positions along one dimension
     * @param n_data Number of data points
     * @param bins Number of histogram bins
     * @param range_min Minimum of histogram range
     * @param range_max Maximum of histogram range
     * @param density Output density array (must be pre-allocated with size bins)
     */
    void compute_histogram_1d(const double *data, std::size_t n_data,
                              int bins, double range_min, double range_max,
                              double *density)
    {
        const double bin_width = (range_max - range_min) / bins;
        const double inv_bin_width = 1.0 / bin_width;

        // Zero out density array
        std::fill_n(density, bins, 0.0);

        // Count occurrences in each bin
        std::vector<std::size_t> counts(bins, 0);

#pragma omp parallel
        {
            // Thread-local counts to avoid race conditions
            std::vector<std::size_t> local_counts(bins, 0);

#pragma omp for nowait
            for (std::size_t i = 0; i < n_data; ++i)
            {
                const double val = data[i];
                if (val >= range_min && val < range_max)
                {
                    const int bin_idx = static_cast<int>((val - range_min) * inv_bin_width);
                    if (bin_idx >= 0 && bin_idx < bins)
                    {
                        local_counts[bin_idx]++;
                    }
                }
            }

// Reduction: combine local counts into global counts
#pragma omp critical
            {
                for (int i = 0; i < bins; ++i)
                {
                    counts[i] += local_counts[i];
                }
            }
        }

        // Normalize to density
        std::size_t total_count = 0;
        for (int i = 0; i < bins; ++i)
        {
            total_count += counts[i];
        }

        if (total_count > 0)
        {
            const double norm = 1.0 / (static_cast<double>(total_count) * bin_width);
            for (int i = 0; i < bins; ++i)
            {
                density[i] = static_cast<double>(counts[i]) * norm;
            }
        }
    }

    /**
     * Fast centering algorithm.
     *
     * Iteratively shifts positions until boundary densities fall below threshold.
     *
     * @param pos_group Positions along centering direction (modified in-place)
     * @param range_min Minimum of histogram range
     * @param range_max Maximum of histogram range
     * @param delta Density threshold for convergence
     * @param shift Amount to shift per iteration
     * @param max_shift Maximum total shift allowed
     * @param bins Number of histogram bins (default: 10)
     * @param max_iter Maximum iterations (default: 100)
     * @return Total shift applied
     */
    double center_fast(py::array_t<double> &pos_group,
                       double range_min,
                       double range_max,
                       double delta,
                       double shift,
                       double max_shift,
                       int bins = 10,
                       int max_iter = 100)
    {

        // Validate input
        if (pos_group.ndim() != 1)
        {
            throw std::invalid_argument("pos_group must be 1D array");
        }
        if (bins <= 0 || bins > 1000)
        {
            throw std::invalid_argument("bins must be in range (0, 1000]");
        }
        if (shift <= 0.0)
        {
            throw std::invalid_argument("shift must be positive");
        }
        if (max_iter <= 0)
        {
            throw std::invalid_argument("max_iter must be positive");
        }
        auto buf = pos_group.request();
        auto *data = static_cast<double *>(buf.ptr);
        const std::size_t n_data = static_cast<std::size_t>(buf.shape[0]);

        std::vector<double> density(bins);
        double total_shift = 0.0;

        // Release GIL for the computation
        {
            py::gil_scoped_release release;

            for (int iteration = 0; iteration < max_iter; ++iteration)
            {
                // Compute histogram
                compute_histogram_1d(data, n_data, bins, range_min, range_max, density.data());

                // Check convergence: boundary densities below threshold
                if (!(density[0] > delta || density[bins - 1] > delta))
                {
                    break;
                }

                // Apply shift
                total_shift += shift;
                if (total_shift >= max_shift)
                {
                    throw std::runtime_error("Centering failure: maximum shift exceeded");
                }

// Shift all positions in-place
#pragma omp parallel for
                for (std::size_t i = 0; i < n_data; ++i)
                {
                    data[i] += shift;
                }
            }
        }

        return total_shift;
    }

    /**
     * Compute mean of positions array.
     *
     * @param pos_group 1D array of positions
     * @return Mean value
     */
    double compute_mean(const py::array_t<double> &pos_group)
    {
        auto buf = pos_group.request();
        const auto *data = static_cast<const double *>(buf.ptr);
        const std::size_t n_data = static_cast<std::size_t>(buf.shape[0]);

        if (n_data == 0)
        {
            return 0.0;
        }

        double sum = 0.0;

        {
            py::gil_scoped_release release;

#pragma omp parallel for reduction(+ : sum)
            for (std::size_t i = 0; i < n_data; ++i)
            {
                sum += data[i];
            }
        }

        return sum / static_cast<double>(n_data);
    }

    /**
     * Apply a shift to all elements of an array in-place.
     *
     * @param pos_array 1D array of positions
     * @param shift Value to add to all elements
     */
    void apply_shift(py::array_t<double> &pos_array, double shift)
    {
        auto buf = pos_array.request();
        auto *data = static_cast<double *>(buf.ptr);
        const std::size_t n_data = static_cast<std::size_t>(buf.shape[0]);

        {
            py::gil_scoped_release release;

#pragma omp parallel for
            for (std::size_t i = 0; i < n_data; ++i)
            {
                data[i] += shift;
            }
        }
    }

} // namespace

PYBIND11_MODULE(center_fast, m)
{
    m.doc() = "Fast C++ implementation of pywc's centering algorithm";

    m.def("center_fast", &center_fast,
          py::arg("pos_group"),
          py::arg("range_min"),
          py::arg("range_max"),
          py::arg("delta"),
          py::arg("shift"),
          py::arg("max_shift"),
          py::arg("bins") = 10,
          py::arg("max_iter") = 100,
          R"pbdoc(
          Fast centering of positions along one dimension.

          Iteratively shifts positions until the density at the boundaries
          falls below a threshold. This is the performance-critical part of
          Interface._center().

          Parameters
          ----------
          pos_group : ndarray (N,)
              Positions along centering direction (modified in-place)
          range_min : float
              Minimum of histogram range
          range_max : float
              Maximum of histogram range
          delta : float
              Density threshold for convergence
          shift : float
              Amount to shift per iteration
          max_shift : float
              Maximum total shift allowed (raises RuntimeError if exceeded)
          bins : int, optional
              Number of histogram bins (default: 10)
          max_iter : int, optional
              Maximum iterations (default: 100)

          Returns
          -------
          total_shift : float
              Total shift applied to positions

          Notes
          -----
          - The pos_group array is modified in-place
          - Uses OpenMP for parallelization if available
          - Releases the GIL during computation
          )pbdoc");

    m.def("compute_mean", &compute_mean,
          py::arg("pos_group"),
          "Compute mean of 1D position array (parallel, GIL-free)");

    m.def("apply_shift", &apply_shift,
          py::arg("pos_array"),
          py::arg("shift"),
          "Apply shift to all elements in-place (parallel, GIL-free)");
}
