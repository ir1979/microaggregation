#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include "myTiming.h"

namespace py = pybind11;

double calculate_sse_dynamic_cpp(py::array_t<double> input_data, int min_cluster_size) {
    py::buffer_info buf = input_data.request();
    
    if (buf.ndim != 2) {
        throw std::runtime_error("Input array must be 2-dimensional");
    }
    
    int n_records = buf.shape[0];
    int n_features = buf.shape[1];
    
    double* data = static_cast<double*>(buf.ptr);
    
    // Initialize DP table
    std::vector<double> dp_table(n_records + 1, std::numeric_limits<double>::infinity());
    dp_table[0] = 0.0;
    
    // Pre-compute cumulative sums for efficiency
    std::vector<std::vector<double>> cumulative_sum(n_records, std::vector<double>(n_features, 0.0));
    std::vector<std::vector<double>> cumulative_sum_squared(n_records, std::vector<double>(n_features, 0.0));
    
    // Fill cumulative sums
    for (int i = 0; i < n_records; i++) {
        for (int j = 0; j < n_features; j++) {
            double val = data[i * n_features + j];
            cumulative_sum[i][j] = val + (i > 0 ? cumulative_sum[i-1][j] : 0.0);
            cumulative_sum_squared[i][j] = val * val + (i > 0 ? cumulative_sum_squared[i-1][j] : 0.0);
        }
    }
    
    // Dynamic programming
    for (int i = 1; i <= n_records; i++) {
        int max_cluster_size = std::min(i, 2 * min_cluster_size-1);
        // int loop_exclusive_upper_bound = std::min(i + 1, 2 * min_cluster_size);

        for (int cluster_size = min_cluster_size; cluster_size <= max_cluster_size; cluster_size++) {
        // for (int cluster_size = min_cluster_size; cluster_size < loop_exclusive_upper_bound; cluster_size++) {
            if (i - cluster_size >= 0) {
                double group_sse = 0.0;
                
                // Calculate SSE for this group
                for (int j = 0; j < n_features; j++) {
                    double group_sum = cumulative_sum[i-1][j];
                    double group_sum_sq = cumulative_sum_squared[i-1][j];
                    
                    if (i - cluster_size > 0) {
                        group_sum -= cumulative_sum[i - cluster_size - 1][j];
                        group_sum_sq -= cumulative_sum_squared[i - cluster_size - 1][j];
                    }
                    
                    group_sse += group_sum_sq - (group_sum * group_sum) / cluster_size;
                }
                
                dp_table[i] = std::min(dp_table[i], dp_table[i - cluster_size] + group_sse);
            }
        }
    }
    
    return dp_table[n_records];
}

double calculate_total_variance_cpp(py::array_t<double> input_data) {
    py::buffer_info buf = input_data.request();
    
    if (buf.ndim != 2) {
        throw std::runtime_error("Input array must be 2-dimensional");
    }
    
    int n_records = buf.shape[0];
    int n_features = buf.shape[1];

    //std::cout << "Number of records: " << n_records << ", Number of features: " << n_features << std::endl;
    
    double* data = static_cast<double*>(buf.ptr);
    
    // Calculate means
    std::vector<double> means(n_features, 0.0);
    for (int i = 0; i < n_records; i++) {
        for (int j = 0; j < n_features; j++) {
            means[j] += data[i * n_features + j];
            // if (i<3) { // Only print first 3 records for debugging
            //     std::cout << "Record " << i << ", Feature " << j << ": " << std::fixed << std::setprecision(8) 
            //               << data[i * n_features + j] << std::endl;
            // }
        }
    }
    for (int j = 0; j < n_features; j++) {
        means[j] /= n_records;
        // std::cout << "Mean of feature " << j << ": " << means[j] << std::endl;
    }
    
    // Calculate total variance
    double total_variance = 0.0;
    for (int i = 0; i < n_records; i++) {
        for (int j = 0; j < n_features; j++) {
            double diff = data[i * n_features + j] - means[j];
            total_variance += diff * diff;
        }
    }
    
    return total_variance;
}

// Timing utilities
double get_wall_time_wrapper() {
    return get_wall_time();
}

double get_cpu_time_wrapper() {
    return get_cpu_time();
}

PYBIND11_MODULE(microaggregation_cpp, m) {
    m.doc() = "Fast C++ implementation of microaggregation algorithms";
    
    m.def("calculate_sse_dynamic", &calculate_sse_dynamic_cpp, 
          "Calculate Sum of Squared Errors using dynamic programming",
          py::arg("data"), py::arg("min_cluster_size"));
    
    m.def("calculate_total_variance", &calculate_total_variance_cpp,
          "Calculate total variance of the dataset",
          py::arg("data"));
          
    m.def("get_wall_time", &get_wall_time_wrapper, "Get wall clock time");
    m.def("get_cpu_time", &get_cpu_time_wrapper, "Get CPU time");
}
