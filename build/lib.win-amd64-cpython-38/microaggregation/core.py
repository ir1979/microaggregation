import numpy as np
from typing import Union

try:
    from . import microaggregation_cpp
    _cpp_available = True
    # print("C++ extension available for microaggregation")
except ImportError:
    _cpp_available = False
    import warnings
    warnings.warn("C++ extension not available, falling back to Python implementation")
    # print("If you want to use the C++ implementation, ensure you have built the extension correctly.")

def calculate_sse_dynamic(data: np.ndarray, min_cluster_size: int, use_cpp: bool = True) -> float:
    """
    Calculate Sum of Squared Errors using dynamic programming.
    
    Parameters:
    -----------
    data : np.ndarray
        2D array of shape (n_records, n_features)
    min_cluster_size : int
        Minimum cluster size for microaggregation
    use_cpp : bool, default=True
        Whether to use C++ implementation if available
        
    Returns:
    --------
    float
        Sum of squared errors
    """
    if use_cpp and _cpp_available:
        data_for_cpp = np.ascontiguousarray(data, dtype=np.float64) # Ensure C-contiguity
        return microaggregation_cpp.calculate_sse_dynamic(data_for_cpp.astype(np.float64), min_cluster_size)
    else:
        return _calculate_sse_dynamic_python(data, min_cluster_size)

def _calculate_sse_dynamic_python(data: np.ndarray, min_cluster_size: int) -> float:
    """Python fallback implementation"""
    n_records, n_features = data.shape
    dp_table = np.full((n_records + 1,), float('inf'), dtype=np.float64)
    dp_table[0] = 0.0
    
    cumulative_sum = np.cumsum(data, axis=0)
    cumulative_sum_squared = np.cumsum(data ** 2, axis=0)
    
    for i in range(1, n_records + 1):
        for cluster_size in range(min_cluster_size, min(i + 1, 2 * min_cluster_size)):
            if i - cluster_size >= 0:
                group_sum = cumulative_sum[i - 1] - (cumulative_sum[i - cluster_size - 1] if i - cluster_size > 0 else 0)
                group_sum_squared = cumulative_sum_squared[i - 1] - (cumulative_sum_squared[i - cluster_size - 1] if i - cluster_size > 0 else 0)
                group_sse = np.sum(group_sum_squared - (group_sum ** 2) / cluster_size)
                dp_table[i] = min(dp_table[i], dp_table[i - cluster_size] + group_sse)
    
    return float(dp_table[n_records])

def calculate_total_variance(data: np.ndarray, use_cpp: bool = True) -> float:
    """
    Calculate total variance of the dataset.
    
    Parameters:
    -----------
    data : np.ndarray
        2D array of shape (n_records, n_features)
    use_cpp : bool, default=True
        Whether to use C++ implementation if available
        
    Returns:
    --------
    float
        Total variance
    """
    if use_cpp and _cpp_available:
        data_for_cpp = np.ascontiguousarray(data, dtype=np.float64)
        return microaggregation_cpp.calculate_total_variance(data_for_cpp.astype(np.float64))
    else:
        return float(np.sum((data - np.mean(data, axis=0)) ** 2))

def get_timing():
    """Get timing utilities if C++ extension is available"""
    if _cpp_available:
        return {
            'wall_time': microaggregation_cpp.get_wall_time,
            'cpu_time': microaggregation_cpp.get_cpu_time
        }
    else:
        import time
        return {
            'wall_time': time.time,
            'cpu_time': time.process_time
        }
