# Microaggregation

Fast C++ implementation of microaggregation algorithms for data privacy protection.

## Features

- High-performance C++ implementation of SSE dynamic programming
- Python fallback when C++ extension is not available
- Compatible with NumPy arrays
- Timing utilities for performance measurement

## Installation

bash pip install microaggregation

## Usage

```python
import numpy as np
from microaggregation import calculate_sse_dynamic, calculate_total_variance, get_timing

# Example usage:
my_data = np.random.rand(100, 5)
min_k = 3

sse_value = calculate_sse_dynamic(my_data, min_k)
# sse_value_py_only = calculate_sse_dynamic(my_data, min_k, use_cpp=False) # To force Python

total_var = calculate_total_variance(my_data)

timers = get_timing()
start_wall = timers['wall_time']()
# ... some operation ...
end_wall = timers['wall_time']()
print(f"Wall time: {end_wall - start_wall}")

print(f"SSE: {sse_value}, Total Variance: {total_var}")
print(f"IL: {sse_value / total_var * 100 if total_var != 0 else float('inf')}")
```


