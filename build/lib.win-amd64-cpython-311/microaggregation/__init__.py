"""
Microaggregation: Fast algorithms for data privacy protection through microaggregation
"""

__version__ = "0.1.7"
__author__ = "Reza Mortazavi"
__email__ = "ir1979@gmail.com"

from .core import calculate_sse_dynamic, calculate_total_variance, get_timing

__all__ = ['calculate_sse_dynamic', 'calculate_total_variance', 'get_timing']
