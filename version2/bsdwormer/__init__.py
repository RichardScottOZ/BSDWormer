"""
BSDWormer - BSD Licensed Worm Detection for Geophysical Potential Fields

A Python package for detecting multiscale edges ("worms") in geophysical 
potential field data using Poisson wavelet transforms.

Copyright (c) 2013-2017, Franklin G. Horowitz
Licensed under the BSD 2-Clause License.
"""

__version__ = "2.0.0"
__author__ = "Franklin G. Horowitz"
__email__ = "frank@horow.net"
__license__ = "BSD-2-Clause"

# Core functionality
from bsdwormer.core.wormer import Wormer
from bsdwormer.core.fourier_domain_grid import FourierDomainGrid
from bsdwormer.core.fourier_domain_ops import FourierDomainOps

# Configuration
from bsdwormer.config import load_config, Config

__all__ = [
    "Wormer",
    "FourierDomainGrid",
    "FourierDomainOps",
    "load_config",
    "Config",
    "__version__",
]
