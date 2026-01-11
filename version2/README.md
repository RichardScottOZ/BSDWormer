# BSDWormer Version 2 - Improvements

This is an improved, modernized version of the BSDWormer package. **This version does not modify the original code** which remains in the `src/` directory.

## What is BSDWormer?

BSDWormer is a BSD-licensed implementation of the Poisson wavelet multiscale edge ("worm") detection algorithm for geophysical potential fields. It uses Fourier domain operations to detect edges in gravity and magnetic field data.

## Key Improvements in Version 2

### 1. **Proper Package Structure**
- Organized as a proper Python package with clear module hierarchy
- Installable via pip with `setup.py`
- Clear separation of concerns: core algorithms, utilities, I/O, and visualization

### 2. **Enhanced Documentation**
- Comprehensive README with examples
- Improved docstrings with type hints
- API documentation
- Tutorial examples

### 3. **Modern Python Practices**
- Removed deprecated `future` imports (Python 3+ native)
- Type hints for better IDE support and code clarity
- Proper error handling with custom exceptions
- Configuration file support

### 4. **Better Code Quality**
- Removed hardcoded paths
- Extracted magic numbers into constants
- Improved error messages
- Added logging support

### 5. **Testing Infrastructure**
- Unit tests using pytest
- Extracted doctests into proper test files
- Test fixtures for reproducibility
- CI/CD ready

### 6. **Dependency Management**
- Clear requirements.txt with pinned versions
- Separated dev dependencies
- Documentation of system dependencies

## Installation

### System Dependencies

BSDWormer requires GDAL/OGR for raster I/O:

```bash
# Ubuntu/Debian
sudo apt-get install gdal-bin python3-gdal

# macOS with Homebrew
brew install gdal

# Windows
# Use OSGeo4W or install from https://www.gisinternals.com/
```

### Python Package Installation

```bash
cd version2
pip install -e .  # Development install
# or
pip install .     # Regular install
```

## Quick Start

### Basic Usage

```python
from bsdwormer import Wormer
import numpy as np

# Create a wormer instance
wormer = Wormer()

# Load a raster file
wormer.import_raster('path/to/magnetic_data.tif')

# Process at a specific depth/height level
dz = 100  # meters
wormer.build_padded_raster(padded_shape=(2048, 2048))
wormer.worm_level_as_image(dz)

# Export results
wormer.export_raster(wormer.worm_image, 'worms_100m.tif')
```

### Processing Multiple Levels

```python
from bsdwormer import Wormer

wormer = Wormer()
wormer.import_raster('magnetic_data.tif')
wormer.build_padded_raster(padded_shape=(2048, 2048))

# Process multiple depth levels
depths = [50, 100, 200, 500, 1000]  # meters
for dz in depths:
    wormer.worm_level_as_points(dz)
    wormer.build_worm_segs(dz=dz)
    wormer.build_level_for_vtk(dz)

# Export to VTK for visualization
from bsdwormer.utils import write_vtk_worm_levels
write_vtk_worm_levels(
    'output_worms',
    wormer.all_points,
    wormer.all_lines,
    wormer.all_vals
)
```

## Module Organization

```
version2/
├── bsdwormer/           # Main package
│   ├── __init__.py      # Package initialization
│   ├── core/            # Core algorithms
│   │   ├── __init__.py
│   │   ├── wormer.py    # Main Wormer class
│   │   ├── fourier_domain_grid.py
│   │   └── fourier_domain_ops.py
│   ├── utils/           # Utility functions
│   │   ├── __init__.py
│   │   ├── fft_utils.py
│   │   └── geometry.py
│   ├── io/              # Input/Output operations
│   │   ├── __init__.py
│   │   ├── raster.py    # GDAL raster I/O
│   │   └── database.py  # PostGIS database I/O
│   └── visualization/   # Visualization utilities
│       ├── __init__.py
│       └── vtk_writer.py
├── tests/               # Test suite
│   ├── test_wormer.py
│   ├── test_fourier.py
│   └── fixtures/
├── examples/            # Example scripts and notebooks
│   ├── basic_usage.py
│   └── notebooks/
├── docs/                # Documentation
└── setup.py             # Package installation
```

## Configuration

BSDWormer v2 supports configuration files for common parameters:

```yaml
# config.yaml
padding:
  rolloff_size: 100
  pad_type: 'hann'

processing:
  nodata_value: -100
  log_vals: true

output:
  format: 'GeoTIFF'
  compression: 'LZW'
```

Load configuration:

```python
from bsdwormer import Wormer, load_config

config = load_config('config.yaml')
wormer = Wormer(config=config)
```

## Examples

See the `examples/` directory for:
- Basic usage examples
- Jupyter notebooks with visualizations
- Processing workflows for different data types

## Development

### Running Tests

```bash
pytest tests/
```

### Code Style

We use:
- `black` for code formatting
- `flake8` for linting
- `mypy` for type checking

```bash
black bsdwormer/
flake8 bsdwormer/
mypy bsdwormer/
```

## Citation

If you use BSDWormer in your research, please cite:

```
@software{bsdwormer,
  author = {Horowitz, Franklin G.},
  title = {BSDWormer: Poisson Wavelet Multiscale Edge Detection for Geophysical Potential Fields},
  year = {2013-2017},
  license = {BSD-2-Clause}
}
```

## License

BSDWormer is licensed under the BSD 2-Clause License. See the COPYRIGHT.txt file in the root directory for the full license text.

## Contributing

Contributions are welcome! Please see CONTRIBUTING.md for guidelines.

## Support

For questions or issues:
- Open an issue on GitHub
- Contact: Frank Horowitz <frank@horow.net>

## References

The "worm" algorithm is based on multiscale edge detection using Poisson wavelet transforms. Key references:

1. Poisson wavelet analysis for potential field data
2. Multiscale edge detection in geophysical potential fields
3. Fourier domain upward continuation operators

## Changelog

See CHANGELOG.md for version history and changes.
