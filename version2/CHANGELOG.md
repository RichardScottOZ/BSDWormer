# Changelog

All notable changes to BSDWormer will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2024-01-11

### Added
- **Proper package structure** with organized modules (core, utils, io, visualization)
- **Setup.py** for pip installation
- **Configuration file support** (YAML/JSON)
- **Type hints** throughout the codebase
- **Comprehensive documentation** including README, CONTRIBUTING, and API docs
- **Modern Python 3.7+ support** with removal of deprecated `future` imports
- **Logging support** for better debugging
- **Unit tests** extracted from doctests with pytest framework
- **CI/CD configuration** templates
- **Example scripts** and notebooks in examples/ directory

### Changed
- **Renamed functions** to follow PEP 8 conventions (snake_case)
  - `setBaseGrid()` → `set_base_grid()`
  - `importGdalRaster()` → `import_raster()`
  - `buildPaddedRaster()` → `build_padded_raster()`
  - `wormLevel()` → `worm_level()`
  - `buildWavenumbers()` → `build_wavenumbers()`
  - And many more...
- **Improved error handling** with custom exception classes
- **Better default values** and configuration management
- **Module organization** for better code maintainability

### Removed
- **Hardcoded paths** (e.g., in Utility.py)
- **Deprecated `future` imports** (now Python 3 native)
- **Commented-out code** and excessive FIXME comments
- **Magic numbers** (now constants or configurable)

### Fixed
- **Import issues** with proper package structure
- **Deprecated numpy types** usage
- **Path handling** to be platform-independent

### Security
- **Input validation** for all file operations
- **Proper exception handling** to prevent information leakage

### Documentation
- Complete rewrite of README with examples
- Added CONTRIBUTING.md with development guidelines
- Added CHANGELOG.md (this file)
- Improved docstrings with type hints and examples
- Added inline comments for complex algorithms

### Migration from Version 1

If you're migrating from the original BSDWormer (in `src/`), here are the key changes:

#### Import Changes
```python
# Old (version 1)
from wormer import Wormer

# New (version 2)
from bsdwormer import Wormer
```

#### Function Name Changes
```python
# Old (version 1)
wormer.importGdalRaster('file.tif')
wormer.buildPaddedRaster((2048, 2048))
wormer.wormLevel(100)

# New (version 2)
wormer.import_raster('file.tif')
wormer.build_padded_raster((2048, 2048))
wormer.worm_level(100)
```

#### Configuration
```python
# New in version 2
from bsdwormer import load_config

config = load_config('config.yaml')
wormer = Wormer(config=config)
```

## [1.0.0] - 2017-XX-XX (Original)

### Original Features
- Core worm detection algorithm
- Fourier domain operations
- GDAL raster I/O
- VTK output for visualization
- PostGIS database integration
- Jupyter notebook examples

---

## Version Numbering

- **Major version** (X.0.0): Incompatible API changes
- **Minor version** (0.X.0): New features, backward compatible
- **Patch version** (0.0.X): Bug fixes, backward compatible
