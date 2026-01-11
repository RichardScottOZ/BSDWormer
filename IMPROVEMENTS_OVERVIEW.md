# BSDWormer - Repository Improvement Summary

## Overview

This repository now contains **two versions** of BSDWormer:

1. **Original (src/)** - The original implementation, preserved unchanged
2. **Version 2 (version2/)** - Modernized implementation with improvements

## Quick Links

### For New Users
👉 **Start here: [version2/README.md](version2/README.md)**

### Version 2 Documentation
- **[version2/INDEX.md](version2/INDEX.md)** - Navigation guide
- **[version2/README.md](version2/README.md)** - Complete documentation
- **[version2/QUICKREF.md](version2/QUICKREF.md)** - Quick reference
- **[version2/IMPROVEMENTS.md](version2/IMPROVEMENTS.md)** - Detailed improvements list
- **[version2/examples/](version2/examples/)** - Usage examples

### Original Code
- **[src/](src/)** - Original implementation (unchanged)
- **[README.md](README.md)** - Original README

## What's New in Version 2?

### 🎯 Key Improvements

1. **Modern Package Structure**
   - Organized modules (core, utils, io, visualization)
   - Installable via pip
   - Type hints throughout

2. **Comprehensive Documentation**
   - Complete README with examples
   - Contribution guidelines
   - Quick reference guide
   - API documentation

3. **Better Code Quality**
   - PEP 8 compliant (snake_case functions)
   - Removed deprecated imports
   - Better error handling
   - Configuration file support

4. **Testing Infrastructure**
   - pytest-based test suite
   - Example tests included
   - CI/CD ready

5. **Easy Installation**
   ```bash
   cd version2
   pip install -e .
   ```

### 📊 Comparison

| Feature | Original (src/) | Version 2 (version2/) |
|---------|----------------|----------------------|
| Package Structure | Scripts | Proper Python package |
| Function Naming | camelCase | snake_case (PEP 8) |
| Type Hints | No | Yes |
| Documentation | Minimal | Comprehensive |
| Configuration | Hardcoded | Config files |
| Testing | Doctests only | pytest suite |
| Installation | Manual | pip install |

## Which Version Should I Use?

### Use Version 2 if you want:
- ✅ Modern Python package
- ✅ Easy installation
- ✅ Better documentation
- ✅ Configuration files
- ✅ Type hints and IDE support
- ✅ Active development

### Use Original if you need:
- ✅ Exact original implementation
- ✅ Existing workflows compatibility
- ✅ Reference implementation

**Recommendation:** Start with Version 2 for new projects.

## Migration Guide

Migrating from original to Version 2 is straightforward:

### Import Changes
```python
# Original
from wormer import Wormer

# Version 2
from bsdwormer import Wormer
```

### Function Names (camelCase → snake_case)
```python
# Original → Version 2
importGdalRaster()  → import_raster()
buildPaddedRaster() → build_padded_raster()
wormLevel()         → worm_level()
```

See [version2/IMPROVEMENTS.md](version2/IMPROVEMENTS.md) for complete details.

## Repository Structure

```
BSDWormer/
├── version2/              ⭐ NEW: Improved implementation
│   ├── bsdwormer/        # Main package
│   ├── tests/            # Test suite
│   ├── examples/         # Usage examples
│   ├── docs/             # Documentation
│   ├── README.md         # Complete documentation
│   ├── IMPROVEMENTS.md   # Improvements summary
│   ├── QUICKREF.md       # Quick reference
│   └── setup.py          # Installation script
│
├── src/                  # Original implementation (unchanged)
│   ├── wormer.py
│   ├── FourierDomainGrid.py
│   ├── FourierDomainOps.py
│   └── ...
│
├── test_data/            # Test datasets
├── README.md             # This file
└── COPYRIGHT.txt         # License
```

## Getting Started

### Option 1: Version 2 (Recommended)

```bash
# Clone repository
git clone https://github.com/RichardScottOZ/BSDWormer.git
cd BSDWormer/version2

# Install
pip install -e .

# Use
python examples/basic_usage.py
```

### Option 2: Original

```bash
# Clone repository
git clone https://github.com/RichardScottOZ/BSDWormer.git
cd BSDWormer/src

# Use existing notebooks
jupyter notebook
```

## Documentation

### Version 2 (Comprehensive)
- [README.md](version2/README.md) - Complete guide
- [INDEX.md](version2/INDEX.md) - Navigation
- [QUICKREF.md](version2/QUICKREF.md) - Quick reference
- [IMPROVEMENTS.md](version2/IMPROVEMENTS.md) - All improvements
- [CONTRIBUTING.md](version2/CONTRIBUTING.md) - How to contribute
- [CHANGELOG.md](version2/CHANGELOG.md) - Version history

### Original
- [README.md](README.md) - Original README
- Jupyter notebooks in `src/`

## Contributing

We welcome contributions! Please see:
- [version2/CONTRIBUTING.md](version2/CONTRIBUTING.md) for guidelines
- Open issues on GitHub for bugs or feature requests

## License

BSD 2-Clause License - See [COPYRIGHT.txt](COPYRIGHT.txt)

## Contact

- **Author:** Franklin G. Horowitz
- **Email:** frank@horow.net
- **Repository:** https://github.com/RichardScottOZ/BSDWormer

## Acknowledgments

- Original BSDWormer implementation: Franklin G. Horowitz
- Version 2 improvements: Modernization while preserving the excellent scientific algorithm

---

## Quick Start Examples

### Version 2 - Simple Usage
```python
from bsdwormer import Wormer

wormer = Wormer()
wormer.import_raster('magnetic_data.tif')
wormer.build_padded_raster((2048, 2048))
wormer.worm_level_as_image(dz=100)
wormer.export_raster(wormer.worm_image, 'output.tif')
```

### Version 2 - With Configuration
```python
from bsdwormer import Wormer, load_config

config = load_config('config.yaml')
wormer = Wormer(config=config)
wormer.import_raster('data.tif')
```

### Version 2 - Multiple Levels
```python
from bsdwormer import Wormer
from bsdwormer.visualization import write_vtk_worm_levels

wormer = Wormer()
wormer.import_raster('data.tif')
wormer.build_padded_raster((2048, 2048))

for dz in [50, 100, 200, 500]:
    wormer.worm_level_as_points(dz)
    wormer.build_worm_segs(dz=dz)
    wormer.build_level_for_vtk(dz)

write_vtk_worm_levels('output', wormer.all_points, 
                      wormer.all_lines, wormer.all_vals)
```

---

**Ready to get started?** Head to [version2/README.md](version2/README.md) for complete documentation!
