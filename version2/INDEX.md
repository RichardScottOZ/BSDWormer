# Welcome to BSDWormer Version 2!

## 📋 Quick Navigation

- **[README.md](README.md)** - Complete overview, installation, and usage
- **[IMPROVEMENTS.md](IMPROVEMENTS.md)** - Detailed list of all improvements  
- **[QUICKREF.md](QUICKREF.md)** - Quick reference guide
- **[CONTRIBUTING.md](CONTRIBUTING.md)** - How to contribute
- **[CHANGELOG.md](CHANGELOG.md)** - Version history

## 🎯 What's Here

### Documentation
- `README.md` - Main documentation with examples
- `IMPROVEMENTS.md` - Summary of improvements over version 1
- `QUICKREF.md` - Quick reference for common tasks
- `CONTRIBUTING.md` - Contribution guidelines
- `CHANGELOG.md` - Version history and changes

### Code
- `bsdwormer/` - Main package with organized modules
  - `core/` - Core algorithms (Wormer, FourierDomainGrid, FourierDomainOps)
  - `utils/` - Utility functions (FFT utils, geometry)
  - `io/` - Input/output operations (raster, database)
  - `visualization/` - Visualization tools (VTK writer)
- `tests/` - Unit test suite
- `examples/` - Example scripts and notebooks

### Configuration
- `config.example.yaml` - Example configuration file
- `setup.py` - Package installation script
- `requirements.txt` - Python dependencies
- `requirements-dev.txt` - Development dependencies
- `pytest.ini` - Test configuration

## 🚀 Quick Start

1. **Install:**
   ```bash
   cd version2
   pip install -e .
   ```

2. **Use:**
   ```python
   from bsdwormer import Wormer
   
   wormer = Wormer()
   wormer.import_raster('data.tif')
   wormer.build_padded_raster((2048, 2048))
   wormer.worm_level_as_image(dz=100)
   wormer.export_raster(wormer.worm_image, 'output.tif')
   ```

3. **Learn More:**
   - See `examples/basic_usage.py` for a complete workflow
   - Read `README.md` for detailed documentation
   - Check `QUICKREF.md` for common patterns

## 🔑 Key Features

### ✨ Modern Python Package
- Proper package structure with `__init__.py` files
- Installable via pip
- Type hints for better IDE support
- PEP 8 compliant naming

### 📚 Comprehensive Documentation
- Detailed README with examples
- Contribution guidelines
- Quick reference guide
- Inline code documentation

### 🧪 Testing Infrastructure
- Unit tests with pytest
- Test fixtures
- Coverage reporting
- CI/CD ready

### ⚙️ Configuration Support
- YAML/JSON configuration files
- Centralized parameter management
- Easy parameter tuning
- Reproducible workflows

### 🛠️ Development Tools
- Code formatting (black)
- Linting (flake8)
- Type checking (mypy)
- Import sorting (isort)

## 📊 Improvements Summary

| Category | Improvements |
|----------|-------------|
| **Structure** | Organized modules, proper package hierarchy |
| **Naming** | PEP 8 compliant (snake_case functions) |
| **Documentation** | Type hints, docstrings, examples |
| **Configuration** | YAML/JSON config files |
| **Code Quality** | Removed deprecated code, better error handling |
| **Testing** | Unit tests, fixtures, coverage |
| **Dependencies** | Clear requirements, version pinning |
| **Installation** | One-command pip install |

## 🎨 What's Different from Version 1?

### Original (src/)
```python
from wormer import Wormer
wormer.importGdalRaster('file.tif')
wormer.buildPaddedRaster((2048, 2048))
wormer.wormLevel(100)
```

### Version 2 (version2/)
```python
from bsdwormer import Wormer
wormer.import_raster('file.tif')
wormer.build_padded_raster((2048, 2048))
wormer.worm_level(100)
```

**Same algorithm, better interface!**

## 📦 What's Included

```
version2/
├── bsdwormer/              # Main package
│   ├── __init__.py
│   ├── config.py
│   ├── core/               # Core algorithms
│   ├── utils/              # Utilities
│   ├── io/                 # I/O operations
│   └── visualization/      # VTK output
├── tests/                  # Test suite
│   └── test_*.py
├── examples/               # Usage examples
│   ├── basic_usage.py
│   └── notebooks/
├── docs/                   # Documentation
├── README.md               # Main documentation
├── IMPROVEMENTS.md         # Improvements summary
├── QUICKREF.md            # Quick reference
├── CONTRIBUTING.md         # Contribution guide
├── CHANGELOG.md            # Version history
├── setup.py                # Installation
├── requirements.txt        # Dependencies
├── pytest.ini              # Test config
└── config.example.yaml     # Config example
```

## 🤝 Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## 📄 License

BSD 2-Clause License - See COPYRIGHT.txt in the root directory.

## 👤 Author

Franklin G. Horowitz <frank@horow.net>

## 🔗 Links

- Original Code: `../src/`
- GitHub: https://github.com/RichardScottOZ/BSDWormer
- Issues: https://github.com/RichardScottOZ/BSDWormer/issues

## ❓ Need Help?

1. Check the [README.md](README.md) for detailed documentation
2. See [QUICKREF.md](QUICKREF.md) for common patterns
3. Look at examples in `examples/`
4. Open an issue on GitHub

---

**Note:** The original BSDWormer code (version 1) remains unchanged in the `src/` directory. This version 2 is a modernized implementation with the same core algorithm but improved code organization and usability.
