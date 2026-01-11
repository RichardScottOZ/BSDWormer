# Contributing to BSDWormer

Thank you for your interest in contributing to BSDWormer! This document provides guidelines for contributing to the project.

## Code of Conduct

Please be respectful and constructive in all interactions. We aim to maintain a welcoming and inclusive community.

## How to Contribute

### Reporting Bugs

If you find a bug, please create an issue on GitHub with:
- A clear, descriptive title
- Steps to reproduce the problem
- Expected behavior vs actual behavior
- Your environment (OS, Python version, GDAL version)
- Sample data or code if possible

### Suggesting Enhancements

Enhancement suggestions are welcome! Please create an issue with:
- A clear description of the enhancement
- Use cases and benefits
- Any potential drawbacks or alternatives considered

### Pull Requests

1. **Fork the repository** and create a branch from `main`
2. **Make your changes** following our coding standards
3. **Add tests** for any new functionality
4. **Update documentation** as needed
5. **Run the test suite** to ensure everything passes
6. **Submit a pull request** with a clear description

## Development Setup

### Prerequisites

- Python 3.7 or higher
- GDAL system library (3.0+)
- Git

### Setting Up Your Environment

```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/BSDWormer.git
cd BSDWormer/version2

# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install development dependencies
pip install -e ".[dev]"

# Install pre-commit hooks (optional but recommended)
pip install pre-commit
pre-commit install
```

## Coding Standards

### Style Guide

We follow PEP 8 with some modifications:
- Line length: 88 characters (Black default)
- Use type hints for function signatures
- Write docstrings for all public functions/classes

### Code Formatting

We use `black` for automatic code formatting:

```bash
black bsdwormer/
```

### Linting

Run `flake8` to check code quality:

```bash
flake8 bsdwormer/
```

### Type Checking

Run `mypy` for static type checking:

```bash
mypy bsdwormer/
```

### Import Sorting

Use `isort` to sort imports:

```bash
isort bsdwormer/
```

## Testing

### Running Tests

Run the full test suite:

```bash
pytest tests/
```

Run with coverage:

```bash
pytest --cov=bsdwormer tests/
```

Run specific test files:

```bash
pytest tests/test_wormer.py
```

### Writing Tests

- Place tests in the `tests/` directory
- Name test files `test_*.py`
- Name test functions `test_*`
- Use descriptive test names
- Include docstrings explaining what each test validates
- Use fixtures for shared test data

Example:

```python
def test_fourier_domain_grid_initialization():
    """Test that FourierDomainGrid initializes with correct defaults."""
    fdg = FourierDomainGrid(dx=1.0, dy=1.0)
    assert fdg.dx == 1.0
    assert fdg.dy == 1.0
    assert fdg.spatial_grid is None
```

## Documentation

### Docstrings

Use Google-style docstrings:

```python
def function_name(param1: str, param2: int) -> bool:
    """Short description.
    
    Longer description if needed.
    
    Args:
        param1: Description of param1
        param2: Description of param2
        
    Returns:
        Description of return value
        
    Raises:
        ValueError: When something goes wrong
        
    Example:
        >>> function_name("hello", 42)
        True
    """
    pass
```

### Building Documentation

```bash
cd docs
make html
```

View documentation at `docs/_build/html/index.html`

## Git Workflow

### Branches

- `main`: Stable, production-ready code
- `develop`: Integration branch for features
- `feature/feature-name`: New features
- `bugfix/bug-name`: Bug fixes
- `docs/topic`: Documentation improvements

### Commit Messages

Write clear, concise commit messages:

```
Short (50 chars or less) summary

More detailed explanatory text, if necessary. Wrap it to about 72
characters. The blank line separating the summary from the body is
critical.

- Bullet points are okay
- Use present tense ("Add feature" not "Added feature")
- Reference issues and pull requests
```

Examples:
- `Add configuration file support`
- `Fix upward continuation operator scaling`
- `Update documentation for FourierDomainGrid`

### Pull Request Process

1. Update the CHANGELOG.md with your changes
2. Update documentation if needed
3. Ensure all tests pass
4. Get at least one code review
5. Squash commits if requested
6. Maintainer will merge when ready

## Questions?

If you have questions about contributing, please:
- Check existing issues and documentation
- Create a new issue with the "question" label
- Contact the maintainers directly

## License

By contributing, you agree that your contributions will be licensed under the BSD 2-Clause License.

## Acknowledgments

Thank you for helping improve BSDWormer! Your contributions make this project better for everyone.
