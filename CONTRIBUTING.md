# Contributing to pyWC

Thank you for considering contributing to pyWC! This document provides guidelines for contributing to the project.

## About pyWC

pyWC is a specialized fork of [pytim](https://github.com/Marcello-Sega/pytim) that focuses exclusively on the Willard-Chandler method with performance optimizations for large molecular systems.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue on [GitHub Issues](https://github.com/octupole/pyWC/issues) with:

- A clear, descriptive title
- Steps to reproduce the issue
- Expected behavior vs. actual behavior
- Your environment (OS, Python version, GPU/CPU info)
- Minimal code example that reproduces the problem

### Suggesting Enhancements

We welcome suggestions for improvements! Please open an issue with:

- A clear description of the enhancement
- Use cases and benefits
- Any implementation ideas you have

### Pull Requests

1. **Fork the repository**
   ```bash
   git clone https://github.com/octupole/pyWC.git
   cd pyWC
   ```

2. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

3. **Set up development environment**
   ```bash
   # Create virtual environment
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate

   # Install in development mode with dev dependencies
   pip install -e .[dev,gpu]
   ```

4. **Make your changes**
   - Write clear, documented code
   - Follow existing code style
   - Add tests if applicable
   - Update documentation as needed

5. **Test your changes**
   ```bash
   # Run basic tests
   python -m pytest pywc/

   # Test C++ extensions
   python -c "from pywc import WillardChandler; print('Import successful')"

   # Test GPU backend (if CuPy installed)
   python scripts/compare_cpu_gpu.py
   ```

6. **Commit your changes**
   ```bash
   git add .
   git commit -m "Brief description of changes"
   ```

7. **Push and create pull request**
   ```bash
   git push origin feature/your-feature-name
   ```

   Then open a pull request on GitHub with:
   - Clear description of changes
   - Reference to related issues (if applicable)
   - Any breaking changes noted

## Development Guidelines

### Code Style

- **Python**: Follow PEP 8 style guide
  - Use 4 spaces for indentation
  - Maximum line length: 100 characters (flexible for readability)
  - Use descriptive variable names

- **C++**: Follow existing style in `.cpp` files
  - Use consistent indentation (4 spaces)
  - Comment complex algorithms
  - Use modern C++17 features

- **Documentation**:
  - Use NumPy-style docstrings
  - Document all public functions and classes
  - Include parameter types and return values

### Example Docstring

```python
def compute_surface(positions, alpha, mesh):
    """
    Compute Willard-Chandler surface from atomic positions.

    Parameters
    ----------
    positions : np.ndarray, shape (N, 3)
        Atomic coordinates in Ångströms
    alpha : float
        Gaussian width parameter in Ångströms
    mesh : float
        Grid spacing in Ångströms

    Returns
    -------
    vertices : np.ndarray, shape (M, 3)
        Surface vertex coordinates
    faces : np.ndarray, shape (K, 3)
        Triangle indices
    normals : np.ndarray, shape (M, 3)
        Vertex normal vectors

    Raises
    ------
    ValueError
        If alpha or mesh are non-positive
    """
```

### Testing

- Add tests for new features in appropriate test files
- Ensure existing tests pass
- Test both CPU and GPU backends if applicable
- Include edge cases and error conditions

### Performance

Since pyWC focuses on performance:

- **Profile your code** before optimizing
- **Document performance improvements** in pull requests
- **Benchmark** against existing implementation
- Consider both **CPU and GPU** paths

### Commit Messages

Write clear commit messages:

```
Short (50 chars or less) summary

More detailed explanatory text, if necessary. Wrap at 72 characters.
Explain the problem that this commit is solving, not just what changed.

- Use bullet points for multiple changes
- Reference issues: Fixes #123
```

## Areas for Contribution

We especially welcome contributions in these areas:

### Performance Improvements
- Further GPU optimizations
- CPU SIMD vectorization
- Memory usage reduction
- Faster grid generation

### Features
- Additional surface analysis metrics
- Better visualization tools
- Trajectory analysis utilities
- Support for new output formats

### Documentation
- Usage examples and tutorials
- Performance benchmarking guides
- API reference improvements
- Visualization examples with ParaView/VMD

### Testing
- Expanded test coverage
- Benchmark suite
- Integration tests
- Cross-platform testing

### Build System
- Improved CUDA detection
- Better error messages for build failures
- Support for more compilers
- Conda package recipe

## Building C++ Extensions

The C++ extensions require:
- C++17 compatible compiler
- OpenMP support
- pybind11

### Linux/macOS

```bash
# Install build dependencies
pip install pybind11 numpy scipy

# Build extensions
python setup.py build_ext --inplace

# Or install in development mode
pip install -e .
```

### Windows

```bash
# Requires Visual Studio 2017 or newer
pip install pybind11 numpy scipy
python setup.py build_ext --inplace
```

### Troubleshooting Build Issues

**OpenMP not found:**
- Linux: Install `libomp-dev`
- macOS: Install via Homebrew: `brew install libomp`
- Windows: Use Visual Studio with OpenMP support

**CUDA not detected:**
```bash
# Set CUDA_HOME environment variable
export CUDA_HOME=/usr/local/cuda  # Linux/macOS
set CUDA_HOME=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.0  # Windows
```

**Skip CuPy installation:**
```bash
export PYWC_SKIP_CUPY=1
pip install .
```

## Code of Conduct

### Our Standards

- Be respectful and inclusive
- Welcome newcomers
- Accept constructive criticism
- Focus on what's best for the community
- Show empathy toward others

### Our Responsibilities

Maintainers are responsible for clarifying standards and taking appropriate action in response to unacceptable behavior.

## Questions?

- Open an issue for general questions
- Email the maintainer: massimo@octupole.org
- Check the [README](README.md) for basic usage

## Attribution

Contributors will be acknowledged in [AUTHORS.md](AUTHORS.md).

## License

By contributing, you agree that your contributions will be licensed under the GNU General Public License v3.0 (GPL-3.0), the same license as pyWC and its upstream project pytim.

---

Thank you for contributing to pyWC!
