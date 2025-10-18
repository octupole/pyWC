# Installation

This guide covers installing pyWC on different platforms and configurations.

## Requirements

### System Requirements
- **Python**: ≥ 3.10
- **C++ Compiler**: Supporting C++17 (GCC 7+, Clang 5+, MSVC 2017+)
- **OpenMP**: For CPU parallelization

### Platform Support
- Linux (tested on Ubuntu 20.04+, Fedora, Arch)
- macOS (tested on 11.0+, both Intel and Apple Silicon)
- Windows (tested on Windows 10/11 with Visual Studio 2017+)

## Basic Installation (CPU Only)

### From Source

Clone the repository and install:

```bash
git clone https://github.com/octupole/pyWC.git
cd pyWC
pip install .
```

This will:
- Build C++ extensions with OpenMP support
- Install all required dependencies
- Make command-line tools available

### Dependencies Installed

The basic installation includes:

- `numpy` ≥ 2.1.3 - Numerical arrays
- `scipy` ≥ 1.11.3 - Scientific computing
- `scikit-image` ≥ 0.24.0 - Marching cubes algorithm
- `MDAnalysis` ≥ 2.8.0 - Trajectory I/O
- `packaging` ≥ 23.0 - Version handling

### Build Dependencies

The following are needed at build time:

- `pybind11` ≥ 2.11 - C++/Python bindings
- `Cython` - For DBSCAN extension
- `setuptools` ≥ 61.0
- `wheel`

These are automatically installed during the build process.

## GPU Acceleration (Optional)

For GPU support on NVIDIA GPUs:

### Step 1: Install CUDA Toolkit

Download and install from [NVIDIA CUDA Downloads](https://developer.nvidia.com/cuda-downloads).

Verify installation:
```bash
nvcc --version
nvidia-smi
```

### Step 2: Install CuPy

Choose the appropriate CuPy package for your CUDA version:

=== "CUDA 12.x"
    ```bash
    pip install cupy-cuda12x
    ```

=== "CUDA 11.x"
    ```bash
    pip install cupy-cuda11x
    ```

=== "Auto-detect"
    ```bash
    pip install cupy
    ```

### Step 3: Install pyWC with GPU Support

```bash
pip install .[gpu]
```

Or if installing from source with CUDA detected:
```bash
pip install .  # Automatically includes CuPy if CUDA is detected
```

### Skip GPU Installation

If you have CUDA but don't want GPU support:

```bash
export PYWC_SKIP_CUPY=1
pip install .
```

## Development Installation

For contributing or development:

```bash
git clone https://github.com/octupole/pyWC.git
cd pyWC

# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install in editable mode with dev dependencies
pip install -e .[dev,gpu]
```

This installs additional tools:
- `nose` - Testing framework
- `coverage` - Code coverage analysis

## Platform-Specific Instructions

### Linux

#### Ubuntu/Debian
```bash
# Install build dependencies
sudo apt-get update
sudo apt-get install build-essential python3-dev libomp-dev

# Install pyWC
pip install .
```

#### Fedora/RHEL
```bash
# Install build dependencies
sudo dnf install gcc-c++ python3-devel libomp-devel

# Install pyWC
pip install .
```

#### Arch Linux
```bash
# Install build dependencies
sudo pacman -S base-devel python openmp

# Install pyWC
pip install .
```

### macOS

#### Install Homebrew (if not installed)
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

#### Install OpenMP
```bash
brew install libomp
```

#### Install pyWC
```bash
pip install .
```

!!! note "Apple Silicon (M1/M2/M3)"
    pyWC works on Apple Silicon Macs. Use Python from Homebrew or miniforge for best compatibility:
    ```bash
    brew install python@3.11
    # or
    brew install miniforge
    ```

### Windows

#### Install Visual Studio

1. Download [Visual Studio 2022 Community](https://visualstudio.microsoft.com/downloads/)
2. During installation, select:
   - "Desktop development with C++"
   - Ensure "MSVC v143" and "Windows 10 SDK" are checked

#### Install Python

Download from [python.org](https://www.python.org/downloads/) or use Anaconda.

#### Install pyWC

```cmd
pip install .
```

!!! tip "Windows Subsystem for Linux (WSL)"
    For best experience on Windows, consider using WSL2 with Ubuntu:
    ```bash
    wsl --install -d Ubuntu-22.04
    ```
    Then follow Linux installation instructions.

## Verifying Installation

### Test Basic Import

```python
python -c "from pywc import WillardChandler; print('pyWC installed successfully!')"
```

### Test C++ Extensions

```python
from pywc import WillardChandler
from pywc._wc_kde import evaluate_pbc_fast_auto
print("C++ extensions loaded successfully!")
```

### Test GPU Support (if installed)

```python
from pywc.center_gpu import center_gpu
print("GPU support available!")
```

### Run Command-Line Tools

```bash
pywc-compare-wc-backends --help
pywc-wc-area --help
pywc-bending-rigidity --help
```

## Troubleshooting

### OpenMP Not Found

=== "Linux"
    ```bash
    sudo apt-get install libomp-dev  # Ubuntu/Debian
    sudo dnf install libomp-devel    # Fedora/RHEL
    ```

=== "macOS"
    ```bash
    brew install libomp
    ```

=== "Windows"
    Visual Studio includes OpenMP support. Ensure you selected "Desktop development with C++".

### CUDA Not Detected

```bash
# Set CUDA_HOME environment variable
export CUDA_HOME=/usr/local/cuda         # Linux
export CUDA_HOME=/opt/cuda              # Some Linux systems
set CUDA_HOME=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.0  # Windows
```

### Compilation Errors

If you encounter compilation errors:

1. **Update build tools:**
   ```bash
   pip install --upgrade pip setuptools wheel
   ```

2. **Install specific pybind11 version:**
   ```bash
   pip install pybind11==2.11.1
   ```

3. **Clean build:**
   ```bash
   rm -rf build/ dist/ *.egg-info
   pip install .
   ```

### Import Errors

If you get import errors after installation:

```bash
# Ensure you're not in the source directory
cd ~
python -c "import pywc; print(pywc.__version__)"
```

### Performance Issues

If C++ extensions are not being used:

```python
import pywc
print(pywc._center_impl.HAS_CENTER_FAST)  # Should be True
print(pywc._center_impl.HAS_CENTER_GPU)   # True if CuPy installed
```

## Next Steps

- [Quick Start Guide](quickstart.md) - Learn basic usage
- [API Reference](api/willard_chandler.md) - Detailed API documentation
- [Examples](examples.md) - Real-world usage examples
