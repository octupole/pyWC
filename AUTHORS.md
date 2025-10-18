# Authors and Contributors

## pyWC (This Fork)

### Maintainer
- **Massimo Marchi** (massimo@octupole.org)
  - Fork creator and maintainer
  - Performance optimizations (C++, GPU backends)
  - Streamlined Willard-Chandler-focused implementation

## Original pytim Authors

pyWC is a fork of [pytim](https://github.com/Marcello-Sega/pytim). The original pytim authors are:

### Original Authors
- **Marcello Sega** (m.sega@ucl.ac.uk) - Original author and maintainer
  - University College London
  - Core architecture and implementation

- **Balazs Fabian**
  - Co-author
  - Algorithm development

- **Gyorgy Hantal**
  - Co-author
  - Implementation and testing

- **Pál Jedlovszky**
  - Co-author
  - Scientific methodology

## Acknowledgments

### Scientific Foundation
The Willard-Chandler method implemented in this package is based on:
- **Adam P. Willard** and **David Chandler** - Original method developers
  - Willard, A. P., & Chandler, D. (2010). "Instantaneous liquid interfaces."
    *The Journal of Physical Chemistry B*, 114(5), 1954-1958.

### Software Dependencies
This project builds upon:
- **MDAnalysis Team** - Molecular dynamics analysis framework
- **NumPy/SciPy Developers** - Scientific computing foundations
- **scikit-image Team** - Marching cubes implementation
- **pybind11 Contributors** - C++/Python binding framework
- **CuPy Team** - GPU acceleration framework

## Fork History

pyWC was created as a specialized fork of pytim to:
1. Focus exclusively on the Willard-Chandler method
2. Implement significant performance improvements for large molecular systems
3. Add GPU acceleration support via CuPy
4. Streamline the codebase by removing unused ITIM/GITIM modules

## Contributing

We welcome contributions! If you contribute to pyWC, you will be added to this file.

To contribute:
1. Fork the repository: https://github.com/octupole/pyWC
2. Create a feature branch
3. Submit a pull request

## Copyright Notice

pyWC:
- Copyright (C) 2025 Massimo Marchi

Original pytim components:
- Copyright (C) 2017-2018 Marcello Sega, Balazs Fabian, Gyorgy Hantal, Pál Jedlovszky

Both are distributed under the GNU General Public License v3.0 (GPL-3.0-only).
See the [LICENSE](LICENSE) file for full details.
