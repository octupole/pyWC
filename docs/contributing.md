# Contributing to pyWC

Thank you for considering contributing to pyWC! This page provides guidelines for contributing to the project.

For the complete contribution guide, please see [CONTRIBUTING.md](https://github.com/octupole/pyWC/blob/main/CONTRIBUTING.md) in the repository.

## Quick Links

- **Report Bugs**: [GitHub Issues](https://github.com/octupole/pyWC/issues)
- **Suggest Features**: [GitHub Issues](https://github.com/octupole/pyWC/issues)
- **Repository**: [github.com/octupole/pyWC](https://github.com/octupole/pyWC)

## How to Contribute

### 1. Fork and Clone

```bash
git clone https://github.com/octupole/pyWC.git
cd pyWC
```

### 2. Set Up Development Environment

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -e .[dev,gpu]
```

### 3. Make Changes

- Write clear, documented code
- Follow existing code style (PEP 8)
- Add tests if applicable
- Update documentation

### 4. Submit Pull Request

```bash
git checkout -b feature/your-feature
git add .
git commit -m "Description of changes"
git push origin feature/your-feature
```

Then open a pull request on GitHub.

## Areas for Contribution

We especially welcome contributions in:

- **Performance improvements**: GPU optimizations, SIMD vectorization
- **Documentation**: Examples, tutorials, API docs
- **Testing**: Expanded test coverage, benchmarks
- **Features**: New analysis metrics, output formats

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Accept constructive criticism
- Focus on what's best for the community

## Questions?

- Open an issue for questions
- Email: massimo@octupole.org
- See full guide: [CONTRIBUTING.md](https://github.com/octupole/pyWC/blob/main/CONTRIBUTING.md)
