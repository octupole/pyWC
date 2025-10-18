# pyWC Documentation

This directory contains the source files for the pyWC documentation website.

## Building Locally

### Install Dependencies

```bash
pip install -r requirements.txt
```

### Preview Documentation

```bash
mkdocs serve
```

Then visit http://127.0.0.1:8000 in your browser.

### Build Static Site

```bash
mkdocs build
```

The built site will be in `site/` directory.

## Deployment

Documentation is automatically built and deployed to GitHub Pages when changes are pushed to the `main` branch.

The deployed site is available at: https://octupole.github.io/pyWC/

## Structure

```
docs/
├── index.md              # Homepage
├── installation.md       # Installation guide
├── quickstart.md         # Quick start tutorial
├── examples.md           # Usage examples
├── api/
│   ├── willard_chandler.md  # Main API reference
│   └── backends.md          # Backend documentation
├── citation.md           # How to cite
├── contributing.md       # Contributing guide
├── stylesheets/
│   └── extra.css        # Custom CSS
├── javascripts/
│   └── mathjax.js       # MathJax configuration
└── requirements.txt      # Build dependencies
```

## Updating Documentation

1. Edit the relevant `.md` files in `docs/`
2. Preview changes with `mkdocs serve`
3. Commit and push to `main` branch
4. GitHub Actions will automatically deploy

## Theme

We use [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) for a modern, responsive design with:

- Light/dark mode toggle
- Code syntax highlighting
- Search functionality
- Navigation tabs
- Mobile-friendly layout

## Markdown Extensions

Available extensions:

- **Admonitions**: Info boxes with `!!! note`, `!!! warning`, etc.
- **Code blocks**: Syntax highlighting with line numbers
- **Tables**: Standard Markdown tables
- **Math**: LaTeX math with MathJax (`$...$` or `$$...$$`)
- **Tabs**: Content tabs with `=== "Tab 1"`
- **Icons/Emojis**: `:material-icon:` syntax

## Contributing to Docs

See the main [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

For documentation-specific contributions:
- Follow existing page structure
- Use clear, concise language
- Include code examples where applicable
- Test locally before pushing
