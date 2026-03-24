# Website and Documentation (GitHub Pages)

This folder contains the MkDocs source for the public documentation site published with GitHub Pages.

## Local preview

```bash
pip install -e ".[docs]"
pip install -r website/requirements.txt
mkdocs serve -f website/mkdocs.yml
```

Then open `http://127.0.0.1:8000/`.

## Structure

```text
website/
├── docs/                 # Markdown source pages
├── mkdocs.yml            # Site navigation and theme configuration
└── requirements.txt      # Extra packages for local docs builds
```

## Publishing

GitHub Actions builds this folder and deploys it to GitHub Pages using `.github/workflows/docs.yml`.

## Editing guidance

- Keep tutorials aligned with the notebooks in `tutorials/`.
- Prefer current package docstrings as the source of truth for API behavior.
- When adding a page, register it in `mkdocs.yml`.
