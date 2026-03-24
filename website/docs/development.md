# Development

## Documentation source

The documentation site is maintained under `website/` and published with GitHub Pages through `.github/workflows/docs.yml`.

## Local workflow

```bash
pip install -e ".[docs]"
pip install -r website/requirements.txt
mkdocs serve -f website/mkdocs.yml
```

## Content policy

- Prefer current package docstrings over duplicated narrative text.
- Keep example references pointed at notebooks in `tutorials/`.
- Avoid reintroducing generated Sphinx files under `website/source/`.
