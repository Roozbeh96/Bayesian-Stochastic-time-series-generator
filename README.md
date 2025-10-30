# Bayesian-Stochastic generation of (spatio/temporal)-series signal

Tool to stochastically generate (spatio/temporal)-series velocity signal. This velocity signal can be used in stochastic generation of turbulent boundary layer velocity field.

## Contents
- [python package](#python-package)

## Python package
<!-- brief blurb or link to docs/install/usage -->
<details>
  <summary>Click to see the documentation for Python files</summary>

### Python files
<details>
<summary> Installation:</summary>

1. 🛠️ Installing Poetry

To install [Poetry](https://python-poetry.org/) (Python dependency management and packaging tool), run the following command in your terminal:

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

After installation, make sure Poetry is in your `PATH`
- macOS/Linux
```bash
export PATH="$HOME/.local/bin:$PATH
```

Verify installation:

```bash
poetry version
```

Keep venv inside the project (works great with VS Code) poetry 

```bash
config virtualenvs.in-project true
```

2. Installing environment: 

This environment is set with Python 3.13. Change the requires-python = ">=3.13" in pyproject.tmol file if you have other versions on your PC.

run:
```bash
poetry install
```

In case if you want to make environment from scratch, run:(Do not recommended)

```bash
poetry new project_name
```
</details>

</details>