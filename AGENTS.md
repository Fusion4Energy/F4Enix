# F4Enix — AI Agent Instructions

F4Enix is a Python API for parsing and manipulating Monte Carlo simulation input/output files (primarily MCNP and related formats). Developed by the F4E Radiation Transport team.

## Essential Commands

```bash
# Install in editable mode with test dependencies
pip install -e .[tests]

# Run full test suite with coverage
pytest --cov=src/f4enix

# Run a specific test module
pytest tests/materials_test.py

# Run sub-package tests (e.g. weight windows)
pytest tests/test_ww_gvr/
```

> Linux only: PyVista tests require a display. CI uses a virtual framebuffer (`Xvfb`). Set `DISPLAY` or use `pyvista.start_xvfb()` if running headlessly.

## Project Layout

```
src/f4enix/
  input/      # Input file parsers (MCNP, materials, ACE, d1suned, ww_gvr/, …)
  output/     # Output file parsers (MCNP .o, mctal, meshtal/, rssa/, cdgs/, …)
  core/       # Shared: constants, spectra, irradiation, material_library, egroups
tests/
  resources/  # Test fixtures organized as importable packages (one sub-dir per module)
  test_ww_gvr/ # Sub-package tests for ww_gvr/
```

The package `__init__.py` re-exports the four most-used classes: `Input`, `Output`, `Mctal`, `Meshtal`.

## Architecture

- **`input/`**: Parsers model the hierarchical structure of file formats. E.g. `MatCardsList → Material → SubMaterial → Element → Zaid` mirrors MCNP material card structure.
- **`output/`**: Parsers accept a filepath in `__init__` and parse immediately (`Output(path)`).
- **`core/constants.py`**: Centralizes all shared regex patterns (`PAT_COMMENT`, `PAT_MAT`, `SCIENTIFIC_PAT`, etc.). Import from here rather than defining patterns locally.
- **Sub-packages** (`ww_gvr/`, `meshtal/`, `rssa/`, `cdgs/`) have dedicated test directories alongside the flat `tests/` files.

## Key Conventions

### Class Construction
There are two patterns — check the class before instantiating:
- **Direct init**: `Output(filepath)`, `Mctal(filepath)` — parses in `__init__`
- **Classmethods**: `MatCardsList.from_input(path)`, `Zaid.from_string(text)` — preferred for complex parsers

### Test Conventions
- Tests are **class-based**: `class TestMaterial:`, `class TestZaid:` grouping related assertions.
- Resources use `importlib.resources.files(...)` — **never** raw `__file__`-relative paths.
- Expensive fixtures (e.g. `LibManager`) are constructed **once at module level** and reused across all test classes.
- Two test filename conventions coexist (`materials_test.py` vs `test_cuv_sampling_error.py`). Prefer `test_<module>.py` for new files.

### License Header
Every new source file must include the standard **EUPL-1.2** copyright header (copy from any existing module).

### Dependencies
- `numjuggler` — MCNP input tokenizer; `chardet < 7.0.0` is pinned because newer versions break it — do not bump `chardet`.
- `pyvista` — mesh visualization; only needed for `plotter.py`-related code.
- `polars` and `pandas` both used for tabular data — check which a module uses before adding dataframe operations.

### Versioning
Version is managed by `setuptools_scm` and written to `src/_version.py`. Always install with `pip install -e .` to avoid version warnings.

## Branching Model
- `developing` — integration branch; feature branches target here
- PRs to `developing` require an independent reviewer approval before merge

## Useful Links
- [Full documentation](https://f4enix.readthedocs.io/en/latest/)
- [Installation guide](docs/source/usage/installation.rst)
- [Contributing guide](docs/source/developers/contributing.rst)
- [CI workflows](.github/workflows/)
