# pmotools

A collection of tools to interact with the [portable microhaplotype object (PMO) file format](https://github.com/PlasmoGenEpi/portable-microhaplotype-object).

Documentation on the format and pmotools-python, including detailed tutorials can be found [here](https://plasmogenepi.github.io/PMO_Docs/).

The manual for pmotools-python can be found [here](https://plasmogenepi.github.io/pmotools-python/index.html).

# Installation

Requires Python 3.11+.

```bash
pip install pmotools
```

This installs the `pmotools-python` command-line tool and the `pmotools` Python library.

# Usage

pmotools-python can be used from the command line or imported as a Python library.

## Command line

After installation, the `pmotools-python` command is available:

```bash
pmotools-python --help
pmotools-python validate_pmo --help
```

Examples:

```bash
pmotools-python validate_pmo --pmo my_data.pmo.json
pmotools-python combine_pmos --pmo_files run1.pmo.json run2.pmo.json --output combined.pmo.json
```

Commands are grouped by task:

| Group | Purpose |
| --- | --- |
| `convertors_to_json` | Convert tables and metadata into PMO inputs |
| `extractors_from_pmo` | Subset PMOs by metadata, samples, targets, or read counts |
| `pmo_to_table` | Export PMO contents to spreadsheets or other formats |
| `extract_basic_info_from_pmo` | Inspect metadata and counts |
| `working_with_multiple_pmos` | Combine PMOs |
| `validation` | Validate PMO files against the JSON schema |

See the [command-line manual](https://plasmogenepi.github.io/pmotools-python/commands/index.html) for full documentation of each command.

## Python library

Import pmotools in your own code to read, validate, write, and build PMO files:

```python
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.pmo_engine.pmo_checker import PMOChecker

pmo = PMOReader.read_in_pmo("example.pmo.json")

checker = PMOChecker()
checker.validate_pmo_json(pmo)
```

Core modules:

- `pmotools.pmo_engine` — read, write, validate, and export PMOs
- `pmotools.pmo_builder` — build PMOs from tables and metadata
- `pmotools.scripts` — the same logic used by the CLI

See the [API reference](https://plasmogenepi.github.io/pmotools-python/api/modules.html) for module-level documentation.

## Auto completion

To enable bash tab completion for `pmotools-python`, append the script from `etc/bash_completion` to your `~/.bash_completion`, or generate it with:

```bash
pmotools-python --bash-completion >> ~/.bash_completion
source ~/.bash_completion
```

## Developer Setup

To contribute to `pmotools`, follow these steps:

1. **Clone the repository** and switch to the develop branch:
```bash
git clone git@github.com:your-org/pmotools.git
cd pmotools
git checkout develop
```

2. **Create your feature branch**:
```bash
git checkout -b feature/my-feature
```

3. **Install and set up UV.** This creates .venv/ and installs everything from pyproject.toml:
```bash
pip install -U uv
uv sync --dev
```

4. **Install pre-commit hooks** (for formatting & linting):
```bash
uv run pre-commit install
```

5. **Run pre-commit** manually on all files (first time):
```bash
uv run pre-commit run --all-files
```

6. **Develop your code**. Pre-commit will automatically run on staged files before each commit, checking:
* Formatting (Ruff)
* Linting (Ruff)
* Trailing whitespace, YAML syntax, large files

7. **Run tests**:
```bash
uv run pytest
```

8. **Commit and push** your changes:
```bash
git add .
git commit -m "Your message"
git push origin feature/my-feature
```

### Documentation updating

Documentation for pmotools is automatically generated from docstring under `man/`. This is automatically built and deployed through GitHub actions on merging to main.

You should check documentation before deployment. To update the documentation locally make sure you have pmotools-python installed.

Next, from your development environment, install sphinx and its dependencies using the following command:

```bash
pip install -r man/requirements.txt
```

Build the documentation using the following commands
```bash
cd man
make update_autodocs
make html
```
You can open the html to review changes.

### Releasing pmotools-python

To release pmotools-python you need to update the version in two places:
* update `version` in `pyproject.toml`, adhering to semantic versioning conventions.
* update `release` in the documentation by updating `man/source/conf.py`

**Note:** It is not always the case, but sometimes a release of pmotools-python will coincide with a new schema version. See the section below for notes on updating the schema version.

Once you have merged these updates into the `main` branch, create a new release with notes. GitHub actions will automatically deploy to PyPi.

### Updating the schema version

You can release pmotools-python without updating the schema version. However, for the schema to be updated and become the default within pmotools-python, you must update the schema and then do a release of pmotools-python by following the instructions above.

Update `__schema_version__` in `src/pmotools/__init__.py` to the new version of the schema. Do not include the `v` here (e.g. 1.0.0 not v1.0.0).
