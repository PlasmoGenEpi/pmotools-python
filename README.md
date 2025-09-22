# pmotools

A collection of tools to interact with [portable microhaplotype object (pmo) file format](https://github.com/PlasmoGenEpi/portable-microhaplotype-object)

# Setup

Install using pip
```bash
pip install .
```

# Usage

This package is built to either be used as a library in python projects and then are also command line already created which can be found within the scripts directory. All scripts are also all wrapped into the main runner file `scripts/pmotools-runner.py`. The scripts do require the pmotools-python package to be installed like with pip above.


## Auto completion

If you want to add auto-completion to the scripts master function [scripts/pmotools-runner.py](scripts/pmotools-runner.py) you can add the following to your `~/.bash_completion`. This can also be found in etc/.bash_completion in current directory.

```bash
_comp_pmotools_runner()
{
    local cur prev opts base
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    if [[ $COMP_CWORD -lt 2 ]] ; then
    	opts=$(for x in `${COMP_WORDS[0]} | grep "^\s.*-" | sed 's/ -.*//g' | tr -d '[:blank:]'`; do echo ${x} ; done )
		COMPREPLY=($(compgen -W "${opts}" -- ${cur}))
    elif [[ ${cur} == -* ]]; then
    	opts=$(for x in `${COMP_WORDS[0]} ${COMP_WORDS[1]} -h | grep " -" | sed "s/^. *-/-/g" | sed "s/   .*//g" | sed "s/, / /g"`; do echo ${x} ; done )
		COMPREPLY=($(compgen -W "${opts}" -- ${cur}))
    else
    	_filedir
    fi
   return 0
}


complete -F _comp_pmotools_runner pmotools-runner.py

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
