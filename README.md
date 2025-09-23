# pmotools

A collection of tools to interact with [portable microhaplotype object (pmo) file format](https://github.com/PlasmoGenEpi/portable-microhaplotype-object)

# Setup

Install using pip
```bash
pip install .
```

# Usage

This package is built to either be used as a library in python projects and a command line interface already created which can be called from the commandline `pmotools-python` which will install with `pip install .`.


## Auto completion

If you want to add auto-completion to the scripts master function [pmotools-python](scripts/pmotools-runner.py) you can add the following to your `~/.bash_completion`. This can also be found in etc/bash_completion in the current directory. Or can be generated with `pmotools-python --bash-completion`

```bash
_pmotools_python_complete()
{
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # 1) Completing the command name (1st arg): list all commands
    if [[ ${COMP_CWORD} -eq 1 ]]; then
        # Our CLI prints machine-friendly list via --list-plain:
        # "<command>\t<group>\t<help>"
        local lines cmds
        lines="$(${COMP_WORDS[0]} --list-plain 2>/dev/null)"
        cmds="$(printf '%s\n' "${lines}" | awk -F'\t' '{print $1}')"
        COMPREPLY=( $(compgen -W "${cmds}" -- "${cur}") )
        return 0
    fi

    # 2) Completing flags for a leaf command: scrape leaf -h
    if [[ "${cur}" == -* ]]; then
        local helps opts
        helps="$(${COMP_WORDS[0]} ${COMP_WORDS[1]} -h 2>/dev/null)"
        # Pull out flag tokens and split comma-separated forms
        opts="$(printf '%s\n' "${helps}" \
            | sed -n 's/^[[:space:]]\{0,\}\(-[-[:alnum:]][-[:alnum:]]*\)\(, *-[[:alnum:]][-[:alnum:]]*\)\{0,\}.*/\1/p' \
            | sed 's/, / /g')"
        COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
        return 0
    fi

    # 3) Otherwise, fall back to filename completion for positional args
    COMPREPLY=( $(compgen -f -- "${cur}") )
    return 0
}

complete -F _pmotools_python_complete pmotools-python
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
