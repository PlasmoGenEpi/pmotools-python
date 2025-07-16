# pmotools


A collection of tools to interact with [portable microhaplotype object (pmo) file format](https://github.com/PlasmoGenEpi/portable-microhaplotype-object)

# Setup 

Download repo 
```bash
git clone git@github.com:PlasmoGenEpi/pmotools-python.git
```

Recommend using the conda environment file to ensure all python modules tested with tool are installed 
```bash
cd pmotools-python
conda env create -f envs/pmotools-env.yml 

conda activate pmotools
```

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

# Development  
Development is done the development branch while the master branch will point to the most recent released version following [git-flow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) standards

```bash
# navigate to git directory
cd pmotools-python

# install the environment 
conda env create -f envs/pmotools-env.yml 

# install the current state of the tool in edit mode to continue to work and then update afterwards
conda activate pmotools
pip install -e ./ 

# run main runner to list all currently loaded tools
./scripts/pmotools-runner.py

## or add to path (must be done when within the git directory, or replace $(pwd) with the full path nmae to the git repo
export PATH="$(pwd)/scripts:$PATH"

```
