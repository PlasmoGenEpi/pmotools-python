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
conda env activate -f envs/pmotools-env.yml 
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
