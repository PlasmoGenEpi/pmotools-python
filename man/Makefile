# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Paths relative to this Makefile (which lives in man/)
PACKAGEDIR    = ../src/pmotools
CLIGEN        = gen_cli_docs.py

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Regenerate all auto-generated .rst sources: the API stubs (sphinx-apidoc)
# and the per-command CLI pages (gen_cli_docs.py). Run this after adding or
# renaming modules or CLI commands, then build with `make html`.
APIDIR        = $(SOURCEDIR)/api
CLIDIR        = $(SOURCEDIR)/commands

update_autodocs:
	rm -rf "$(APIDIR)" "$(CLIDIR)"
	sphinx-apidoc -f -e -o "$(APIDIR)/" "$(PACKAGEDIR)/"
	python "$(CLIGEN)"

.PHONY: help Makefile update_autodocs

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
