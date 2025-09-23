# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "pmotools-python"
copyright = "2024, Plasmogenepi"
author = "Plasmogenepi"
release = "1.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.duration",  # will report duration while building
    "sphinx.ext.doctest",  # can run doc checks
    "sphinx.ext.autodoc",  # allows for autofucntion documentation from documented code
    "sphinx.ext.autosummary",  # allows for doing a summary of code with autodoc from code
    "sphinx.ext.githubpages",  # add .nojekyll to gh-pages
    "sphinx_copybutton",  # add copy button to code chunks
    "sphinx_toolbox.github",  # link to github
    "sphinx_licenseinfo",  # add license information
    "notfound.extension",  # 404 page
    "sphinx.ext.autosectionlabel",  # reference sections using their title
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.linkcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
html_static_path = ["_static"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_context = {
    "display_github": True,
    "github_user": "PlasmoGenEpi",
    "github_repo": "pmotools-python",
    "github_version": "main/docs/",
}


# -- Sphinx Toolbox configuration-----------------------------------------------
github_username = "PlasmoGenEpi"
github_repository = "pmotools-python"


# -- 404 Page configuration-----------------------------------------------------
notfound_urls_prefix = "/pmotools-python/"


# -- Auto Section configuration-------------------------------------------------
# Make sure the target is unique
autosectionlabel_prefix_document = True


# -- linkcode configuration ----------------------------------------------------
def linkcode_resolve(domain, info):
    if domain != "py":
        return None
    if not info["module"]:
        return None
    filename = info["module"].replace(".", "/")
    return (
        "https://github.com/PlasmoGenEpi/pmotools-python/tree/develop/%s.py" % filename
    )


# code blocks
pygments_style = "sphinx"
