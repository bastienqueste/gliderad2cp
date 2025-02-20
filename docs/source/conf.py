# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../.."))


# -- Project information -----------------------------------------------------
import datetime

year = datetime.date.today().year
project = "gliderad2cp"
author = " Bastien Queste, Callum Rollo and Estel Font"
copyright = f"{year}, {author}"  # noqa

# The full version, including alpha/beta/rc tags
# import gliderad2cp
# version = release = gliderad2cp.__version__
version = "v0.0.7"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "nbsphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "IPython.sphinxext.ipython_console_highlighting",
]


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["**.ipynb_checkpoints"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

source_suffix = [".rst", ".md"]

html_logo = "_static/logo.png"

html_context = {
    'display_github': True,
    'github_user': 'bastienqueste',
    'github_repo': 'gliderad2cp',
    'github_version': 'main/docs/source/'
}