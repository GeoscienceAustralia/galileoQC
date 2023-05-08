# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.

import pathlib
# import os

import sys
# sys.path.insert(0, os.path.abspath(os.path.join('..', '..', '..', 'AirGravQC')))
# sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'AirGravQC')))
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())
# sys.path.insert(0, pathlib.Path("/Users/markdransfield/Documents/GitHub/AirGravQC/AirGravQC/qualityAnalysys.py").as_posix())

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AirGravQC'
copyright = '2023, GA'
author = 'Mark Dransfield'
release = '0.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'myst_parser',
    'sphinx.ext.napoleon'
]
templates_path = ['_templates']
exclude_patterns = []

myst_enable_extensions = [
    'colon_fence',
    'dollarmath',
    'attrs_block'
]

html_theme_options = {
    'show_powered_by': 'false',
    'body_text': 'Black',
    'font_family': 'Helvetica',
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

numfig = True
