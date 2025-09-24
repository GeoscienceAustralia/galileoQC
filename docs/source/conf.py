# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.

import pathlib
import importlib.metadata

import sys
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pegasusQC'
copyright = '2023, Mark Helm Dransfield'
author = 'Mark Helm Dransfield'
# version = '0.0.1'
# release = '0.0.1'
version = release = importlib.metadata.version('pegasusQC')
today_fmt = '%d %b %Y' 

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'myst_nb',
]

templates_path = ['_templates']

exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "Thumbs.db",
    ".DS_Store",
    ".env",
    ".venv",
]

myst_enable_extensions = [
    'colon_fence',
    'dollarmath',
    'attrs_block'
]

myst_footnote_transition = False

source_suffix = {
    '.rst': 'restructuredtext',
    '.ipynb': 'myst-nb',
}

nbsphinx_execute = "auto"

nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'png2x'}",
    "--InlineBackend.rc=figure.dpi=96",
]

nbsphinx_kernel_name = "python3"


# html_sidebars = { '**': ['globaltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html'] }

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_static_path = ['static']
# html_logo = 'static/pegasus_logo.png'
# html_favicon = "static/pegasus_icon.png"
html_title = "pegasusQC"
html_theme_options = {
    "light_logo": "pegasus_logo_light.png",
    "dark_logo": "pegasus_logo_dark.png",
    # "sidebar_hide_name": True,
    # 'show_powered_by': 'false',
    # 'body_text': 'Black',
    # 'font_family': 'Palatino',
}

numfig = True
