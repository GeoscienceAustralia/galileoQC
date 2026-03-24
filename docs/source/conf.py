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

project = 'galileoQC'
copyright = '2025, Geoscience Australia'
author = 'Mark Helm Dransfield'

try:
    version = release = importlib.metadata.version('galileoQC')
except importlib.metadata.PackageNotFoundError:
    version = release = '0.0.4'  # fallback version


# version = release = importlib.metadata.version('galileoQC')
today_fmt = '%d %b %Y' 

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # 'sphinx.ext.autodoc',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'myst_nb',
    'autodoc2',
]

# Napoleon settings
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# autodoc2 settings
autodoc2_packages = [
    {
        "path": "../../galileoQC",
        "auto_mode": True,
    },
]
autodoc2_skip_module_regexes = [r'^_.*'] # This regex matches any name that starts with an underscore
autodoc2_render_plugin = "myst"

# myst settings
myst_enable_extensions = ["fieldlist"]

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

# jupyter notebook settings
nbsphinx_execute = "auto"
nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'png2x'}",
    "--InlineBackend.rc=figure.dpi=96",
]
nbsphinx_kernel_name = "python3"
nb_execution_timeout = 180

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

if False:
    html_theme = "furo"
    html_theme_options = {
        # "light_logo": "pegasus_logo_light.png",
        # "dark_logo": "pegasus_logo_dark.png",
        "sidebar_hide_name": False,
        "version_selector": True,
    }

if True:
    html_theme = "alabaster"
    html_theme_options = {
        # 'logo': './pegasus_logo_light.png',
        # 'logo_name': True,
        'description': 'Quality control for airborne gravity surveys.',
        'fixed_sidebar': 'true',
        'show_relbars': 'true',
        'show_relbar_top': 'false',
    }

numfig = True
html_title = "galileoQC"
html_static_path = ['static']

