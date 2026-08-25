# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path
import numpydoc
import myst_parser

srcdir = Path('../../src').resolve()
sys.path.insert(0, str(srcdir))

# -- Project information -----------------------------------------------------

project = "ribbon"
copyright = "2026, Sarit Dutta"
author = "Sarit Dutta"
version = "0.5.0"
release = "0.5.0"

# -- General configuration ---------------------------------------------------
extensions = [
        "sphinx.ext.duration",
        "sphinx.ext.autodoc",
        "sphinx.ext.autosummary",
        "sphinx.ext.viewcode",
        "sphinx.ext.intersphinx",
        "myst_parser",
        "numpydoc",
        ]

exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "inherited-members": True,
    "show-inheritance": True,
    "member-order": "bysource",
}
autosummary_generate = True
toc_object_entries_show_parents = "hide"

#Global TeX macros
macros = {
        'bs': [r'\boldsymbol{#1}', 1],
        'bm': [r'\mathbf{#1}', 1],
        'omag': [r'\lVert \boldsymbol{\upomega} \rVert',0],
         }
mathjax4_config = {
    'loader': {'load': ['[tex]/upgreek']},
    'tex': {
        'packages': {'[+]': ['upgreek']},
        'macros': macros
        }
}

intersphinx_mapping = {
        "numpy": ("https://numpy.org/doc/stable/", None),
        "scipy [latest]": ("https://docs.scipy.org/doc/scipy/", None),
        }

# MyST parser extensions
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

# -- Options for HTML output -------------------------------------------------
html_theme = "pydata_sphinx_theme"

html_static_path = ["_static"]
html_css_files = ["custom.css"]

# SmartyPants will be used to convert quotes and dashes to typographically
# correct entities.
html_use_smartypants = True
html_theme_options = {
    "primary_sidebar_end": ["indices.html"]
}
#html_theme_options = {
#    "use_edit_page_button": True,
#}
