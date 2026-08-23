# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path
import numpydoc

srcdir = Path('../../src').resolve()
sys.path.insert(0, str(srcdir))

# -- Project information -----------------------------------------------------

project = 'Spirauliya'
copyright = '2026, Sarit Dutta'
author = 'Sarit Dutta'
release = '0.5.0'

# -- General configuration ---------------------------------------------------
extensions = [
        'sphinx.ext.duration',
        'sphinx.ext.autodoc',
        'sphinx.ext.autosummary',
        'sphinx.ext.viewcode',
        'sphinx.ext.intersphinx',
        'numpydoc',
        ]

exclude_patterns = []

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "inherited-members": True,
    "show-inheritance": True,
    "member-order": 'bysource',
}
autosummary_generate = True
toc_object_entries_show_parents = "hide"

intersphinx_mapping = {'numpy': ('https://numpy.org/doc/stable/', None)}

# -- Options for HTML output -------------------------------------------------
html_theme = 'pydata_sphinx_theme'

html_static_path = ['_static']
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
