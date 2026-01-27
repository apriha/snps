#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# snps documentation build configuration file
#

import os
import sys

import snps

sys.path.insert(0, os.path.abspath("../src"))

# -- General configuration ------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx.ext.intersphinx",
    "myst_parser",
    "sphinx_copybutton",
]

# Napoleon settings for NumPy-style docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
napoleon_attr_annotations = True

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autodoc_typehints_format = "short"

# Mock imports for modules that may not be available during doc build
autodoc_mock_imports = ["numpy", "pandas"]

# Intersphinx mappings for cross-referencing external documentation
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}

# Support both Markdown and RST files
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# MyST configuration for Markdown files
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "linkify",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 3

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "snps"
copyright = "2019, Andrew Riha"
author = "Andrew Riha"

# Version info
version = snps.__version__
release = version

language = "en"

# List of patterns to ignore
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Pygments style
pygments_style = "sphinx"

# -- Options for HTML output ----------------------------------------------

html_theme = "furo"

html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#2166ac",
        "color-brand-content": "#2166ac",
    },
    "dark_css_variables": {
        "color-brand-primary": "#92c5de",
        "color-brand-content": "#92c5de",
    },
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "top_of_page_button": "edit",
    "source_repository": "https://github.com/apriha/snps/",
    "source_branch": "main",
    "source_directory": "docs/",
}

html_static_path = []

# -- Options for HTMLHelp output ------------------------------------------

htmlhelp_basename = "snpsdoc"

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {}

latex_documents = [
    (master_doc, "snps.tex", "snps Documentation", "Andrew Riha", "manual")
]

# -- Options for manual page output ---------------------------------------

man_pages = [(master_doc, "snps", "snps Documentation", [author], 1)]

# -- Options for Texinfo output -------------------------------------------

texinfo_documents = [
    (
        master_doc,
        "snps",
        "snps Documentation",
        author,
        "snps",
        "Tools for reading, writing, merging, and remapping SNPs.",
        "Miscellaneous",
    )
]


# Document __init__ methods
def skip(app, what, name, obj, would_skip, options):
    if name == "__init__":
        return False
    return would_skip


def setup(app):
    app.connect("autodoc-skip-member", skip)
