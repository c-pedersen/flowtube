# Configuration file for the Sphinx documentation builder.
import os
import sys
import re

# Make the project root importable
sys.path.insert(0, os.path.abspath("../.."))

extensions = [
    "sphinx.ext.napoleon",  # For Google/NumPy style docstrings
    "sphinx.ext.viewcode",  # Add source code links
    "sphinx_autodoc_typehints",
    "sphinx_mdinclude",
    "autoapi.extension",
    "nbsphinx",
]

project = "flowtube"
copyright = "2025, Corey Pedersen"
author = "Corey Pedersen"

with open("../../pyproject.toml") as f:
    text = f.read()
match = re.search(r'version\s*=\s*"([^"]+)"', text)
if match:
    release = match.group(1)
else:
    raise ValueError("Version not found in pyproject.toml")

version = release.rsplit(".", 1)[0]  # Short version (e.g., "1.2")

autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
]
autoapi_dirs = ["../../flowtube/"]

napoleon_google_docstring = True
napoleon_numpy_docstring = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
