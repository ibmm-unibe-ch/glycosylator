# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "glycosylator"
copyright = "2023, Noah Kleinschmidt"
author = "Noah Kleinschmidt, Boris Gusev, and Thomas Lemmin"
release = "4.4.0"


import plotly.io as pio

pio.renderers.default = "sphinx_gallery"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx.ext.autodoc",
    "sphinx_design",
    "nbsphinx",
    "sphinxmermaid",
    # "sphinx_panels",
    # "sphinx_gallery.gen_gallery",
]

# sphinx_gallery_conf = {
#     "examples_dirs": "../examples",  # path to your example scripts
#     "gallery_dirs": "./examples",  # path to where to save gallery generated output
# }

napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_admonition_for_notes = True
nbsphinx_allow_errors = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
# html_logo = "_resources/logo_small.svg"
html_static_path = ["_static"]

html_theme_options = {
    "pygment_light_style": "tango",
    "pygment_dark_style": "monokai",
    "logo": {
        "image_light": "_resources/logo_small.svg",
        "image_dark": "_resources/logo_small_dark.svg",
    },
}
