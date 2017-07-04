#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import sphinx_rtd_theme

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

extensions = ['breathe', 'sphinx.ext.mathjax']
breathe_projects = { 'exaFMM': '../xml'}

source_suffix = '.rst'
master_doc = 'index'
project = 'exaFMM'

exclude_patterns = []
highlight_language = 'c++'
pygments_style = 'sphinx'
todo_include_todos = False
htmlhelp_basename = 'exafmmdoc'
