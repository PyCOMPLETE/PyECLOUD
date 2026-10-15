"""Documentation builds independently of the simulation dependencies."""

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
project = "PyECLOUD"
author = "PyECLOUD contributors"
copyright = "CERN,  Geneva  2011-Present"
version_module = ast.parse((ROOT / "PyECLOUD" / "_version.py").read_text())
release = next(
    ast.literal_eval(node.value)
    for node in version_module.body
    if isinstance(node, ast.Assign)
    and any(isinstance(target, ast.Name) and target.id == "__version__"
            for target in node.targets)
)
version = release

extensions = ["myst_parser", "sphinx.ext.mathjax"]
source_suffix = {".md": "markdown"}
myst_enable_extensions = ["dollarmath"]
myst_heading_anchors = 3
exclude_patterns = ["_build", ".venv-docs", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
html_title = "PyECLOUD documentation"
html_logo = "pyecl_logo_darkbg.png"
html_theme_options = {"logo_only": True, "navigation_depth": 3}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sourcelink = True
linkcheck_ignore = [r"http://localhost.*", r"http://127\.0\.0\.1.*"]
