"""Hook for the mkdocs `macros` plugin: defines the macros the documentation pages can call.

`output_variable_catalogue()` writes the output variable catalogue page. It is rebuilt from
the source files (the YAML catalogue and the Fortran bindings) each time the documentation is
built, so the page cannot go out of date. See `catalogue_docs.py`.
"""
import sys
from pathlib import Path

# The macros plugin loads this file by path, so its own directory is not on the import path.
# Add it, otherwise the import below fails and the plugin reports "no main module" instead.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from catalogue_docs import render_catalogue  # noqa: E402


def define_env(env):
    """Register the macros. `env.project_dir` is the directory that contains mkdocs.yml."""

    @env.macro
    def output_variable_catalogue():
        return render_catalogue(env.project_dir)
