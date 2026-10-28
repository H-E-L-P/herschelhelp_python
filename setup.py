# All metadata lives in pyproject.toml. This file only exists to run
# the custom step that builds herschelhelp/data.db before packaging.
import os
import sys

from setuptools import setup
from setuptools.command.build_py import build_py

HERE = os.path.abspath(os.path.dirname(__file__))


class CustomBuild(build_py):
    """Build the database, then do the standard build."""

    def run(self):
        sys.path.insert(0, HERE)
        import database_builder
        database_builder.build_base()
        super().run()


setup(cmdclass={"build_py": CustomBuild})
