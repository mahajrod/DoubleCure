__author__ = 'mahajrod'
import os
from setuptools import setup

dependencies = ["pandas", "routoolpa"]

setup(name='doublecure',
      version='0.1a',
      author='mahajrod',
      install_requires=dependencies,
      scripts=["workflow/scripts/unify_agps.py", "workflow/scripts/correct_breakpoints.py"],
      long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),)
