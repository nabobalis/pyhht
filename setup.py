import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "pyhht",
    version = "0.0.1",
    author = "Nabil Freij",
    description = ("A package for the manipulation of SST data."),
    license = "BSD",
    url = "https://github.com/nabobalis/pyhht",
    packages=['pyhht'],
    install_requires=['numpy','scipy'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)