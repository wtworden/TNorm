
from __future__ import print_function

import os

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

this_directory = os.path.dirname(__file__)
source_directory = os.path.join(this_directory, 'tnorm')
exec(open(os.path.join(source_directory, 'version.py')).read())  # Load in the variable __version__.


setuptools.setup(
    name='tnorm',  
    version=__version__,
    author="William Worden",
    author_email="wtworden@gmail.com",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wtworden/TNorm",
    packages=['tnorm', 'tnorm.kernel', 'tnorm.GUI'],
    package_data={
        'tnorm.GUI': ['images/spinpoly/*','images/icon/*','data/*'],
        'tnorm': ['data/x3d_templates/*', 'data/x3dom/*'],
        },
    classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
    ],
 )