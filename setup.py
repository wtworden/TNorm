
from __future__ import print_function

import os
import platform

import setuptools

from sage.version import version as sage_version


with open("README.md", "r") as fh:
    long_description = fh.read()

this_directory = os.path.dirname(__file__)
source_directory = os.path.join(this_directory, 'tnorm')
exec(open(os.path.join(source_directory, 'version.py')).read())  # Load in the variable __version__.

mac_ver = [int(x) for x in platform.mac_ver()[0].split('.')]

if mac_ver[0] == 10 and mac_ver[1] <= 13:
    if float(sage_version) >= 9:
        dependencies = ['networkx>=2.4', 'snappy>=3.0', 'sageRegina>=6.0.1']
    else:
        dependencies = ['snappy','sageRegina==5.1.5']

else:
    if float(sage_version) >= 9:
        dependencies = ['networkx>=2.4', 'snappy>=3.0', 'regina>=6.1.0.dev0']
    else:
        dependencies = ['snappy','sageRegina==5.1.5']


setuptools.setup(
    name='tnorm',  
    version=__version__,
    author="William Worden",
    author_email="wtworden@gmail.com",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wtworden/TNorm",
    packages=['tnorm', 'tnorm.kernel', 'tnorm.GUI', 'tnorm.test', 'tnorm.utilities', 'tnorm.extras'],
    package_data={
        'tnorm.GUI': ['images/spinpoly/*','images/icon/*','data/*'],
        'tnorm.utilities': ['data/x3d_templates/*', 'data/x3dom/*'],
        },
    install_requires=dependencies,
    classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
    ],
 )


## sage -python setup.py sdist bdist_wheel --universal
## twine upload dist/*

## to install dev version, cd to TNorm directory, then:
## sage -pip install -e .