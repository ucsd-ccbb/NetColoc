#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
from setuptools import setup, find_packages

with open(os.path.join('netcoloc', '__init__.py')) as ver_file:
    for line in ver_file:
        if line.startswith('__version__'):
            version=re.sub("'", "", line[line.index("'"):])

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'click>=6.0',
    'ndex2',
    'networkx>=2.0',
    'mygene>=3.2.2',
    'scipy>=1.5.3',
    'numpy',
    'pandas',
    'tqdm',
    'matplotlib',
    'seaborn',
    'statsmodels',
    'gprofiler-official>=1.0.0',
    'ipycytoscape',
    'ipywidgets'
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='netcoloc',
    version=version,
    description="",
    long_description=readme + '\n\n' + history,
    author="Brin Rosenthal, Sophie Liu",
    author_email='sbrosenthal@health.ucsd.edu, sol015@ucsd.edu',
    url='https://github.com/ucsd-ccbb/netcoloc',
    packages=find_packages(include=['netcoloc']),
    package_dir={'netcoloc':
                 'netcoloc'},
    entry_points={
        'console_scripts': [
            'netcoloc=netcoloc.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='netcoloc',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=['wheel']
)
