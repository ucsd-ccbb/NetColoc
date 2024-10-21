#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open(os.path.join('netcoloc', '__init__.py')) as ver_file:
    for line in ver_file:
        line = line.rstrip()
        if line.startswith('__version__'):
            version = re.sub("'", "", line[line.index("'"):])
        elif line.startswith('__author__'):
            author = re.sub("'", "", line[line.index("'"):])
        elif line.startswith('__email__'):
            email = re.sub("'", "", line[line.index("'"):])

requirements = [
    'click>=6.0',
    'ndex2',
    'networkx>=2.0,<3.0',
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
    'ipywidgets',
    'cdapsutil',
    'obonet'
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='netcoloc',
    version=version,
    description="Analyze network colocalization of gene sets",
    long_description=readme + '\n\n' + history,
    author=author,
    author_email=email,
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
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13'
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=['wheel']
)
