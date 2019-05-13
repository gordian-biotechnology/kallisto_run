#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
'''
with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()
'''
requirements = [
    "numpy>=1.0, <1.15",
    "torch>=0.4.1",
    "matplotlib>=2.0",
    "scikit-learn>=0.18, <0.20.0",
    "scipy>=1.2",
    "h5py>=2.8",
    "pandas>=0.24",
    "loompy>=2.0",
    "tqdm >= 4",
    "anndata >= 0.6",
    "xlrd >= 1.0",
    "jupyter>=1.0.0",
    "nbconvert>=5.4.0",
    "nbformat>=4.4.0",
    "ipython>=7",
    "umap-learn>=0.3.7",
    "seaborn>=0.9.0",
    "leidenalg>=0.7.0",
    "python-igraph>=0.7.1",
    "hyperopt>=0.1.2",
    "scanpy>=1.4"
]


author = 'Ian Driver'

setup(
    author=author,
    author_email='ian@gordian.bio',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    description="Kallisto single cell bus with scvi/scanpy analysis",
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='kallisto single-cell',
    name='kallisto_run',
    packages=find_packages(),
    url='https://github.com/iandriver/kallisto_run',
    package_dir={'kallisto_run':
                 'kallisto_run'},
    entry_points={
          'console_scripts': ['kallisto_run = kallisto_run.kallisto_run:main']
        },
    version='0.2.0',
    zip_safe=False,
)
