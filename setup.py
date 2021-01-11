#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools
from msa2vcf import version

setuptools.setup(
    name="msa2vcf",
    version=version.__version__,
    url="https://github.com/connor-lab/msa2vcf",

    description="",
    long_description="",

    author="Matt Bull",
    author_email="bullmj2@gmail.com",

    maintainer="Matt Bull",
    maintainer_email="bullmj2@gmail.com",

    packages=setuptools.find_packages(),
    #install_requires=requirements,

    entry_points = {
        'console_scripts': [
            'msa2vcf = msa2vcf.cli:main',
        ]
    },

    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
    ],

)
