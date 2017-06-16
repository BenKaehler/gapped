#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2017--, Ben Kaehler
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

setup(
    name='gapped',
    version='0.0.0-dev',
    license='BSD-3-Clause',
    packages=find_packages(),
    install_requires=['cogent'],
    author="Ben Kaehler",
    author_email="kaehler@gmail.com",
    description="Monkey patches to make cogent codon models work with gaps",
    url="https://github.com/BenKaehler/gapped"
)
