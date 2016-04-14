#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##
##  Copyright 2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##############################################################################


#from distutils.core import setup
from setuptools import setup

setup(
    name="inphest",
    version="0.1.0",
    author="Jeet Sukumaran",
    author_email="jeetsukumaran@gmail.com",
    packages=[
        "inphest",
        # "test"
        ],
    # package_data = {
    #         'inphest': ['R/*.R'],
    #         'inphest': ['libexec/*.R'],
    #     },
    include_package_data = True,
    scripts=[
    #         "bin/inphest-classify.py",
    #         "bin/inphest-profile-trees.py",
            "bin/inphest-simulate.py",
    #         "bin/inphest-summarize.py",
    #         "bin/inphest-generate-data-files-from-tip-labels.py",
            ],
    url="http://pypi.python.org/pypi/inphest/",
    license="LICENSE.txt",
    description="A Project",
    long_description=open("README.md").read(),
    # install_requires=[ ],
)
