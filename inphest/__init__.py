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

import os

__project__ = "Inphest"
__version__ = "0.1.0"
__inphest_revision__ = None
__inphest_description__ = None

from inphest.simulate import repeat_run as run

INPHEST_HOME = os.path.dirname(os.path.abspath(__file__))

def libexec_filepath(filename):
    try:
        import pkg_resources
        filepath = pkg_resources.resource_filename("inphest", "libexec/{}".format(filename))
    except:
        filepath = os.path.normpath(os.path.join(os.path.dirname(__file__), filename))
    return filepath

def revision():
    global __inphest_revision__
    if __inphest_revision__ is None:
        from dendropy.utility import vcsinfo
        try:
            try:
                __homedir__ = os.path.dirname(os.path.abspath(__file__))
            except IndexError:
                __homedir__ = os.path.dirname(os.path.abspath(__file__))
        except OSError:
            __homedir__ = None
        except:
            __homedir__ = None
        __inphest_revision__ = vcsinfo.Revision(repo_path=__homedir__)
    return __inphest_revision__

def description():
    global __inphest_description__
    if __inphest_description__ is None:
        inphest_revision = revision()
        if inphest_revision.is_available:
            revision_text = " ({})".format(inphest_revision)
        else:
            revision_text = ""
        __inphest_description__  = "{} {}{}".format(__project__, __version__, revision_text)
    return __inphest_description__

