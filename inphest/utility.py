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


import decimal
import sys
import os
import logging
import inspect
import collections
import dendropy
import json

_LOGGING_LEVEL_ENVAR = "INPHEST_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "INPHEST_LOGGING_FORMAT"

TEST_DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "test", "data")

def open_output_file_for_csv_writer(filepath, append=False):
    if filepath is None or filepath == "-":
        out = sys.stdout
    elif sys.version_info >= (3,0,0):
        out = open(filepath,
                "a" if append else "w",
                newline='')
    else:
        out = open(filepath, "ab" if append else "wb")
    return out

def is_almost_equal(x, y, ndigits=7):
    return round(x-y, ndigits) == 0

def assert_in_collection(
        items,
        collection,
        at_least_one=True,
        no_more_than_one=True,
        error_type=TypeError):
    found = set()
    for item in items:
        if item in collection:
            found.add(item)
    if len(found) == 0:
        if at_least_one and no_more_than_one:
            raise error_type("Exactly one of the following must be specified: ".format(items))
        elif at_least_one:
            raise error_type("At least one of the following must be specified: ".format(items))
    elif len(found) > 1 and no_more_than_one:
        raise error_type("Only one of the following can be specified: ".format(found))
    return found

def dump_stack():
    for frame, filename, line_num, func, source_code, source_index in inspect.stack()[2:]:
        if source_code is None:
            print("{}: {}".format(filename, line_num))
        else:
            print("{}: {}: {}".format(filename, line_num, source_code[source_index].strip()))

def is_in_range(x, a, b):
    """
    Checks if $x$ is in $[a, b]$, taking into account floating-point error.
    """
    x = decimal.Decimal("{:0.8f}".format(x))
    a = decimal.Decimal("{:0.8f}".format(a))
    b = decimal.Decimal("{:0.8f}".format(b))
    return a <= x <= b

class IndexGenerator(object):

    def __init__(self, start=0):
        self.start = start
        self.index = start

    def __next__(self):
        c = self.index
        self.index += 1
        return c
    next = __next__

    def reset(self, start=None):
        if start is None:
            start = self.start
        self.index = start

class RunLogger(object):

    def __init__(self, **kwargs):
        self.name = kwargs.get("name", "RunLog")
        self._log = logging.getLogger(self.name)
        self._log.setLevel(logging.DEBUG)
        self.handlers = []
        if kwargs.get("log_to_stderr", True):
            handler1 = logging.StreamHandler()
            stderr_logging_level = self.get_logging_level(kwargs.get("stderr_logging_level", logging.INFO))
            handler1.setLevel(stderr_logging_level)
            handler1.setFormatter(self.get_default_formatter())
            self._log.addHandler(handler1)
            self.handlers.append(handler1)
        if kwargs.get("log_to_file", True):
            if "log_stream" in kwargs:
                log_stream = kwargs.get("log_stream")
            else:
                log_stream = open(kwargs.get("log_path", self.name + ".log"), "w")
            handler2 = logging.StreamHandler(log_stream)
            file_logging_level = self.get_logging_level(kwargs.get("file_logging_level", logging.DEBUG))
            handler2.setLevel(file_logging_level)
            handler2.setFormatter(self.get_default_formatter())
            self._log.addHandler(handler2)
            self.handlers.append(handler2)
        self._system = None

    def _get_system(self):
        return self._system

    def _set_system(self, system):
        self._system = system
        if self._system is None:
            for handler in self.handlers:
                handler.setFormatter(self.get_default_formatter())
        else:
            for handler in self.handlers:
                handler.setFormatter(self.get_simulation_generation_formatter())

    system = property(_get_system, _set_system)

    def get_logging_level(self, level=None):
        if level in [logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING,
            logging.ERROR, logging.CRITICAL]:
            return level
        elif level is not None:
            level_name = str(level).upper()
        elif _LOGGING_LEVEL_ENVAR in os.environ:
            level_name = os.environ[_LOGGING_LEVEL_ENVAR].upper()
        else:
            level_name = "NOTSET"
        if level_name == "NOTSET":
            level = logging.NOTSET
        elif level_name == "DEBUG":
            level = logging.DEBUG
        elif level_name == "INFO":
            level = logging.INFO
        elif level_name == "WARNING":
            level = logging.WARNING
        elif level_name == "ERROR":
            level = logging.ERROR
        elif level_name == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
        return level

    def get_default_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_simulation_generation_formatter(self):
        # f = logging.Formatter("[%(asctime)s] t = %(elapsed_time)10.6f: %(message)s")
        f = logging.Formatter("[%(asctime)s] %(simulation_time)s%(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_logging_formatter(self, format=None):
        if format is not None:
            format = format.upper()
        elif _LOGGING_FORMAT_ENVAR in os.environ:
            format = os.environ[_LOGGING_FORMAT_ENVAR].upper()
        if format == "RICH":
            logging_formatter = self.get_rich_formatter()
        elif format == "SIMPLE":
            logging_formatter = self.get_simple_formatter()
        elif format == "NONE":
            logging_formatter = self.get_raw_formatter()
        else:
            logging_formatter = self.get_default_formatter()
        if logging_formatter is not None:
            logging_formatter.datefmt='%H:%M:%S'

    def supplemental_info_d(self):
        if self._system is None or self._system.elapsed_time == 0 :
            return {
                    "simulation_time" : "Setup: ",
                    }
        else:
            return {
                    "simulation_time" : "[t = {:13.6f}] ".format(self._system.elapsed_time),
                    }

    def debug(self, msg):
        self._log.debug("[DEBUG] {}".format(msg), extra=self.supplemental_info_d())

    def info(self, msg):
        self._log.info(msg, extra=self.supplemental_info_d())

    def warning(self, msg):
        self._log.warning(msg, extra=self.supplemental_info_d())

    def error(self, msg):
        self._log.error(msg, extra=self.supplemental_info_d())

    def critical(self, msg):
        self._log.critical(msg, extra=self.supplemental_info_d())

