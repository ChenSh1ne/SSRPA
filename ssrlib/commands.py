#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""Population Analysis With SSR Markers"""

__version__ = '0.1.0.0210_alpha'

import os
import sys
import re
import argparse
import logging
import subprocess
import random
import multiprocessing

from . import runStructure, numericalize, utils
from .clumpp import autoCLUMPP
from .harvester import runHarvester


AP = argparse.ArgumentParser(
                description="SSRPA, a toolkit for population analysis with SSR markers."
        )
AP_subparsers = AP.add_subparsers(
                help="Sub-commands (use with -h for more info)"
        )

# batch -----------------------------------------------------------------------

def _batch_command(args):
    pass


# structure -------------------------------------------------------------------

def _structure_command(args):
    """Perform Structure Analysis."""
    strfile = args.strfile
    if not os.path.exists(strfile):
        sys.exit('ERROR: File %s not found. Please perform [formatting] first.' % strfile)
    config = utils.read_config(args.cfg)
    runStructure.run_in_parellel(config)
    runHarvester.run_harverter(config)
    autoCLUMPP.str_visualize(config)
    return

P_struct = AP_subparsers.add_parser('structure', help=_structure_command.__doc__)
P_struct.add_argument('-file--struct', metavar='structure format file', dest='strfile', required=True)
P_struct.add_argument('-cfg', metavar='config file', required=True)
P_struct.set_defaults(func=_structure_command)

# formatting ------------------------------------------------------------------

def _format_command(args):
    """File format transformer. Subpopulation or provenance information
    should be provided if they were explicit."""

    config = utils.read_config(args.cfg)
    numericalize.data_prepare(config)
    logging.info('raw ssr data transfered to formatted data')
    return

P_format = AP_subparsers.add_parser('formatting', help=_format_command.__doc__)
P_format.add_argument('-cfg', metavar='config file', required=True)
P_format.set_defaults(func=_format_command)

# pholygenetic ----------------------------------------------------------------

# genetic diversity -----------------------------------------------------------

# _____________________________________________________________________________
# Shim for command-line execution

def parse_args(args=None):
    """Parse the command line."""
    return AP.parse_args(args=args)






