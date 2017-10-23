#!/usr/bin/env python
import os
import re
import sys
import logging



def check_file(checking_file):
    if not os.path.isfile(checking_file):
        #logging ...
        logging.error('File %s did not exit!', checking_file)
        sys.exit('ERROR: File %s did not exist!' % checking_file)

def check_dir(checking_dir):
    if not os.path.exists(checking_dir):
        #logging ...
        logging.error('Dir %s did not exist!', checking_dir)
        sys.exit('ERROR: Dir %s did not exist!' % checking_dir)

def make_dir(dirname):
    if not os.path.exists(dirname):
        logging.info('Dir %s did not exist, automatically make dir', dirname)
        os.makedirs(dirname)


def read_config(configfile):
    print 'Reading config...'
    config = {}
    for line in file(configfile):
        if line.startswith('#'):
            continue
        if re.match(r'^\s$', line):
            continue
        line = re.sub(r'[\r\n\s+]','',line)
        arr = line.split('=')
        config[arr[0]] = arr[1]
    logging.info('configuration info loaded.')
    return config


def _help():
    # print 'SSRPA version: %s' % __version__
    if len(sys.argv) < 2:
        os.system("python %s -h" % sys.argv[0])
        #logging...
        sys.exit("\n[ERROR] NO config file provided!")
