#!/usr/bin/env python
import os
import re
import time
import signal
import random
import logging
import subprocess

def make_dir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def check_Xvfb():
    """Check if Xvfb is running, if not, it'll be started."""
    status = is_running('xvfb')
    if status:
        return status
    # if not is_running('xvfb'):
    else:
        server_num = random.randint(1,20) # generate a random num for Xvfb
        print 'server num:', server_num
        # ps = subprocess.Popen("Xvfb -ac :%d -screen 0 1280x1024x8 > /dev/null 2>&1" % server_num, shell=True)
        subprocess.call("Xvfb -ac :%d -screen 0 1280x1024x8 > /dev/null 2>&1" % server_num, shell=True)
        os.environ["DISPLAY"] = ":%d" % server_num
        print 'wait for process starting...'
        time.sleep(5)
        # pspid = is_running('xvfb')
        # return ps.pid

def is_running(proc):
    """Check if a specific process is running."""
    ps = os.popen("ps a")
    # ps = subprocess.Popen(["ps","a"], shell=True, stdout=subprocess.PIPE)
    for x in ps.readlines():
        if re.search(proc,x,re.I):
            pid = x.split()[0]
            return pid
    return False

def check_file(checking_file):
    if not os.path.isfile(checking_file):
        #logging ...
        sys.exit('ERROR: File %s did not exist!' % checking_file)

if __name__ == '__main__':
    proc = check_Xvfb()
    # print proc
    # time.sleep(5)
    # os.kill(proc, signal.SIGKILL)
    # print 'killed.'
