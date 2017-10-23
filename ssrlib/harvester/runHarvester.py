import os
import subprocess
basedir = os.path.abspath(os.path.dirname(__file__))


def run_harverter(config):
    str_result_dir = os.path.join(config['report'],'structure')
    outdir = os.path.join(config['report'],'strharvest')
    subprocess.call(['python', os.path.join(basedir, 'structureHarvester/structureHarvester.py'),'--dir=%s' %str_result_dir, '--out=%s' %outdir, '--evanno', '--clumpp'])


