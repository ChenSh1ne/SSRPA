#!/usr/bin/python
import re
import os
import sys
import logging
import utils
import subprocess
from collections import namedtuple
from multiprocessing import Pool



def cal_sample_info(config):
    data_ssr = os.path.join(config['data'],'ssr_structure')
    utils.check_file(data_ssr) # check if file data_ssr exists
    ploidy = int(config['ploidy'])
    samples = []
    # A namedtuple to contain parameters
    project_info = namedtuple('project_info',['samples','loci','filestr'])
    count = 0
    for line in file(data_ssr):
        count += 1
        if count == 1:
            header = line.strip().split()
        else:
            samp_name = line.strip().split()[0]
            # samples.add(line.strip().split()[0])
            if samp_name not in samples:
                samples.append(samp_name)
    header = filter(lambda x: x.strip(),header)
    loci_num = len(header)
    samp_num = len(samples)
    INFO = project_info(samples=samp_num,loci=loci_num,filestr=data_ssr)
    write_names(config, samples)
    if (count - 1)%ploidy != 0:
        #logging...
        logging.warning('sample num incorrect in %s', data_ssr)
        print '[WARNING] sample num incorrect in %s' % data_ssr
    return INFO

def write_names(config,samples):
    data_dir = config['data']
    with open("%s/samples.num" %data_dir,'w') as out:
        for i, j in enumerate(samples):
            out.write('%s %s\n' %(i+1,j))

def write_param_files(PROJ_INFO,config):
    print 'writing structure param files...'
    tmp_dir = os.path.join(config['routine'],'tmp')
    utils.make_dir(tmp_dir)
    with open('%s/mainparams' %tmp_dir,'w') as handle:
        handle.write('#define INFILE %s\n' %PROJ_INFO.filestr)
        handle.write('#define NUMINDS %s\n' %PROJ_INFO.samples)
        handle.write('#define NUMLOCI %s\n' %PROJ_INFO.loci)
        handle.write('#define PLOIDY %s\n' %config['ploidy'])
        handle.write('#define MISSING -9\n')
        handle.write('#define ONEROWPERIND 0\n')
        handle.write('#define LABEL 1\n')
        handle.write('#define POPDATA 1\n')
        handle.write('#define POPFLAG 0\n')
        handle.write('#define LOCDATA 0\n')
        handle.write('#define PHENOTYPE 0\n')
        handle.write('#define BURNIN 10000\n')
        handle.write('#define NUMREPS 50000\n')
        handle.write('#define EXTRACOLS 0\n')
        handle.write('#define MARKERNAMES 1\n')
        handle.write('#define RECESSIVEALLELES 0\n')
        handle.write('#define MAPDISTANCES 0\n')
        handle.write('#define PHASED 0\n')
        handle.write('#define PHASEINFO 0\n')

    with open('%s/extraparams' %tmp_dir,'w') as ehandle:
        ehandle.write(
        '''#define NOADMIX 0
        #define LINKAGE
        #define USEPOPINFO 0
        #define LOCPRIOR 0
        #define ONEFST 0
        #define INFERALPHA 1
        #define POPALPHAS 0
        #define ALPHA 1.0
        #define INFERLAMBDA 0
        #define POPSPECIFICLAMBDA 0
        #define LAMBDA''')
    print 'structure paramfile OK.'
    return None

def run_structure(mainparams, extraparams,k,config,run_per_k=3):
    print 'run structure'
    outdir = os.path.join(config['report'],'structure')
    utils.make_dir(outdir)
    output = os.path.join(outdir,'output_k%s' %k)
    utils.check_file(mainparams)
    utils.check_file(extraparams)
    if re.match(r'\d+',config.get('replicate','')):
        print 'matched'
        run_per_k = int(config['replicate'])
        logging.info('Got replicate %s for each K', run_per_k)
    for r in range(run_per_k):
        logging.info('Running structure at K %s for round %s',k,r)
        subprocess.call("/home/pub/software/Structure/bin/structure -m %s -e %s -K %d -o %s_run%s" %(mainparams, extraparams,k, output,r),shell=True)
    return None

def run_in_parellel(config,processes=20):
    PINFO = cal_sample_info(config)
    print PINFO
    write_param_files(PINFO, config)
    mainparams = os.path.join(config['routine'],'tmp','mainparams')
    extraparams = os.path.join(config['routine'],'tmp','extraparams')
    maxpop =int(config['maxpop'])
    pool = Pool(processes)
    for k in xrange(1,maxpop+1):
        pool.apply_async(run_structure,(mainparams,extraparams,k,config))
    pool.close()
    pool.join()
    return None


if __name__ == '__main__':
    configfile = sys.argv[1]
    config = utils.read_config(configfile)
    run_in_parellel(config)


