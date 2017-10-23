#!/usr/bin/env python
import os
import re
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
import subprocess
from .. import utils
import logging

basedir = os.path.abspath(os.path.dirname(__file__))


def read_evanno(config):
    fevanno = os.path.join(config['report'],'strharvest/evanno.txt')
    utils.check_file(fevanno)
    dic_evanno = defaultdict(list)
    for line in file(fevanno):
        if line.startswith('#') or re.match(r'^\s+$',line):
            continue
        arr = line.strip().split()
        dic_evanno['K'].append(int(arr[0]))
        dic_evanno['Mean_LnP'].append(float(arr[2]))
        dic_evanno['stdev'].append(float(arr[3]))
        dic_evanno["Ln'(K)"].append(arr[4])
        dic_evanno["Ln''(K)"].append(arr[5])
        dic_evanno["Delta(K)"].append(arr[6])
    return dic_evanno

def plot_LK(x_vector, y_vector,err,outdir):
    fig, ax = plt.subplots()
    ax.errorbar(x_vector, y_vector, yerr=err,fmt='o')
    ax.set_ylabel('Mean of est. Ln prob. of data')
    ax.set_xlabel('K')
    ax.set_title('L(K) (mean +- SD)')
    fig.savefig('%s/meanLnPlot.png'%outdir,dpi=500)
    plt.close(fig)
    logging.info('--- %s/meanLnPlot.png plotted---',outdir)

def plot_line(x_vector, y_vector, title, ylabel, filename,outdir,style='scatter'):
    fig, ax = plt.subplots()
    if style == 'scatter':
        ax.scatter(x_vector,y_vector)
    else:
        ax.plot(x_vector,y_vector,marker='o')
    ax.set_xlim(0,len(y_vector)+2)
    ax.set_xlabel('K')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.savefig('%s/%s.png'%(outdir,filename), dpi=500)
    plt.close(fig)
    logging.info('--- %s/%s.png plotted ---',outdir, filename)

def main_plot(dic_evanno,config):
    report = config['report']
    K = dic_evanno['K']
    Mean_LnP = dic_evanno['Mean_LnP']
    stdev = dic_evanno['stdev']
    Ln_K_prime = dic_evanno["Ln'(K)"][1:]
    Ln_K_2prime = dic_evanno["Ln''(K)"][1:-1]
    Ln_K_deltaK = dic_evanno["Delta(K)"][1:-1]
    name_of_K_prime = 'Rate of change of the likelihood distribution(mean)'
    name_of_K_2prime= 'Absolute value of the 2nd order rate of change of the likelihood distribution(mean)'
    name_of_deltaK = 'Deltak = mean(|L''(K)|)/sd(L(K))'
    plot_LK(K, Mean_LnP,stdev,report)
    plot_line(K[1:], Ln_K_prime,name_of_K_prime,"L'(K)",'LnPK',report)
    plot_line(K[1:-1], Ln_K_2prime,name_of_K_2prime,"L''(K)",'LnPPK',report)
    plot_line(K[1:-1], Ln_K_deltaK,name_of_deltaK,"Delta K",'deltaK',report,style='line')

def get_est_K(dic_evanno):
    K = dic_evanno['K'][1:-1]
    Ln_K_deltaK = dic_evanno["Delta(K)"][1:-1] # first and last element of this list is 'NA', which is meanless
    data_delta = map(lambda x: float(x), Ln_K_deltaK)
    tmpdic = dict(zip(K, data_delta)) # a dict of K: Ln_K_deltaK
    # est_K = data_delta.index(max(data_delta))+2
    # reverse sort of K by Ln_K_deltaK, so that first K of the sorted list has greatest Ln_K_deltaK
    decrease_lnk = sorted(tmpdic.keys(), key=lambda k: tmpdic[k], reverse=True)
    est_K = decrease_lnk[0]
    logging.info('[NOTICE] BEST K ESTIMATED WAS %s', est_K)
    return est_K

def sample_num(config):
    fspagedi = os.path.join(config['data'], 'spagedi.txt')
    line1 = open(fspagedi).readline().strip().split()
    nsam, npop, _, loci = line1[0:4]
    return loci, nsam, npop

def write_clumpp_parafile(config,est_K,dic_evanno, nsam, npop,datatype=0):
    tmpdir = os.path.join(config['routine'],'tmp')
    harvest_dir = os.path.join(config['report'],'strharvest')
    report = config['report']
    r = config.get('replicate') if config.get('replicate') else 3
    c = nsam if datatype == 0 else npop
    outfile = 'K%s.indq'% est_K if datatype == 0 else 'K%s.popq'% est_K
    popfile = os.path.join(harvest_dir, 'K%s.popfile' % est_K)
    if datatype == 1:
        popfile = polish_popfile(harvest_dir, est_K, r)
    logging.info('Prepare clumpp paramfile for %s', outfile)
    print 'Prepare clumpp paramfile for %s' % outfile
    with open('%s/clumpparam' % tmpdir,'w') as output:
        output.write('DATATYPE %s\n' % datatype)
        output.write('INDFILE %s/K%s.indfile\n' %(harvest_dir,est_K))
        output.write('POPFILE %s\n' % popfile)
        output.write('OUTFILE %s/%s\n' %(report,outfile))
        output.write('MISCFILE %s/K%s.miscfile\n' %(report,est_K))
        output.write('K %s\n' % est_K)
        output.write('C %s\n' % c)
        output.write('R %s\n' % r)
        output.write('M 1\n')
        output.write('W 1\n')
        output.write('S 2\n')
        output.write('GREEDY_OPTION 2\n')
        output.write('REPEATS 1000\n')
        output.write('PRINT_PERMUTED_DATA 1\n')
        output.write('PRINT_PERMUTED_DATA 1\n')
        output.write('OVERRIDE_WARNINGS 1\n')
        output.write('ORDER_BY_RUN 1\n')
        output.write('PRINT_EVERY_PERM 0\n')
        output.write('PRINT_RANDOM_INPUTORDER 0\n')
    return '%s/clumpparam' % tmpdir

def polish_popfile(harvest_dir, est_K, replicate):
    polished = os.path.join(harvest_dir, 'K%s.check.popfile' % est_K)
    current_line = ''
    count  = 0
    with open(polished, 'w') as oh:
        for line in file(os.path.join(harvest_dir, 'K%s.popfile' % est_K)):
            arr = line.strip().split()
            if count <= int(replicate):
                if not arr:
                    continue
                count += 1
                line_info = [arr[0], arr[-1]]
                if count == 1:
                    current_line = line_info
                    line1 = line
                    continue
                if line_info == current_line:
                    continue
                else:
                    oh.write(line1)
                    count = int(replicate) + 1
            oh.write(line)
    return polished

def run_clumpp(config,est_K, dic_evanno, nsam, npop):
    paramfile = write_clumpp_parafile(config, est_K, dic_evanno, nsam, npop)
    subprocess.call([os.path.join(basedir, 'CLUMPP'), paramfile])
    paramfile = write_clumpp_parafile(config, est_K, dic_evanno, nsam, npop,datatype=1)
    subprocess.call([os.path.join(basedir, 'CLUMPP'), paramfile])

def distruct_prepare(config,est_K,nsam, npop):
    """Prepare drawparams for distruct."""
    import math
    report = config['report']
    tmpdir = os.path.join(config['routine'],'tmp')
    subprocess.call('cp -r /home/wuj/bin/SSRPop/distruct %s' %tmpdir,shell=True)
    subprocess.call('cp %s/K%s.* %s/distruct' %(report,est_K,tmpdir),shell=True)
    popfile = os.path.join(tmpdir, 'distruct', 'K%s.popq' % est_K)
    sort_pops(config, popfile)
    bin_width = min(25/math.log(float(nsam))**2, 5.0)
    with open('%s/distruct/drawparams' %tmpdir,'a') as output:
        output.write('#define INFILE_POPQ   %s\n' % popfile)
        output.write('#define INFILE_INDIVQ   K%s.indq\n' % est_K)
        output.write('#define INFILE_LABEL_BELOW infile_sort.names\n')
        output.write('#define INDIVWIDTH %s\n' % bin_width)
        output.write('#define OUTFILE   K%s.ps\n' % est_K)
        output.write('#define K   %s\n' % est_K)
        output.write('#define NUMPOPS   %s\n' % npop)
        output.write('#define NUMINDS   %s\n' % nsam)
    return

def sort_pops(config, popfile):
    tmpdir = os.path.join(config['routine'],'tmp')
    datadir= config['data']
    tmpdic = {}
    for line in file(popfile):
        line = line.replace(':', '')
        arr = line.split()
        if not arr:
            continue
        tmpdic[arr[0]] = map(float, arr[1:])
    decrease = sorted(tmpdic.keys(), key=lambda k: tmpdic[k][0], reverse=True)

    tmpdic2 = {}
    for line in file(os.path.join(datadir, 'infile.names')):
        arr = line.split()
        if not arr:
            continue
        tmpdic2[arr[0]] = line
    with open(os.path.join(tmpdir, 'distruct', 'infile_sort.names'), 'w') as oh:
        for n in decrease:
            line = tmpdic2.get(n, '_\n')
            oh.write(line)
    return

def run_distruct(config,est_K):
    report = config['report']
    tmpdir = os.path.join(config['routine'],'tmp')
    os.chdir('%s/distruct' %tmpdir)
    subprocess.call('./distructLinux1.1', shell=True)
    subprocess.call('ps2pdf K%s.ps %s/K%s.pdf' %(est_K,report,est_K),shell=True)


def str_visualize(config):
    dic_evanno = read_evanno(config)
    main_plot(dic_evanno,config)
    est_K = get_est_K(dic_evanno)
    loci, nsam, npop = sample_num(config)
    run_clumpp(config, est_K, dic_evanno, nsam, npop)
    distruct_prepare(config, est_K, nsam, npop)
    run_distruct(config, est_K)



if __name__ == '__main__':
    configfile = sys.argv[1]
    config = utils.read_config(configfile)
    str_visualize(config)



