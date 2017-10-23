#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import re
import logging
import utils
import xlrd
import xlsxwriter
from collections import defaultdict

def load_workbook(str_xlsx):
    return xlrd.open_workbook(str_xlsx, on_demand = True)

def read_depth(workbook):
    sample_depth = defaultdict(lambda: defaultdict(int))
    sheet_depth = workbook.sheet_by_index(0)
    nrow = sheet_depth.nrows
    ncol = sheet_depth.ncols
    header = sheet_depth.row_values(0)
    targets= sheet_depth.col_values(0)
    for row in range(1,nrow):
        for col in range(3,ncol):
            value = sheet_depth.cell(row, col).value
            if re.match(r'\d+',str(value)):
                sample_depth[header[col]][targets[row]] += int(value) # header[col] = sample, targets[row] = target
    return sample_depth, set(targets[1:]), header[3:]

def read_typing(workbook):
    dic_type = defaultdict(lambda: defaultdict())
    target_str = {}
    typing = workbook.sheet_by_name('TypingAbsolute')
    nrow = typing.nrows
    header = typing.row_values(0)
    # samples= sheet_depth.col_values(0)
    # targets= sheet_depth.col_values(1)
    for row in range(1, nrow):
        rowvals = typing.row_values(row)
        target_str[rowvals[1]] = rowvals[2]
        dic_type[rowvals[0]][rowvals[1]] = parse_row(rowvals, header)  # rowvals[0] = sample, rowvals[1] = target
    return dic_type, target_str

def parse_row(rowvals, header):
    alleles = []
    for i, j in enumerate(rowvals[4:]):
        if re.match(r'\d+',str(j)):
            motif = header[4:][i]
            tmp = [int(motif)] * int(j)
            alleles.extend(tmp)
    return alleles

def write_alleles_xlsx(dic_type, target_str, dic_depth):
    # ploidy = int(config['ploidy'])
    # outdir = config['data']
    ploidy = 6
    workbook = xlsxwriter.Workbook('STRAlleles.xlsx')
    worksheet1 = workbook.add_worksheet('Mark')
    worksheet2 = workbook.add_worksheet('Mark2')
    row = 0
    col = 0
    title = ['allele%s' %i for i in range(1, ploidy + 1)]
    title = ['samplename', 'mark', 'ref'] + title
    for i in title:
        worksheet1.write(row, col, i)
        worksheet2.write(row, col, i)
        col += 1
    row += 1
    new_dic = dict_transpose(dic_type)

    row_tmp1 = []
    row_tmp2 = []
    for target in new_dic:
        for sample in new_dic[target]:
            depth = dic_depth[sample][target]
            row_tmp1.extend([sample, target, target_str[target]])
            row_tmp2 = row_tmp1[::]
            alleles = new_dic[target][sample]
            if alleles and depth >= 100:
                row_tmp1.extend(alleles)
                for i in set(alleles):
                    row_tmp2.append(i)
            else:
                row_tmp1.append(-9)
                row_tmp2.append(-9)
            for i,item in enumerate(row_tmp1):
                worksheet1.write(row, i, item)
            for i, j in enumerate(row_tmp2):
                worksheet2.write(row, i, j)
            row_tmp1 = []
            row_tmp2 = []
            row += 1


def parse_type(genotype):
    """Treat combination of repeats as genotype of an individule."""
    if genotype.strip():
        arr = genotype.strip().split('/')
        rep_pairs = map(lambda x: x.split('*'),arr)
        return set([pair[1] for pair in rep_pairs])
    else:
        return genotype.strip()


def cal_reps(genotype):
    """Turn genotype of an individule at a loci into a list of motif repeats,
    regardless of redundant."""
    if genotype.strip():
        arr = genotype.split('/')
        rep_pairs = map(lambda x: x.split('*'),arr)
        return [int(pair[1]) for pair in rep_pairs for i in range(int(pair[0]))]
    else:
        return [-9]

def read_popinfo(fpopinfo):
    """Read subpop information from a external file and save it in a dictionary"""
    dic_pop = {}
    subpop  = set()
    count = 0
    for line in file(fpopinfo):
        count += 1
        if re.match(r'^\s+$',line) or count == 1:
            continue
        arr = line.strip().split()
        dic_pop[arr[0]] = arr[1]
        subpop.add(arr[1])
    return dic_pop, subpop


def load_data(fssr, func=parse_type):
    """ Diferent `func` for different purpose,
        `parse_type` for binary_data,
        `cal_reps` for structure.
    """
    logging.debug('Data loading method was %s', func)
    utils.check_file(fssr)
    dic_type = defaultdict(lambda: defaultdict())
    targets = []
    count = 0
    for line in file(fssr):
        count += 1
        if count == 1:
            header = line.strip().split('\t')
        else:
            larr = line.strip().split('\t')
            targets.append(larr[0])
            for i, j in enumerate(larr[1:]):
                j_transed = func(j)
                dic_type[larr[0]][header[i+1]] = j_transed #dic_type['target']['sample']
    #targets.insert(0,'Sample')
    return dic_type,header,targets

def allels_of_loci(dic_type,ssr):
    """List all the allels of a single loci"""
    allels_of_loci = []
    for sample in dic_type[ssr]:
        geno = dic_type[ssr][sample]
        if not geno: # In case geno is empty
            continue
        for i in geno: # "i" means allel (ie. repeat times) here
            if i not in allels_of_loci:
                allels_of_loci.append(i)
    allels_of_loci.sort()
    return allels_of_loci

def binary_transfer(all_allels, sample_allels):
    """Transfer genotype to binary data"""
    tarr = []
    for i in all_allels:
        if i in sample_allels:
            tarr.append('1')
        else:
            tarr.append('0')
    return ''.join(tarr)

def to_binary_data(dic_type):
    dic_num = defaultdict(lambda: defaultdict())
    exists_types = []
    for ssr in dic_type:
        # To get all the allels of this ssr
        allels = allels_of_loci(dic_type,ssr)
        for sample in dic_type[ssr]:
            geno = dic_type[ssr][sample]
            binary_data = binary_transfer(allels, geno)
            dic_num[sample][ssr] = binary_data
    return dic_num

def write_out(dic_num,targets):
    with open('ssr_num.txt','w') as out:
        out.write('\t'.join(targets)+'\n')
        for sample in dic_num:
            out.write(sample)
            for ssr in targets[1:]:
                out.write('\t'+str(dic_num[sample][ssr]))
            out.write('\n')
    print 'Result put down.'
    return None

def dict_transpose(dic_motif):
    new_dic = defaultdict(lambda: defaultdict())
    for ssr in dic_motif:
        for sample in dic_motif[ssr]:
            new_dic[sample][ssr] = dic_motif[ssr][sample]
    return new_dic

def spagedi_genotype(l_type):
    """Turn the allel list into spagedi genotype"""
    if len(l_type) > 1:
        return '/'.join(map(lambda x: str(x),l_type))
    else:
        return '/'.join(['0']*6)

def write_spagedi(dic_T,header, targets, config):
    ploidy = int(config['ploidy'])
    outdir = config['data']
    fpopinfo= os.path.join(config['data'],'subpop.info')
    num_sample = len(header)
    num_targets= len(targets)
    if os.path.isfile(fpopinfo):
        popinfo, subpop = read_popinfo(fpopinfo)
        num_pop = len(subpop)
        with open('%s/spagedi.txt' %outdir,'w') as out:
            out.write('%s\t%s\t2\t%s\t1\t%s\n0\n' %(num_sample, num_pop, num_targets,ploidy))
            out.write('Ind\tCat\tX\tY\t') #write first two lines
            out.write('\t'.join(targets)+'\n')
            for sample in dic_T:
                out.write('%s\t%s\t1\t1\t' %(sample,popinfo[sample]))
                for locus in targets:
                    genotype = spagedi_genotype(dic_T[sample][locus])
                    out.write(genotype+'\t')
                out.write('\n')
            out.write('END')
        logging.info('Data for spagedi prepared.')
    else:
        logging.info('No pop.info file, spagedi.txt woule not be created.')

def write_for_poppr(dic_T, header, targets, config):
    ploidy = int(config['ploidy'])
    outdir = config['data']
    fpopinfo= os.path.join(config['data'],'subpop.info')
    # dic_T = dict_transpose(dic_motif)
    num_sample = len(header)
    num_targets= len(targets)
    with open('%s/ssr_poppr.csv' %outdir,'w') as out:
        out.write('%s,%s,1\n\nInd,Pop,' %(num_targets, num_sample))
        for target in targets:
            out_target = target_prepare(target, ploidy)
            out.write('%s,' % out_target)
        out.write('\n')
        for i,sample in enumerate(dic_T):
            if os.path.isfile(fpopinfo):
                popinfo, subpop = read_popinfo(fpopinfo)
                out.write('%s,%s,' % (sample,popinfo[sample]))
            else:
                out.write('%s,%s,' % (sample,1))
            # out.write('%s,%s,' % (sample,i+1))
            for ssr in targets:
                # motifs = dic_motif[ssr][sample]
                motifs = dic_T[sample][ssr]
                for i in range(int(ploidy)):
                    try:
                        r = motifs[i]
                        # r = motifs[i] if motifs[i] != -9 else 0
                    except IndexError:
                        # print 'Index out of range.'
                        r = 0
                    out.write('%s,' % r)
            out.write('\n')
    logging.info('Data for poppr prepared.')
    return None

def target_prepare(target, ploidy):
    """The formatted first line of poppr file."""
    tmp = [' ' for i in range(int(ploidy) -1)]
    tmp.insert(0,target)
    return ','.join(tmp)

def write_for_structure(dic_T,targets,config):
    ploidy = int(config['ploidy'])
    outdir = config['data']
    fpopinfo= os.path.join(config['data'],'subpop.info')
    # dic_T = dict_transpose(dic_motif)
    sample_count = 0
    if os.path.isfile(fpopinfo):
        popinfo, subpop = read_popinfo(fpopinfo)
        dic_pop = {j:(i+1) for i,j in enumerate(subpop)}
        with open(os.path.join(outdir,'infile.names'), 'w') as oh:
            for (sub, num) in dic_pop.items():
                oh.write('%s\t%s\n' %(num, sub))
    with open(os.path.join(outdir,'ssr_structure'), 'w') as out:
        out.write(' \t \t' + '\t'.join(targets)+'\n')
        for sample in dic_T:
            sample_count += 1
            for i in range(ploidy):
                if os.path.isfile(fpopinfo):
                    out.write(sample + '\t%s' %dic_pop[popinfo[sample]])
                else:
                    out.write(sample+'\t1')
                for target in targets:
                    motifs = dic_T[sample][target]
                    try:
                        r = motifs[i]
                    except IndexError:
                        # print 'Index out of range.'
                        r = -9
                    out.write('\t'+str(r))
                out.write('\n')
    logging.info('Structure pre_file put down.')
    return None

def data_prepare(config):
    ssr_xls = config['raw_ssr']
    workbook = load_workbook(ssr_xls)
    dic_depth, targets, samples = read_depth(workbook)
    dic_type, target_str = read_typing(workbook)
    # mydic,header,targets = load_data(fssr,func=cal_reps)
    # dic_T = dict_transpose(mydic)
    write_for_structure(dic_type,targets,config)
    write_for_poppr(dic_type,samples,targets,config)
    write_spagedi(dic_type,samples,targets,config)

if __name__ == '__main__':
    configfile = sys.argv[1]
    config = utils.read_config(configfile)
    # mydic,header,targets = load_data(fssr,func=cal_reps)
    # dic_num = to_binary_data(mydic)
    # write_out(dic_num,targets)
    # write_for_structure(mydic,targets,ploidy)
    # write_for_poppr(mydic,header,targets,ploidy)
    data_prepare(config)

