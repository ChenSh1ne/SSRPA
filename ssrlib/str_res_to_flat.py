## !/usr/bin/python
# -*- coding: utf-8 -*-

import xlrd
import re
import sys
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
    return sample_depth

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
    for i, j in enumerate(rowvals[3:]):
        if re.match(r'\d+',str(j)):
            motif = header[3:][i]
            # target = int(motif)
            tmp = [int(motif)] * int(j)
            alleles.extend(tmp)
            # alleles.append('%s*%s' %(int(j),int(motif)))
    return alleles

def dict_transpose(dic_motif):
    new_dic = defaultdict(lambda: defaultdict())
    for ssr in dic_motif:
        for sample in dic_motif[ssr]:
            new_dic[sample][ssr] = dic_motif[ssr][sample]
    return new_dic

def write_genotype_txt(dic_type):
    samples = dic_type.keys()
    print samples
    new_dic = dict_transpose(dic_type)
    targets = new_dic.keys()
    print targets
    with open('ssr_raw.txt','w') as out:
        out.write('Targets'+'\t'+'\t'.join(samples)+'\n')
        for t in targets:
            out.write(t+'\t')
            for s in samples:
                out.write(new_dic[t][s]+'\t')
            out.write('\n')

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
            if alleles and depth >= 50:
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


if __name__ == '__main__':
    f = sys.argv[1]
    sample_depth = read_depth(load_workbook(f))
    dic_type, target_str = read_typing(load_workbook(f))
    write_alleles_xlsx(dic_type, target_str, sample_depth)


