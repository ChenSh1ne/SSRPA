#!/usr/bin/python
import os
import re
import sys
import time
import logging
import xlsxwriter
import subprocess
import utils
import rpy2.robjects as robjects
from collections import defaultdict
from rpy2.robjects.packages import importr



def load_prepared_ssr_data_to_R(csv_file,ploidy):
    """
    csv_file was generated previously by
    numerize.py from a ssr genotyping data.
    """
    poppr = importr('poppr')
    # logging ...
    logging.debug("R library 'poppr' imported.")
    ssr_genind_obj = poppr.read_genalex(csv_file, ploidy = int(ploidy), geo ='FALSE', region='FALSE',genclone='FALSE')
    func_dist = robjects.r('nei.dist')
    mydist = func_dist(ssr_genind_obj)
    logging.info('dist matrics generated.') #logging
    return mydist

def upgma_tree(dist_matrics,proj_name,outdir):
    """
    Return a UPGMA tree
    """
    logging.info('Phylogenetic analysis within UPGMA method')
    phangorn = importr('phangorn')
    logging.debug("R library 'phangorn' imported.")
    ape = importr('ape')
    logging.debug("R library 'ape' imported.")
    treeUPGMA = phangorn.upgma(dist_matrics)
    outfile = "%s/%s.newick" %(outdir,proj_name)
    ape.write_tree(treeUPGMA, file = outfile) # put tree in a [newick file] via write.tree {ape}
    return outfile

def nj_tree(dist_matrics,proj_name,outdir):
    """
    Return a Nabour-joining tree
    """
    logging.info('Phylogenetic analysis within NJ method')
    ape = importr('ape')
    treeNJ = ape.nj(dist_matrics)
    outfile = "%s/%s.newick" %(outdir,proj_name)
    ape.write_tree(treeNJ, file = outfile)
    return outfile

def tree_anno_ggtree(treefile):
    """
    Tree ornamentation with ggtree
    """
    utils.check_file(treefile)
    grdevices = importr('grDevices')
    Gtree = importr('ggtree')
    ape = importr('ape')
    raw_tree = ape.read_tree(treefile)
    grdevices.pdf(file="output.pdf")
    Gtree.ggtree(raw_tree)
    grdevices.dev_off()
    return None

def tree_visulization(raw_tree_file,outdir,proj_name,scale=100):
    """Tree visualization with ETE3.
    raw_tree_file: a newick format tree.
    """
    from ete3 import Tree, TreeStyle, NodeStyle
    from pyvirtualdisplay import Display
    display = Display(backend='xvfb') # Start Xvfb using PyVisualDisplay
    display.start()
    logging.info('Xvfb started.')
    # ps = subprocess.Popen("Xvfb -ac :%d -screen 0 1280x1024x8 > /dev/null 2>&1" % server_num, shell=True)
    # os.environ["DISPLAY"] = ":%d" % server_num
    t = Tree(raw_tree_file)
    ts= TreeStyle()
    ts.show_leaf_name = True
    # ts.show_branch_support = True
    # ts.scale = scale
    ts.mode = "c"
    t.render('%s/%s_tree.svg'%(outdir,proj_name),tree_style=ts)
    display.popen.terminate()
    return None

def diversity(ssrfile, ploidy, outpath):
    p = subprocess.call("Rscript /home/wuj/bin/SSRPop/pautils/diversity.R %s %s %s" %(ssrfile,ploidy,outpath), shell=True)
    if p != 0:
        logging.warn('diversity.R process wrong.')

def poppr_process(config):
    poppr_file = os.path.join(config['data'],'myssrdata.csv')
    # poppr_file = os.path.join(config['data'],'ssr_poppr.csv')
    ploidy = int(config['ploidy'])
    outdir = os.path.join(config['report'], 'dendrogram')
    outdir_div = os.path.join(config['report'], 'diversity')
    utils.make_dir(outdir)
    utils.make_dir(outdir_div)
    project_name = config['project']
    ssr_dist = load_prepared_ssr_data_to_R(poppr_file,ploidy)
    diversity(poppr_file, ploidy, outdir_div)
    div_writer = divParamXLS(outdir_div)
    div_writer.write_xls()
    mytree = upgma_tree(ssr_dist, project_name,outdir)
    tree_visulization(mytree,outdir,project_name)

class divParamXLS(object):
    """To save the diversity parameters into Excel."""

    def __init__(self,outdir):
        datadic = defaultdict(lambda: defaultdict()) # init a dictionary to save diversity params
        self.out = outdir
        self.dic = datadic

    def file_to_dic(self, resfile, dic):
        """Read diversity result and save information in self.dic"""
        count = 0
        for line in file(resfile):
            count += 1
            arr = line.strip().split()
            if count == 1:
                head = map(lambda x: x.strip('"'), arr)
            else:
                for i,j in enumerate(arr[1:]):
                    self.dic[arr[0].strip('"')][head[i]] = j
        return head # for writing title

    def read_div_res(self):
        divfile = os.path.join(self.out,'diversity.txt')
        sumfile = os.path.join(self.out,'summ_stats.locus')
        utils.check_file(divfile)
        utils.check_file(sumfile)
        param1 = self.file_to_dic(divfile, self.dic)
        param2 = self.file_to_dic(sumfile, self.dic)
        return param1 + param2

    def write_xls(self):
        title = self.read_div_res()
        title.insert(0,'Locus')
        outdir = self.out
        workbook = xlsxwriter.Workbook(os.path.join(outdir,'diversity.xlsx'))
        logging.info('Writing diversity parameters to %s', workbook)
        worksheet = workbook.add_worksheet('diversity')
        row = 0
        col = 0
        for l in title: #writing title
            worksheet.write(row, col, l)
            col += 1

        row += 1 #row grows to 1
        col = 0 #col refresh to 0

        for target in self.dic:
            worksheet.write(row, col, target)
            col += 1
            for d in title[1:]:
                worksheet.write(row, col, self.dic[target][d])
                col += 1
            row += 1
            col = 0




if __name__ == '__main__':
    configfile = sys.argv[1]
    config = utils.read_config(configfile)
    poppr_process(config)


