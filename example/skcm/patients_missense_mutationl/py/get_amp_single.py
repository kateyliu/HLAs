import numpy as np
# import subprocess
import commands
import sys
import os
# sys.path.insert(1,'/usr/lib64/python2.6/site-packages')
import numpy as np
from numpy import genfromtxt
# commands.getstatusoutput('samtools faidx /scratch/yang/TCGA/Homo_sapiens.GRCh38.pep.rename.fa ENSP00000308051')
# os.chdir('/home/dwang/data/MHC_yang/KIDNEY_KIRC/mutec/')
HLA = sys.argv[1]
#HLA = 'HLA-A01:20'# just for test
file = 'patients_id.txt'
# HLA_list = '/home/ylmk2/data/TCGA/HLA-B15:02-09'
folder_name = 'patients_missense_mutation/'

def amplitude(wt, mt):
    return (wt / mt) * (1 / (1 + 0.0003 * wt))

def MHC_pred(id, HLA, folder_name):
    if not os.path.exists(folder_name + 'wt/' + HLA):
        os.makedirs(folder_name + 'wt/' + HLA)
    if not os.path.exists(folder_name + 'mt/' + HLA):
        os.makedirs(folder_name + 'mt/' + HLA)
    if not os.path.exists(folder_name + 'amp/' + HLA):
        os.makedirs(folder_name + 'amp/' + HLA)
    if os.path.exists(folder_name + 'wt/' + id + '_wt.fa'):  # todo check comm
        comm = '/home/tangbo/apps/netMHCpan-4.0/netMHCpan -a ' + HLA + ' -l 9 -BA -xls -xlsfile ' + folder_name + 'wt/' + HLA + '/' + id + '.MHC.xls ' + folder_name + 'wt/' + id + '_wt.fa'
        print(comm)
        tmp = commands.getstatusoutput(comm)
        # os.system(comm)
        comm = '/home/tangbo/apps/netMHCpan-4.0/netMHCpan -a ' + HLA + ' -l 9 -BA -xls -xlsfile ' + folder_name + 'mt/' + HLA + '/' + id + '.MHC.xls ' + folder_name + 'mt/' + id + '_mt.fa'
        tmp = commands.getstatusoutput(comm)
        # os.system(comm)
        if tmp[0] == 0:  # tem re aasinged the results should be tmp1 tmp2 above
            wt_file = folder_name + 'wt/' + HLA + '/' + id + '.MHC.xls'
            mt_file = folder_name + 'mt/' + HLA + '/' + id + '.MHC.xls'
            out_file = open(folder_name + 'amp/' + HLA + '/' + id + '.amp.txt', 'w')
            with open(wt_file) as f1:
                lines_wt = f1.readlines()[2:]  # start 2 lines
            with open(mt_file) as f2:
                lines_mt = f2.readlines()[2:]
            for m, w in zip(lines_wt, lines_mt):
                mt, wt = m.split('\t'), w.split('\t')
                try:
                    if (float(wt[6]) <= 500 and float(mt[6]) > 500):  # todo ???meaning
                        mt_frag, wt_frag, id = mt[3], wt[3], mt[2]
                        mt_nM, wt_nM = mt[6], wt[6]
                        amp = str(amplitude(float(wt[6]), float(mt[6])))
                        content = id + '\t' + mt_frag + '\t' + wt_frag + '\t' + mt_nM + '\t' + wt_nM + '\t' + amp + '\n'
                        # print content
                        out_file.write(content)
                except IndexError:
                    with open("failed.txt", 'a+') as fail:
                        fail.write(mt_file + "\n")
                        fail.write(mt + "\n")
                        fail.write(wt_file + "\n")
                        fail.write(wt + "\n")
                    print("IndexError: list index out of range in", wt_file, "or", mt_file)
            f1.close()
            f2.close()
            out_file.close()
    else:
        print('*.fa   do not exits, make sure step4 finished correctly')       # del tmp


with open(file) as f:
    lines = f.readlines()


# with open(HLA_list) as f:
# HLAs = f.readlines()
HLA = HLA.rstrip("\n")
for i in range(0, len(lines)):#if has header use the 1
    y = lines[i]
    id = y.rstrip("\n")  # remove right \n
    MHC_pred(id, HLA, folder_name)



"""
IOError: [Errno 2] No such file or directory: 'patients_missense_mutation/wt/HLA-A01:01/26a92e1b-08d$
"""

