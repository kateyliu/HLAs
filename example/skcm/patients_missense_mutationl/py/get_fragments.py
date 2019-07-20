import numpy as np
# import subprocess
import commands
import sys
import os
# sys.path.insert(1,'/usr/lib64/python2.6/site-packages')
import numpy as np
from numpy import genfromtxt

# from pathlib import Path
# commands.getstatusoutput('samtools faidx /scratch/yang/TCGA/Homo_sapiens.GRCh38.pep.rename.fa ENSP00000308051')

# file  = sys.argv[1]
# name = sys.argv[2]

# file = 'ff0f875f-113a-4c0f-8c0a-4259cd496e6d.txt'
# name = 'ff0f875f-113a-4c0f-8c0a-4259cd496e6d'
folder_name = 'patients_missense_mutation'

filelist = os.listdir(folder_name)
for f_index in range(len(filelist)):
    files = filelist[f_index]
    if 'txt' in files:
        name = files.split('.')[0]
        print(name)
        file = files
        print("processing " + file + "\n")
        with open(folder_name + '/' + file) as f:
            lines = f.readlines()
        if len(lines) > 0:
            outfile_wt = open(folder_name + '/wt/' + name + '_wt.fa', 'w')
            outfile_mt = open(folder_name + '/mt/' + name + '_mt.fa', 'w')
            outfile = open(folder_name + '/all/' + name + '_all.txt', 'w')
            for x in lines:
                arr = x.split("\t")
                snpid = arr[0] + '_' + arr[1] + '_' + arr[2] + '_' + arr[7] + '_' + arr[9] + "_" + arr[11] + "_" + arr[
                    14].rstrip("\n")
                # print("snpid:",snpid)
                # ('snpid:', 'ZFHX4_chr8_76853798_G_T_D/Y_9536e32d-2707-48d2-a36d-08c521665bb9')
                # snpid = arr[1]+'_'+arr[2]+'_'+arr[11].split("/")[0]+'_'+arr[11].split("/")[1]+"_"+arr[14].rstrip("\n")
                patient_id = arr[14]
                ensp = arr[12]
                # print('arr[10]\n',arr[10])
                try:
                    pos, length = int(arr[10].split("/")[0]), int(arr[10].split("/")[1])
                except ValueError:
                    print('invalid literal for int() with base 10:', arr[10])
                    print('in this file:', name)
                    with open('invaild.pos', 'w') as p:
                        p.write(name + '.txt' + '\n')
                        p.write(arr[10])
                    pos, length = str(arr[10].split("/")[0]), str(arr[10].split("/")[1])
                if type(pos) is int and type(length) is int:
                    if (pos < 9):  # pos start from 1 ensemble position start from 1
                        startpos = 1
                        idx = pos - 1  # center index start from 0
                    else:
                        startpos = pos - 8
                        idx = 8  # center
                    range = length - pos
                    if (range < 8):
                        endpos = length
                    else:
                        endpos = pos + 8
                    if len(arr[11].split("/")) == 2:
                        wt, mt = arr[11].split("/")[0], arr[11].split("/")[1]
                        comm = 'samtools faidx Homo_sapiens.GRCh38.pep.rename.fa ' + ensp + ':' + str(startpos) + '-' + str(
                            endpos)  # from startpos from 1 to endpos include endpos
                        # print (comm)
                        # samtools faidx Homo_sapiens.GRCh38.pep.rename.fa ENSP00000289877:239-255
                        # comm=['samtools','faidx','Homo_sapiens.GRCh38.pep.rename.f',ensp +':'+str(startpos)+'-'+str(endpos)]
                        # output = subprocess.Popen(comm)
                        output = commands.getstatusoutput(comm)
                        if (output[0] == 0 and ('\n' in output[1])):
                            seq_id = '>' + snpid
                            # print (seq_id)
                            wtseq = ''.join(output[1].split('\n')[1:])
                            outfile_wt.write(seq_id + "\n")
                            outfile_wt.write(wtseq + "\n")
                            # print('wt: '+wtseq)
                            mtseq = wtseq[0:idx] + mt + wtseq[idx + 1:]
                            outfile_mt.write(seq_id + "\n")
                            outfile_mt.write(mtseq + "\n")
                            content = snpid + '\t' + arr[8] + '\t' + wt + '\t' + mt + '\t' + wtseq + '\t' + mtseq + '\n'
                            outfile.write(content)
                # print('mt: '+mtseq)
            outfile_wt.close()
            outfile_mt.close()
            outfile.close()

# module add python/python-2.7.14
# module add samtools/samtools-1.7
# module add bcftools/bcftools-1.7
# module add htslib/htslib-1.7


'''
file  = '/scratch/yang/TCGA/patients_id'
file  = 'tmp_my_id'

with open(file) as f:
    lines = f.readlines()


i = 1
for x in lines:

	print i
	i =i+1
	#id = x.rstrip("\n")
	id = x.split('_')[4].rstrip("\n")
	comm = "grep -A 1 \'"+x.rstrip("\n")+"\' all_wt.fa >> patients/wt/"+id+".fa"
	#comm = '~/netMHCpan-4.0/netMHCpan -l 9 -BA -xls -xlsfile patients/wt/'+id+'.mt.MHC.xls patients/wt/'+id + '.fa'
	#print comm
	tmp = commands.getstatusoutput(comm)
	del tmp


def split_seq(file):
	with open(file) as f:
		lines = f.readlines()
	for x in lines:
		id = x.split('_')[4].rstrip("\n")
		comm = "grep -A 1 \'"+x.rstrip("\n")+"\' all_wt.fa >> patients/wt/"+id+".fa"
		tmp = commands.getstatusoutput(comm)
		del tmp

def MHC_pred(file):
	with open(file) as f:
		lines = f.readlines()
	i = 1
	for x in lines:
		print i 
		i = i+1
		id = x.rstrip("\n")
		comm = '~/netMHCpan-4.0/netMHCpan -l 9 -BA -xls -xlsfile patients/wt/'+id+'.wt.MHC.xls patients/wt/'+id + '.fa'
		tmp = commands.getstatusoutput(comm)
		del tmp


def amplitude(wt,mt):
		return (wt/mt)*(1/(1+0.0003*wt))

for x in lines:
	id = x.rstrip("\n")
	print id
	wt_file = '/scratch/yang/TCGA/patients/wt/'+id+'.wt.MHC.xls'
	mt_file = '/scratch/yang/TCGA/patients/mt/'+id+'.mt.MHC.xls'
	out_file = open('/scratch/yang/TCGA/patients/amp_non/'+id+'amp.non.txt','a')
	with open(wt_file) as f1:		
		lines_wt = f1.readlines()[2:]
	with open(mt_file) as f2:		
		lines_mt = f2.readlines()[2:]
	for m,w in zip(lines_wt,lines_mt):

		mt,wt = m.split('\t'),w.split('\t')
		if(float(wt[6]) <= 500 and float(mt[6]) > 500 ):
			#print mt[6],wt[6]
			mt_frag,wt_frag,id = mt[3],wt[3],mt[2]
			mt_nM,wt_nM = mt[6],wt[6]
			amp = str(amplitude(float(wt[6]),float(mt[6])))
			content = id + '\t'+mt_frag+'\t'+wt_frag+'\t'+mt_nM+'\t'+wt_nM+'\t'+amp+'\n'
			out_file.write(content)
	f1.close()
	f2.close()
	out_file.close()



def get_amp(file):
	data = genfromtxt(file,delimiter = '\t',dtype = '|U50')
	return max(data[:,5].astype(float))

outfile = open('/scratch/yang/TCGA/patients/maxamp.out.txt','a')
i =1
for x in lines:
	id = x.rstrip("\n")
	print i,id
	i =i+1
	infile = '/scratch/yang/TCGA/patients/amp/' + id +'amp.txt'
	if(os.stat(infile).st_size != 0):
		data = genfromtxt(infile,delimiter = '\t',dtype = '|U50',skip_header = 0)
		if (data.ndim != 1):
			max_amp = max(data[:,5].astype(float))
		else:
			max_amp = data[5].astype(float)
	else:
		max_amp = 0

	content =  id + '\t'+str(max_amp)+'\n'
	outfile.write(content)

outfile.close()
'''
