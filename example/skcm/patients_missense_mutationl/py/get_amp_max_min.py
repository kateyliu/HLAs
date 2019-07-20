import sys
import os
#sys.path.insert(1,'/usr/lib64/python2.6/site-packages')
import numpy as np
from numpy import genfromtxt

#os.chdir('/home/wangdu/data/MHC_yang/BREAST_CANCER')
file  = '/home/dwang/data/MHC_yang/KIDNEY_KIRC/mutec/patients_id.txt' #1084 patiences
HLA_list = '/home/dwang/data/MHC_yang/KIDNEY_KIRC/HLA_list'


with open(file) as f:
    lines = f.readlines() #patiences

HLAs=[]
with open(HLA_list) as f:
    for temp in f:
        HLAs.append(temp.strip())




#outfile = open('/home/dwang/data/MHC_yang/BREAST_CANCER/patients_missense_mutation/maxamp_matrix.out.txt','w')
outfile = open('./allamp_matrix.out.txt','w')


outfile.write("id\t"+"\t".join(HLAs)+"\n")

for x in lines:
    id = x.rstrip("\n") #patient id
    print id
    amp_list=[]
    for hla in HLAs:
      if os.path.exists('/home/dwang/data/MHC_yang/KIDNEY_KIRC/mutec/patients_missense_mutation/amp/'+hla+"/"+id+'.amp.txt'):
            infile = '/home/dwang/data/MHC_yang/KIDNEY_KIRC/mutec/patients_missense_mutation/amp/'+hla+"/" + id +'.amp.txt'
            if(os.stat(infile).st_size != 0):
                data = genfromtxt(infile,delimiter = '\t',dtype = '|U50',skip_header = 0)
                if (data.ndim != 1):
                    max_amp = max(1/(data[:,5].astype(float)))
                else:
                    max_amp = data[5].astype(float)
            else:
                max_amp = -1
            amp_list.append(str(max_amp))
    
    if amp_list!=[]:
       content =  id + '\t'+"\t".join(amp_list)+'\n'
       outfile.write(content)

outfile.close()

