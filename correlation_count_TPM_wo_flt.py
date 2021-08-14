#!/usr/bin/python
# -*- coding: UTF-8 -*-
import numpy as np
from scipy import stats
import math
import sys, getopt

def main(argv):
    countfile = ''
    sinputfile = ''
    dinputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hc:s:d:",["cfile=","sfile=","dfile="])
    except getopt.GetoptError:
         print('correlation_count_TPM_flt.py -c <countfile> -s <inputsamedomainfile> -d <inputdifdomainfile>')
         sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
           print('correlation_count_TPM_flt.py -c <countfile> -s <inputsamedomainfile> -d <inputdifdomainfile>')
           sys.exit()
        elif opt in ("-c"):
           countfile = arg
        elif opt in ("-s"):
           sinputfile = arg
        elif opt in ("-d"):
           dinputfile = arg
    htseq=open(countfile)
    matrix={}
    for line in htseq.readlines():
        line=line.strip('\n')
        info=line.split('\t')
        name=info[0]
        del info[0]
        for index, item in enumerate(info):
            info[index] = float(item)
        matrix[name]=info[1:]
    htseq.close
    del info[:]

    same=open(sinputfile)
    if countfile.find('wt') != -1:
        spre=sinputfile.replace(".txt","_wt_TPM_wo_flt.txt")
    elif countfile.find('mut') != -1:
        spre=sinputfile.replace(".txt","_mut_TPM_wo_flt.txt")
    sresfile=spre.replace("domain","result")
   # sfltfile=spre.replace("domain","flt")
    sres=open(sresfile,'w')
    sres.write('geneid1\tgeneid2\tgene1\tgene2\tdistance\tcorrelation\tpvalue\tTAD_ID\n')
    #sflt=open(sfltfile,'w')
    for line in same.readlines():
        line=line.strip('\n')
        info=line.split('\t')
        if info[0] in matrix and info[1] in matrix:
            x=np.array(matrix[info[0]])
            y=np.array(matrix[info[1]])
        #    if sum(x<1)<=(len(x)/2) and sum(y<1)<=(len(y)/2):
            (cof,pv)=stats.pearsonr(x,y)
            cof=float(cof)
            pv=float(pv)
            if math.isnan(cof) == False and cof:
                cob=info[0],info[1],info[3],info[4],info[2],str(cof),str(pv),info[5]
                sres.write('\t'.join(cob)+'\n')
        #   elif sum(x<1)>(len(x)/2):
            #    if sum(y<1)>(len(y)/2):
              #      sflt.write(line+'\tboth\n')
             #   else:
              #      sflt.write(line+'\tgene1\n')
           # elif sum(y<1)>(len(y)/2):
            #    sflt.write(line+'\tgene2\n')
    same.close
    sres.close
    #sflt.close
    del info[:]
    if countfile.find('wt') != -1:
        dpre=dinputfile.replace(".txt","_wt_TPM_wo_flt.txt")
    elif countfile.find('mut') != -1:
        dpre=dinputfile.replace(".txt","_mut_TPM_wo_flt.txt")
    dresfile=dpre.replace("domain","result")
#dfltfile=dpre.replace("domain","flt")
    dif=open(dinputfile)
    dres=open(dresfile,'w')
    dres.write('geneid1\tgeneid2\tgene1\tgene2\tdistance\tcorrelation\tpvalue\tTAD_ID1\tTAD_ID2\n')
    #dflt=open(dfltfile,'w')
    for line in dif.readlines():
        line=line.strip('\n')
        info=line.split('\t')
        if info[0] in matrix and info[1] in matrix:
            x=np.array(matrix[info[0]])
            y=np.array(matrix[info[1]])
          #  if sum(x<1)<=(len(x)/2) and sum(y<1)<=(len(y)/2):
            (cof,pv)=stats.pearsonr(x,y)
            cof=float(cof)
            pv=float(pv)
            if math.isnan(cof) == False and cof:
                cob=info[0],info[1],info[3],info[4],info[2],str(cof),str(pv),info[5],info[6]
                dres.write('\t'.join(cob)+'\n')
           # elif sum(x<1)>(len(x)/2):
         #       if sum(y<1)>(len(y)/2):
          #          dflt.write(line+'\tboth\n')
               # else:
                #    dflt.write(line+'\tgene1\n')
       #     elif sum(y<1)>(len(y)/2):
        #        dflt.write(line+'\tgene2\n')
    dif.close
    dres.close
    #dflt.close
if __name__ == "__main__":
   main(sys.argv[1:])
