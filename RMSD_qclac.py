# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 11:57:08 2019

@author: Nour
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 20:57:49 2019

@author: Nour
"""
import math
import numpy as np
from matplotlib import pyplot as plt
import argparse
import os
def atoms_number(pdb_file,lig):
    atoms=[]
    with open(pdb_file)as pdb:
            for line in pdb:
                for x in lig:
                    if "ATOM" in line and x== (line[17:20]) in line and line[12:16].strip()[0]!= "H": 
                        atoms.append((int(line[7:11])))
    return atoms

def atomsN_number(pdb_file,n):
    atomsN=[]
    with open(pdb_file)as pdb:
            for line in pdb:
                for x in n:
                    lig = x.split(":")[0]
                    nn =x.split(":")[1]
                    if "ATOM" in line and lig== (line[17:20]) in line and line[12:16].strip()[0] != "H" and int(line[23:27]) ==int(nn): 
                        print(line)
                        atomsN.append((int(line[7:11])))
    return atomsN

def qcalc_inp(top,fep,lamb,choice,Atoms,output,dcd):
    f= open("qcalc.inp","w+")
    print(top,file = f)
    print(fep,file = f)
    for x in lamb:
         print(x, sep= " ",end=" ",file = f,)
    print("\n",choice,file = f)
    for y in Atoms:
         print(y, file = f)
    print(".",file = f)
    print(output,file = f)
    print("go",file = f)
    for z in dcd:
        print(z,file = f)
    print(".",file = f)     
    f.close()

parser = argparse.ArgumentParser(description="qcalc RMSD Calcution")
parser.add_argument("-l","--lig", nargs='+',help = "Ligand name")
parser.add_argument("-p","--pdb_file",default = 'top_p.pdb', help = "pdb file")
parser.add_argument("-t","--topology_file",default = 'dualtop.top', help = "Toplogy file")
parser.add_argument("-y","--lambda_states",default = ['1.0','0.0'],nargs='+', help = "lambda states")
parser.add_argument("-f","--FEP_file",default = 'FEP1.fep', help = "FEP file")
parser.add_argument("-o","--output_file",default = 'RMSD.txt', help = "RMSD output file")
parser.add_argument("-d","--dcd_file",nargs='+',default = ['md_0000_1000.dcd'], help = "dcd file")
group=parser.add_mutually_exclusive_group()
group.add_argument("-r","--residue", nargs='+',help = "Residue Sequence Number")
args = parser.parse_args()
if __name__ == "__main__":
    if args.residue is not None and args.lig is None:
        Atoms=atomsN_number(args.pdb_file, args.residue)
        qcalc_inp(args.topology_file,args.FEP_file,args.lambda_states,"1",Atoms,args.output_file,args.dcd_file)
    elif args.residue is None and args.lig is not None:
        Atomss=atoms_number(args.pdb_file,args.lig)
        qcalc_inp(args.topology_file,args.FEP_file,args.lambda_states,"1",Atomss,args.output_file,args.dcd_file)
    else: 
        print("please use"+ " "+ "'RMSD2.py -h'"+" ""for usege ")
qcalc = "~/QligFEP_Q6/Q6/bin/Qcalc6 < qcalc.inp && sed '/R/d' RMSD.txt >RMSD2.txt"
os.system(qcalc)

#Caculate Mean distance and plot
import math
import numpy as np
from matplotlib import pyplot as plt
plt.ioff()
fo = open('RMSD2.txt', "r")
y=[]
for i in fo:
      y.append(float(i))
fo.close()
print(np.mean(y))
print(np.std(y))
f= open("RMSD_mean.txt","w+")
print(np.mean(y),file = f)
print(np.std(y),file = f)
f.close()
plt.plot(y, c='c')
plt.title('RMSD')
plt.xlabel("Trajectory")
plt.ylabel("Distance (Å)")
plt.text(20,max(y), "Mean= "+str("{:.2f}".format(np.mean(y))))
plt.text(20,max(y)-0.2 , "STD  = "+str("{:.2f}".format(np.std(y))))
plt.savefig('RMSD.png',dpi=200)
# 
