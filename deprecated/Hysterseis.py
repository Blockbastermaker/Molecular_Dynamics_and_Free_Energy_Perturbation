#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 19:51:16 2020

@author: nour
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import pandas as pd

import glob
qfep=[]
for filename in glob.glob("qfep*.out"):
        qfep.append(filename)
        
        
dG=[]
for x in qfep:       
    fo = open(x, "r")
    g=[]
    part1 = False 
    for line in fo:
        if line.startswith("# Part 1"):
            part1 = True 
        elif line.startswith("# Min energy-gap"):
            part1 = False 
    
        if part1:
            g.append(line.strip().rstrip().split())
            
            
    del(g[0])
    del(g[-1]) 
    del(g[0][0])
    dG.append(g)
    fo.close()
dfs=[]
for i in dG:
     dfs.append("FEP_" +str(dG.index(i)+1))
     globals()["FEP_" +str(dG.index(i)+1)] = pd.DataFrame(i[1:], columns =i[0], dtype = float) 
     
for df in dfs:     
     eval(df).iloc[:,3]=eval(df).iloc[:,3].values[::-1]
     
     eval(df).iloc[:,4]=eval(df).iloc[:,4].values[::-1]

for df in dfs:
    f = df+ "-ΔGf"
    fc= "C" + str(dfs.index(df)+5)
    rc= "C" + str(dfs.index(df))
    r = df+ "-ΔGr"
    p=plt.plot(eval(df).iloc[:,2],'.',label= f,color=fc)
    p=plt.plot(eval(df).iloc[:,4],'.',label =r,color=rc)
    plt.title('Hysteresis between ΔGf and ΔGr',fontsize=16)
    plt.xlabel("λ",fontsize=14)
    plt.ylabel("ΔG FEP (Kcal/mol)",fontsize=14)
    plt.legend()
plt.savefig('Hysteresis.png',dpi=200)

