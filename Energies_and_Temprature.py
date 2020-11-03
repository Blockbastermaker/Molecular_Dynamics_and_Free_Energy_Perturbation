#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 21:08:19 2020

@author: nour
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools


import glob
lambdas=[]
for filename in glob.glob("/Users/nour/Desktop/Q/17_3steps_2/FEP*.log"):
        lambdas.append(filename)
lambdas.reverse()  
### ele and vdw Potentails     
ele=[]
vdw=[]
for x in lambdas:       
    fo = open(x, "r")
    e=[]
    v=[]
    for i in fo:
        if "Q-any" in i and not "*" in i:
            e.append(float(i.split()[3]))
            v.append(float(i.split()[4]))
    ele.append(e)
    vdw.append(v)
    fo.close()
    
elea=list(itertools.chain(*ele))
vdwa=list(itertools.chain(*vdw))

plt.plot(elea,label=" Electrostatics Energy")
plt.plot(vdwa,label="vdW Energy")
plt.legend(loc="best")
plt.savefig('ele_vdw.png',dpi=200)
plt.close()
#### Total Energies
Ut=[]
Up=[]
Uk=[]
for x in lambdas:       
    fo = open(x, "r")
    t=[]
    p=[]
    k=[]
    for i in fo:
        if "SUM" in i and not "*" in i and not "Q-SUM" in i and not "0.000" in i :
           t.append(float(i.split()[1]))
           p.append(float(i.split()[2]))
           k.append(float(i.split()[3]))
    Ut.append(t)
    Up.append(p)
    Uk.append(k)
    fo.close()
Uta=list(itertools.chain(*Ut))
Upa=list(itertools.chain(*Up))
Uka=list(itertools.chain(*Uk))
plt.plot(Uta,label="Total Energy")
plt.plot(Upa,label="Potentail Energy")
plt.plot(Uka,label="Kinetic Energy")
plt.legend(loc="best")
plt.savefig('Energies.png',dpi=200)
plt.close()

#for i in Ut,Up,Uk:
#    for x in i:
#        plt.plot(x)
  
#plt.plot(list(range(0,len(Up[0]))),Up[0])
#plt.plot(list(range(0,len(Uk[0]))),Uk[0])


### Temperature

Tm=[]
for x in lambdas:       
    fo = open(x, "r")
    m=[]
    for i in fo:
        if "Temperature at step" in i and not "*" in i:
            m.append(float(i.split("=")[2]))
    Tm.append(m)
fo.close()

Tm=list(itertools.chain(*Tm))
plt.plot(Tm)
plt.savefig('Temperature.png',dpi=200)
plt.close()
#
#for i in Tm:
#    if len(i) != 0:
#        print(max(i))


