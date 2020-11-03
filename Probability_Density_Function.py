#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 18:14:22 2019

@author: nour
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.ioff()
import os
import glob

os.chdir("//Users//nour//Desktop//Q//FEP_analysis//propranolol//active")

lambdas=[]
for filename in glob.glob("FEP2*.log"):
        lambdas.append(filename)
lambdas.reverse()  
      
U=[]
for x in lambdas:       
    fo = open(x, "r")
    y=[]
    for i in fo:
        if "Q-SUM" in i and not "*" in i:
            y.append(float(i.split()[3]))
    U.append(y)
    fo.close()




dU=[]
dU_s=[]
for i in U:       
    o=[]
    p=[]
    for x in range (0,len(i),2):
        o.append(i[x+1] -i[x])
        p.append(((i[x+1] -i[x])**2))
    dU.append(o)
    dU_s.append(p)



mean_dU=[]
for i in dU:
   mean_dU.append((np.mean(i)))   
   
mean_dU_s=[]
for i in dU_s:
   mean_dU_s.append((np.mean(i)))  



sigma_s=[]
for i in range(len(mean_dU)):
   sigma_s.append(mean_dU_s[i] -  (mean_dU[i]**2)) 



sigma=[]
for i in sigma_s:
   sigma.append(math.sqrt(i)) 


dt=[]

for i in dU[1:11]:
    dt.append([i])
    
for i in dU[-11:-2]:
    dt.append([i])

for i in dt:
   p=sns.distplot(i, hist = False, kde = True, color="gray",
                 kde_kws = {'lw': 1} )
p=sns.distplot(dU[0], hist = False, kde = True,color="blue",
                 kde_kws = {'shade': True, 'lw': 2},label="State A")

p=sns.distplot(dU[1], hist = False, kde = True, color="gray",
                 kde_kws = {'lw': 1}, label="Intermediate States" )

p=sns.distplot(dU[-1], hist = False, kde = True,color="red",
                 kde_kws = {'shade': True, 'lw': 2},label="State B")



plt.title(' Probability Density Function of ΔU',fontsize=16)
plt.xlabel("ΔU")
plt.ylabel("Probability")
p.figure.savefig('PDF.png',dpi=1000)


#
#