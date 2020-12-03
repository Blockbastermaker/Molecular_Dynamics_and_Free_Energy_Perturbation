#%%
import struct
import numpy as np
import itertools
import os
import glob
from numpy.core.defchararray import array
from numpy.core.overrides import verify_matching_signatures
from numpy.lib.function_base import average
import pandas as pd
import re

from pandas.core.frame import DataFrame

#%%


def ReadBinary(EnergyFiles_Lst):

    State_A_RawEnergies = []
    State_B_RawEnergies = []

    for file in EnergyFiles_Lst:
        with open(file,'rb') as f: fileContent = f.read()

        Version_Check_Str=str(list(struct.unpack("c" * ((len(fileContent[32:112]))//1),fileContent[32:112]))).replace("b'", "").strip("[],' '").replace("'","").replace(",","").replace(" ","")
  
        EnergyFileLength_int=len(fileContent)
        BinaryChankSize_int=120
        NextBinaryChank_int=272
        HeaderSize_int=124
        State_B_Shift_int=132

        ## Check Q_Energies Version !!!

        if "6." in Version_Check_Str: HeaderSize_int += 4
        
        if not '5.' in Version_Check_Str and not "6." in Version_Check_Str:
            print("Pleaes Check the your Qdyn verion in file: "+file +" ----> format is NOT Supported !!! ")
            exit()
            
  

        for Byte in range(HeaderSize_int, EnergyFileLength_int, NextBinaryChank_int):
            State_A_Lst=struct.unpack("d" * ((len(fileContent[Byte:Byte+BinaryChankSize_int]))//8),fileContent[Byte:Byte+BinaryChankSize_int])
            State_B_Lst=struct.unpack("d" * ((len(fileContent[Byte+State_B_Shift_int:Byte+State_B_Shift_int+BinaryChankSize_int]))//8),fileContent[Byte+State_B_Shift_int:Byte+State_B_Shift_int+BinaryChankSize_int])
            State_A_RawEnergies.append(State_A_Lst)
            State_B_RawEnergies.append(State_B_Lst)
        
    return State_A_RawEnergies, State_B_RawEnergies


def createDataFrames(rawEnergy):

    Columns_name=["Lambda","Q_sum","Q_bond","Q_angle","Q_torsion","Q_improper","Q_any_ele","Q_any_vdw","Q_Q_ele","Q_Q_vdw","Q_protein_ele","Q_protein_vdw","Q_water_ele","Q_water_vdw","Restraints"]

    return pd.DataFrame(rawEnergy, columns=Columns_name)

#%%

#if '__name__' == '__main__':

#%%
    #os.chdir("/Users/nour/New_qfep") #MAC
    os.chdir("Z:/jobs/Qfep_NEW/qfep_small")
    EnergyFiles_Lst = [filename for filename in glob.glob("*.en")]  
    State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = ReadBinary(EnergyFiles_Lst)
    State_A_df = createDataFrames(State_A_RawEnergies_Lst)
    State_B_df = createDataFrames(State_B_RawEnergies_Lst)





##Zwnazig
# %%

#Zwanzig_exp= -0.592*(np.log(np.mean(np.exp((State_A_df["Lambda"]*State_A_df["Q_sum"] + State_B_df["Lambda"]*State_B_df["Q_sum"])*0.592))))
#Zwanzig_exp



# %%
# Energies_df=pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })

# %%
# -0.592*(np.log(np.mean(np.exp((-27.05200000000002)*-0.592))))

# %%

# E0= State_A_df["Q_sum"]*State_A_df["Lambda"]+State_B_df["Q_sum"]*State_B_df["Lambda"]
# E1= State_A_df["Q_sum"]*State_B_df["Lambda"]+State_B_df["Q_sum"]*State_A_df["Lambda"]


# dE=E1-E0
# dE= np.exp(-dE/0.592)

# lambdas_dE=[]
# for i in range(0,len(dE),int(len(dE)/2)):
#     lambdas_dE.append(i)

# zz = [sum(dE[i:i+int(len(dE)/2)])/len(dE[i:i+int(len(dE)/2)]) for i in lambdas_dE]
# # for i in lambdas_dE:

# #     z.(sum(dE[i:i+int(len(dE)/2)])/len(dE[i:i+int(len(dE)/2)]))
    
# dg=-0.592*np.log(zz)
# dg
# %%


# %%
# Energies_df=pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })

# %%

# A=pd.DataFrame(dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list)))
# B=pd.DataFrame(dict(Energies_df.groupby('State_B_Lambda',sort=False)['State_B_G'].apply(list)))

# E0=(A*A.columns).values+(B*B.columns).values
# E1=(A*list((A.columns).values)[::-1]).values+(B*list((B.columns).values)[::-1]).values
# dE=pd.DataFrame(E1-E0)

# dG=pd.DataFrame(-0.592*np.log(np.mean(np.exp(-pd.DataFrame(E1-E0)/0.592))))


# %%  
Energies_df=pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })
Energies_df=Energies_df.sort_values('State_A_Lambda')
A_dict=dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list))
B_dict=dict(Energies_df.groupby('State_B_Lambda',sort=False)['State_B_G'].apply(list)) 


dG=[]
dEs=pd.DataFrame()
for i in range(len(A_dict.keys())-1): #reversed()
    # A=pd.DataFrame({((list(A_dict.keys())[i])):(A_dict.get(list(A_dict.keys())[i])),((list(A_dict.keys())[i+1])):(A_dict.get(list(A_dict.keys())[i+1]))})
    # B=pd.DataFrame({((list(B_dict.keys())[i])):(B_dict.get(list(B_dict.keys())[i])),((list(B_dict.keys())[i+1])):(B_dict.get(list(B_dict.keys())[i+1]))})
    A = pd.DataFrame.from_dict({((list(A_dict.keys())[i])):(A_dict.get(list(A_dict.keys())[i])),((list(A_dict.keys())[i+1])):(A_dict.get(list(A_dict.keys())[i+1]))}, orient='index').transpose()
    B = pd.DataFrame.from_dict({((list(B_dict.keys())[i])):(B_dict.get(list(B_dict.keys())[i])),((list(B_dict.keys())[i+1])):(B_dict.get(list(B_dict.keys())[i+1]))}, orient='index').transpose()
    E0=(A*A.columns).values+(B*B.columns).values
    E1=(A*list((A.columns).values)[::-1]).values+(B*list((B.columns).values)[::-1]).values
    #print(((list(A_dict.keys())[i])),((list(A_dict.keys())[i+1])),((list(B_dict.keys())[i])),((list(B_dict.keys())[i+1])))
    #dEs[str((list(A_dict.keys())[i]))+"_"+str((list(A_dict.keys())[i+1]))+"-"+str((list(B_dict.keys())[i]))+"_"+str((list(B_dict.keys())[i+1]))]=E1-E0
    dE=pd.DataFrame(E1-E0,columns=[str(A.columns.values[0])+"_"+str(B.columns.values[0])+"-"+str(A.columns.values[1])+"_"+str(B.columns.values[1]),str(A.columns.values[1])+"_"+str(B.columns.values[1])+"-"+str(A.columns.values[0])+"_"+str(B.columns.values[0])])
    dEs=dEs.append(dE.transpose())
    dG.append(-0.592*np.log(np.mean(np.exp(-dE/0.592))))
dEs=dEs.transpose()

dG=[item for sublist in dG for item in sublist]
dGF=[]
dGR=[]
#dGF=[k for k in range(len(dG)) if k %2]
for i in range(1,len(dG),2):
    dGF.append(dG[i])
    dGR.append(dG[i-1])
   #print( A_dict.get(list(A_dict.keys())[i]))
Zwanzig_df=pd.DataFrame.from_dict({"dGF":dGF,"dGR":dGR})
Zwanzig_exp=-np.mean(abs(np.sum((Zwanzig_df))))
Zwanzig_exp
# %%
#def Zwnazig(dEs_df):
dEs_df=pd.DataFrame(-0.592*np.log(np.mean(np.exp(-dEs/0.592))))
Lambdas=[]
dGF=[]
dGF_sum=[]
dGR=[]
dGR_sum=[]
dG_Average=[]
dGR.append(0.0)
dG_Average.append(0.0)

for i in range(1,len(dEs_df.index),2):
    Lambdas.append(re.split('_|-',dEs_df.index[i-1])[1])
    dGF.append(dEs_df.iloc[i,0])
    dGR.append(dEs_df.iloc[i-1,0])
Lambdas.append(re.split('_|-',dEs_df.index[-1])[1])
dGF.append(0.0)
dGF=dGF[::-1]
for i in range(len(dGF)):
    dGF_sum.append(sum(dGF[:i+1]))
    dGR_sum.append(sum(dGR[:i+1]))

dG_average_raw=(pd.DataFrame(list((Zwanzig_df["dG_Forward"][1:])))-pd.DataFrame(list(Zwanzig_df["dG_Reverse"][:-1])))/2
for i in range(len(list(dG_average_raw.values))):
    dG_Average.append(np.sum(dG_average_raw.values[:i+1]))


Zwanzig_df=pd.DataFrame.from_dict({"Lambda":Lambdas,"dG_Forward":dGF,"SUM_dG_Forward":dGF_sum,"dG_Reverse":dGR[::-1],"SUM_dG_Reverse":dGR_sum[::-1],"dG_Average":dG_Average})
Zwanzig_df

# %%

