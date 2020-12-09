#%%
import struct
import numpy as np
import itertools
import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import re

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
def dE_Calculation2(steps):
    
    dEs=pd.DataFrame()
    Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('State_A_Lambda')
    Energies_df.iloc[:1]
    State_A_Energies_dict=dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list))
    State_B_Energies_dict=dict(Energies_df.groupby('State_B_Lambda',sort=False)['State_B_G'].apply(list)) 

    
    for i in range(len(State_A_Energies_dict.keys())-1):
        State_A_Energies = (pd.DataFrame.from_dict({((list(State_A_Energies_dict.keys())[i])):(State_A_Energies_dict.get(list(State_A_Energies_dict.keys())[i])),((list(State_A_Energies_dict.keys())[i+1])):(State_A_Energies_dict.get(list(State_A_Energies_dict.keys())[i+1]))}, orient='index').transpose()).iloc[:steps]
        State_B_Energies = (pd.DataFrame.from_dict({((list(State_B_Energies_dict.keys())[i])):(State_B_Energies_dict.get(list(State_B_Energies_dict.keys())[i])),((list(State_B_Energies_dict.keys())[i+1])):(State_B_Energies_dict.get(list(State_B_Energies_dict.keys())[i+1]))}, orient='index').transpose()).iloc[:steps]
        
        State_A_Lambda_float=State_A_Energies.columns
        State_B_Lambda_float=State_B_Energies.columns
        
        E0=(State_A_Energies*State_A_Lambda_float).values+(State_B_Energies*State_B_Lambda_float).values
        E1=(State_A_Energies*list((State_A_Lambda_float).values)[::-1]).values+(State_B_Energies*list((State_B_Lambda_float).values)[::-1]).values
        dE=pd.DataFrame(E1-E0,columns=[str(State_A_Lambda_float.values[0])+"_"+str(State_B_Lambda_float.values[0])+"-"+str(State_A_Lambda_float.values[1])+"_"+str(State_B_Lambda_float.values[1]),str(State_A_Lambda_float.values[1])+"_"+str(State_B_Lambda_float.values[1])+"-"+str(State_A_Lambda_float.values[0])+"_"+str(State_B_Lambda_float.values[0])])
        dEs=dEs.append(dE.transpose())
    
    return dEs.transpose()

#%%
def Zwnazig_Estimator(dEs_df):
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

    dG_average_raw=((pd.DataFrame(dGF[1:]))-pd.DataFrame(dGR[1:][::-1]))/2

    for i in range(len(list(dG_average_raw.values))):
        dG_Average.append(np.sum(dG_average_raw.values[:i+1]))

    Zwanzig_df=pd.DataFrame.from_dict({"Lambda":Lambdas,"dG_Forward":dGF,"SUM_dG_Forward":dGF_sum,"dG_Reverse":dGR[::-1],"SUM_dG_Reverse":dGR_sum[::-1],"dG_Average":dG_Average})
    Zwanzig_Final_dG = Zwanzig_df['dG_Average'].iloc[-1]
    return Zwanzig_df, Zwanzig_Final_dG



def Plot_dG(df):
    p=plt.plot(df.iloc[:,2],'.',label= "ΔGf")
    p=plt.plot(df.iloc[:,4][::-1],'.',label ="ΔGr")
    plt.title('Hysteresis between ΔGf and ΔGr',fontsize=16)
    plt.xlabel("λ",fontsize=14)
    plt.ylabel("ΔG FEP (Kcal/mol)",fontsize=14)
    plt.legend()
    plt.savefig('Hysteresis.png',dpi=200)
#%%

#if '__name__' == '__main__':

#%%
    #os.chdir("/Users/nour/New_qfep/qfep_small") #MAC
    os.chdir("Z:/jobs/Qfep_NEW/qfep_small/test")
    #os.chdir("G:/PhD/Project/En")
    EnergyFiles_Lst = [filename for filename in glob.glob("FEP*.en")]  
    State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = ReadBinary(EnergyFiles_Lst)
    State_A_df = createDataFrames(State_A_RawEnergies_Lst)
    State_B_df = createDataFrames(State_B_RawEnergies_Lst)
    dEs =  dE_Calculation2(None)
    Zwanzig_df, Zwanzig_Final_dG= Zwnazig_Estimator(dEs)
    Plot_dG(Zwanzig_df)



##dEs
#%%

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

dG_average_raw=((pd.DataFrame(dGF[1:]))-pd.DataFrame(dGR[1:][::-1]))/2

for i in range(len(list(dG_average_raw.values))):
    dG_Average.append(np.sum(dG_average_raw.values[:i+1]))


#Zwanzig_df=pd.DataFrame.from_dict({"Lambda":Lambdas,"dG_Forward":dGF,"SUM_dG_Forward":dGF_sum,"dG_Reverse":dGR[::-1]})
Zwanzig_df=pd.DataFrame.from_dict({"Lambda":Lambdas,"dG_Forward":dGF,"SUM_dG_Forward":dGF_sum,"dG_Reverse":dGR[::-1],"SUM_dG_Reverse":dGR_sum[::-1],"dG_Average":dG_Average})

Zwanzig_df['dG_Average'].iloc[-1]

# %%
os.chdir("Z:/jobs/Qfep_NEW/qfep_small/test")
State_A_RawEnergies = []
State_B_RawEnergies = []
EnergyFiles_Lst = [filename for filename in glob.glob("FEP*.en")]  

for file in EnergyFiles_Lst:
    with open(file,'rb') as f: fileContent = f.read()

    Version_Check_Str=str(list(struct.unpack("c" * ((len(fileContent[32:112]))//1),fileContent[32:112]))).replace("b'", "").strip("[],' '").replace("'","").replace(",","").replace(" ","")

    EnergyFileLength_int=len(fileContent)
    BinaryChankSize_int=120
    NextBinaryChank_int=272
    HeaderSize_int=124
    State_B_Shift_int=132
    binary_structre=15*"d"+6*"h"+15*"d"+6*"h"
    steps=int((len(fileContent))/388))
    ## Check Q_Energies Version !!!

    if "6." in Version_Check_Str: HeaderSize_int += 4
    
    if not '5.' in Version_Check_Str and not "6." in Version_Check_Str:
        print("Pleaes Check the your Qdyn verion in file: "+file +" ----> format is NOT Supported !!! ")
        exit()
        


    for Byte in range(HeaderSize_int, EnergyFileLength_int, NextBinaryChank_int):

        struct.unpack("="+(binary_structre* steps),fileContent)
        State_A_RawEnergies.append(State_A_Lst)
        State_B_RawEnergies.append(State_B_Lst)
# %%
#YESSSSSSSSSSSSSSSSS
os.chdir("Z:/jobs/Qfep_NEW/qfep_small/test")
State_A_RawEnergies = []
State_B_RawEnergies = []
EnergyFiles_Lst = [filename for filename in glob.glob("FEP*.en")]  

for file in EnergyFiles_Lst:
    with open(file,'rb') as f: fileContent = f.read()
    binary_structre=15*"d"+6*"h"+15*"d"+10*"h"
    steps=int((len(fileContent)-116)/272)-1
    x=struct.unpack("="+(binary_structre* steps),fileContent[124:-264])
print(x)
##########################################
#%%

struct.unpack(15* "d",fileContent[124:244])
struct.unpack(6* "h",fileContent[244:256])
struct.unpack(15* "d",fileContent[256:376])
struct.unpack(10* "h",fileContent[376:396])
struct.unpack(15* "d",fileContent[396:516])
struct.unpack(6* "h",fileContent[516:528])
struct.unpack(15* "d",fileContent[528:648])
struct.unpack(6* "h",fileContent[648:])

# %%
fileContent=(388-116)*5
steps=((fileContent)/272)
steps
# %%
x= list(range(0,388)
len(x[124:244])
x[244:256]
len(x[4:32])

struct.calcsize(fileContent[4:32])
    binary_structre=2*"h"+7*"i"+80*"c"+6*"h"+15*"d"+6*"h"+15*"d"+6*"h"
