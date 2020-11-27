#%%
import struct
import numpy as np
import itertools
import os
import glob
from numpy.core.overrides import verify_matching_signatures
import pandas as pd


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

if '__name__' == '__main__':

#%%
    os.chdir("/Users/nour/New_qfep/qfep_small2") #MAC
    #os.chdir("Z:/jobs/Qfep_NEW/qfep_small")
    EnergyFiles_Lst = [filename for filename in glob.glob("*.en")]  
    State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = ReadBinary(EnergyFiles_Lst)
    State_A_df = createDataFrames(State_A_RawEnergies_Lst)
    State_B_df = createDataFrames(State_B_RawEnergies_Lst)





##Zwnazig
# %%

Zwanzig_exp= -0.592*(np.log(np.mean(np.exp((State_A_df["Lambda"]*State_A_df["Q_sum"] + State_B_df["Lambda"]*State_B_df["Q_sum"])*0.592))))
Zwanzig_exp



# %%
Energies_df=pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Lambda"]*State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Lambda"]*State_B_df["Q_sum"],"U":State_A_df["Q_sum"]-State_B_df["Q_sum"]})

# %%
-0.592*(np.log(np.mean(np.exp((-27.05200000000002)*-0.592))))

# %%

E0= State_A_df["Q_sum"]*State_A_df["Lambda"]+State_B_df["Q_sum"]*State_B_df["Lambda"]
E1= State_A_df["Q_sum"]*State_B_df["Lambda"]+State_B_df["Q_sum"]*State_A_df["Lambda"]


dE=E1-E0
dE= np.exp(-dE/0.592)

lambdas_dE=[]
for i in range(0,len(dE),int(len(dE)/2)):
    lambdas_dE.append(i)

zz = [sum(dE[i:i+int(len(dE)/2)])/len(dE[i:i+int(len(dE)/2)]) for i in lambdas_dE]
# for i in lambdas_dE:

#     z.(sum(dE[i:i+int(len(dE)/2)])/len(dE[i:i+int(len(dE)/2)]))
    
dg=-0.592*np.log(zz)

# %%


# %%
