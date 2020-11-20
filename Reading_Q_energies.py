import struct
import numpy as np
import itertools
import os
import glob
import pandas as pd
from pandas.core.arrays.sparse import dtype

myfile = '/Users/nour/New_qfep/Q6/FEP1_md_1000_0000.en'

                


# def parseInputChunkEnergy(state):

#     return [item for item in state[]]

    # lambda_a.append(state_a[0]),q_sum_a.append(state_a[1]),q_bond_a.append(state_a[2])
    # q_angle_a.append(state_a[3]),q_torsion_a.append(state_a[4]),q_improper_a.append(state_a[5])
    # q_any_ele_a.append(state_a[6]),q_any_vdw_a.append(state_a[7])
    # q_q_ele_a.append(state_a[8]),q_q_vdw_a.append(state_a[9])
    # q_protein_ele_a.append(state_a[10]),q_protein_vdw_a.append(state_a[11])
    # q_water_ele_a.append(state_a[12]),q_water_vdw_a.append(state_a[13]),restraints_a.append(state_a[14])


with open(myfile, "rb") as f:
    byte_list=[]
    byte = f.read(1)
    while byte != b"":
        # Do stuff with byte.
        byte = f.read(1)
        print (byte)
        byte_list.append(byte)
    print(byte_list)



canary, fileheader_arrays, fileheader_totresid,fileheader_types,fleheader_numres,fileheader_resid,fileheader_gcnum,fileheader_version


###header
with open(myfile, "rb") as file:
    fileContent = file.read()
    canary, fileheader_arrays, fileheader_totresid,fileheader_types,fleheader_numres,fileheader_resid,fileheader_gcnum=struct.unpack("i" * ((len(fileContent)-4)//4),fileContent[4:])
    z=struct.unpack("f",fileContent[:4])
    z

with open(myfile, "rb") as file:

    fileContent = file.read()
    canary, fileheader_arrays, fileheader_totresid,fileheader_types,fleheader_numres,fileheader_resid,fileheader_gcnum=struct.unpack("i" * ((len(fileContent[4:32]))//4),fileContent[4:32])
    unk1=struct.unpack("hh",fileContent[:4])
    header=struct.unpack("i" * ((len(fileContent[4:32]))//4),fileContent[4:32])
    version=struct.unpack("c" * ((len(fileContent[32:116]))//1),fileContent[32:116])
    between_version_a=struct.unpack("h" * ((len(fileContent[112:124]))//2),fileContent[112:124])
    state_a=struct.unpack("d" * ((len(fileContent[124:244]))//8),fileContent[124:244])
    between_a_b=struct.unpack("h" * ((len(fileContent[244:256]))//2),fileContent[244:256])
    state_b=struct.unpack("d" * ((len(fileContent[256:-12]))//8),fileContent[256:-12])
    unk2=struct.unpack("hhhhhh",fileContent[-12:])

    lambda_a=[]
    lambda_b=[]
    q_sum_a=[]
    q_sum_b=[]
    q_bond_a=[]
    q_bond_b=[]
    q_angle_a=[]
    q_angle_b=[]
    q_torsion_a=[]
    q_torsion_b=[]
    q_improper_a=[]
    q_improper_b=[]
    q_any_ele_a=[]
    q_any_ele_b=[]
    q_any_vdw_a=[]
    q_any_vdw_b=[]
    q_q_ele_a=[]
    q_q_ele_b=[]
    q_q_vdw_a=[]
    q_q_vdw_b=[]
    q_protein_ele_a=[]
    q_protein_ele_b=[]
    q_protein_vdw_a=[]
    q_protein_vdw_b=[]
    q_water_ele_a=[]
    q_water_ele_b=[]
    q_water_vdw_a=[]
    q_water_vdw_b=[]
    restraints_a=[]
    restraints_b=[]


    stateBstartsIndjvfbnfkb = 132 
    headerSizeInt = 124
    chunkSizesInt = 120

    state_A_Lst = []
    state_B_Lst = []

    energy_files = [filename for filename in glob.glob("*.en")]  

    for file in energy_files:
        
        with open(file,'rb') as f: fileContent = f.read()

        for bin in range(124, len(fileContent), 272) :

            ## State A
            state_a=struct.unpack("d" * ((len(fileContent[bin:bin+120]))//8),fileContent[bin:bin+120])
            lambda_a.append(state_a[0]),q_sum_a.append(state_a[1]),q_bond_a.append(state_a[2])
            q_angle_a.append(state_a[3]),q_torsion_a.append(state_a[4]),q_improper_a.append(state_a[5])
            q_any_ele_a.append(state_a[6]),q_any_vdw_a.append(state_a[7])
            q_q_ele_a.append(state_a[8]),q_q_vdw_a.append(state_a[9])
            q_protein_ele_a.append(state_a[10]),q_protein_vdw_a.append(state_a[11])
            q_water_ele_a.append(state_a[12]),q_water_vdw_a.append(state_a[13]),restraints_a.append(state_a[14])
            ## State B
            state_b=struct.unpack("d" * ((len(fileContent[bin+132:bin+132+120]))//8),fileContent[bin+132:bin+132+120])
            lambda_b.append(state_b[0]),q_sum_b.append(state_b[1]),q_bond_b.append(state_b[2])
            q_angle_b.append(state_b[3]),q_torsion_b.append(state_b[4]),q_improper_b.append(state_b[5])
            q_any_ele_b.append(state_b[6]),q_any_vdw_b.append(state_b[7])
            q_q_ele_b.append(state_b[8]),q_q_vdw_b.append(state_b[9])
            q_protein_ele_b.append(state_b[10]),q_protein_vdw_b.append(state_b[11])
            q_water_ele_b.append(state_b[12]),q_water_vdw_b.append(state_b[13]),restraints_b.append(state_b[14])    



Energies_dict= {"State_A_Lambda":lambda_a,"State_A_G":q_sum_a,"State_B_Lambda":lambda_b,"State_B_G":q_sum_b}

Energies_df=pd.DataFrame(Energies_dict)

Energies_df["dG"]= Energies_df["State_A_G"]-Energies_df["State_B_G"]
with open(myfile, 'rb') as file:
    for line in itertools.islice(file, 33, -1):
        print(line)
###body
with open(myfile, 'rb') as file:
    for line in itertools.islice(file, 33, 34):
        print(struct.unpack("i",line))








with open(myfile, mode='rb') as file: # b is important -> binary
    fileContent = file.read(8)
 
    x=struct.unpack("d" * ((len(fileContent) -28) // 8), fileContent[20:-8])
    y=struct.unpack("d", fileContent[-8:])
x
y



z= np.fromfile(myfile, np.int32)
z

def parseQoutputFiles(args.qfiles):

    pass


def parseArgumments():

    'use argpar'

    pass

def main():

    args = parseArgumments()

    parseQoutputFiles(args.qfiles)




if '__name__' == '__main__':

    main()


###Q6

    version=struct.unpack("c" * ((len(fileContent[32:112]))//1),fileContent[32:112])
    between_version_a=struct.unpack("h" * ((len(fileContent[112:124]))//2),fileContent[112:124])
    state_a=struct.unpack("d" * ((len(fileContent[128:240]))//8),fileContent[128:240])
    state_a

    def ReadBinary(EnergyFiles_Lst):
    
    State_A_RawEnergies = []
    State_B_RawEnergies = []
#%%
with open(myfile,'rb') as f: fileContent = f.read()

EnergyFileLength_int=len(fileContent)
BinaryChankSize_int=120
NextBinaryChank_int=272
HeaderSize_int=116
State_B_Shift_int=132

# %% 
for Byte in range(HeaderSize_int, EnergyFileLength_int, NextBinaryChank_int):
    State_A_Lst=struct.unpack("d" * ((len(fileContent[Byte:Byte+BinaryChankSize_int]))//8),fileContent[Byte:Byte+BinaryChankSize_int])
   # State_B_Lst=struct.unpack("d" * ((len(fileContent[Byte+State_B_Shift_int:Byte+State_B_Shift_int+BinaryChankSize_int]))//8),fileContent[Byte+State_B_Shift_int:Byte+State_B_Shift_int+BinaryChankSize_int])
    
            
        
# %%
State_A_Lst
# %%
len(fileContent[Byte:Byte+BinaryChankSize_int])
# %%
 np.fromfile(myfile,dtype=np.float32)
# %%
version2= str(list(version))
# %%
version2 = [x.strip(' ') for x in version2]
# %%
x2=version2.replace("b'", "").strip("[],' '").replace("'","").replace(",",'').replace(" ","")
# %%
print(x2.find("5.",2,3))
# %%
