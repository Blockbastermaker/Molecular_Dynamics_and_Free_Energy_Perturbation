#%%
import struct
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import transpose
import pandas as pd
import re
import concurrent.futures
import argparse

#%%


#%%
def ReadBinary(EnergyFiles_Lst):
    
    State_A_RawEnergies = []
    State_B_RawEnergies = []

    for file in EnergyFiles_Lst:
        with open(file,'rb') as f: fileContent = f.read()

        Version_Check_Str=str(list(struct.unpack("c" * ((len(fileContent[32:112]))//1),fileContent[32:112]))).replace("b'", "").strip("[],' '").replace("'","").replace(",","").replace(" ","")
        HeaderSize_int=124
        EnergyFileLength_int=len(fileContent)
        NextBinaryChank_int=272
        BinaryStructre_str=15*"d"+6*"h"+15*"d"+10*"h"
        EnergySteps_int=int((EnergyFileLength_int-HeaderSize_int+8)/NextBinaryChank_int)-1 ## +8 is between the headr and state A, -1 is to exclude the last step (has diffrent stucture in the end h*6 not h*10)
        StateUnpackedEnergiesLength_int=15
        StateA_UnpackedEnergiesStart_int=0
        StateB_UnpackedEnergiesStart_int=21
        UnpackedEnergiesNextState_int=46 
        

        ## Check Q_Energies Version !!!

        if "6." in Version_Check_Str: HeaderSize_int += 4
        elif not '5.' in Version_Check_Str: 
            print("Pleaes Check the your Qdyn verion in file: "+file +" ----> format is NOT Supported !!! ")
            exit()
            
  
        UnpackedEnergies_lst=struct.unpack("="+(BinaryStructre_str* EnergySteps_int),fileContent[HeaderSize_int:-264]) #-264 is to exclude the last step (has diffrent stucture in the end h*6 not h*10)
        
        State_A_Lst = [UnpackedEnergies_lst[i:(i + StateUnpackedEnergiesLength_int)] for i in range(StateA_UnpackedEnergiesStart_int, len(UnpackedEnergies_lst), UnpackedEnergiesNextState_int)]
        State_B_Lst = [UnpackedEnergies_lst[i:(i + StateUnpackedEnergiesLength_int)] for i in range(StateB_UnpackedEnergiesStart_int, len(UnpackedEnergies_lst), UnpackedEnergiesNextState_int)]
        for step in State_A_Lst :State_A_RawEnergies.append(step)
        for step in State_B_Lst :State_B_RawEnergies.append(step)
        
    return State_A_RawEnergies, State_B_RawEnergies

def ReadBinaryParallel(file):
        
        with open(file,'rb') as f: fileContent = f.read()

        Version_Check_Str=str(list(struct.unpack("c" * ((len(fileContent[32:112]))//1),fileContent[32:112]))).replace("b'", "").strip("[],' '").replace("'","").replace(",","").replace(" ","")
        HeaderSize_int=124
        EnergyFileLength_int=len(fileContent)
        NextBinaryChank_int=272
        BinaryStructre_str=15*"d"+6*"h"+15*"d"+10*"h"
        EnergySteps_int=int((EnergyFileLength_int-HeaderSize_int+8)/NextBinaryChank_int)-1 ## +8 is between the headr and state A, -1 is to exclude the last step (has diffrent stucture in the end h*6 not h*10)
        StateUnpackedEnergiesLength_int=15
        StateA_UnpackedEnergiesStart_int=0
        StateB_UnpackedEnergiesStart_int=21
        UnpackedEnergiesNextState_int=46 
        

        ## Check Q_Energies Version !!!

        if "6." in Version_Check_Str: HeaderSize_int += 4
        elif not '5.' in Version_Check_Str: 
            print("Pleaes Check the your Qdyn verion in file: "+file +" ----> format is NOT Supported !!! ")
            exit()
            
  
        UnpackedEnergies_lst=struct.unpack("="+(BinaryStructre_str* EnergySteps_int),fileContent[HeaderSize_int:-264]) #-264 is to exclude the last step (has diffrent stucture in the end h*6 not h*10)
        
        State_A_Lst = [UnpackedEnergies_lst[i:(i + StateUnpackedEnergiesLength_int)] for i in range(StateA_UnpackedEnergiesStart_int, len(UnpackedEnergies_lst), UnpackedEnergiesNextState_int)]
        State_B_Lst = [UnpackedEnergies_lst[i:(i + StateUnpackedEnergiesLength_int)] for i in range(StateB_UnpackedEnergiesStart_int, len(UnpackedEnergies_lst), UnpackedEnergiesNextState_int)]

        return State_A_Lst, State_B_Lst

def ReadAndCollectBinariesInParallel(EnergyFiles_Lst):
    State_A_RawEnergies = []
    State_B_RawEnergies = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
            ExtractedBinaries = [executor.submit(ReadBinaryParallel, file) for file in EnergyFiles_Lst]
            [State_A_RawEnergies.extend(y) for y in [ExtractedBinaries[x].result()[0] for x in range(len(ExtractedBinaries))]]
            [State_B_RawEnergies.extend(y) for y in [ExtractedBinaries[x].result()[1] for x in range(len(ExtractedBinaries))]]

    return State_A_RawEnergies, State_B_RawEnergies



def createDataFrames(rawEnergy):

    Columns_name=["Lambda","Q_sum","Q_bond","Q_angle","Q_torsion","Q_improper","Q_any_ele","Q_any_vdw","Q_Q_ele","Q_Q_vdw","Q_protein_ele","Q_protein_vdw","Q_water_ele","Q_water_vdw","Restraints"]

    return pd.DataFrame(rawEnergy, columns=Columns_name)


#%%
def dE_Calculation(steps):
    
    dEs=pd.DataFrame()
    Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('State_A_Lambda')
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
def dE_Calculation2():
    dEs=pd.DataFrame()
    Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('State_A_Lambda')

    State_A_Energies_df=pd.DataFrame(dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list)))
    State_B_Energies_df=pd.DataFrame(dict(Energies_df.groupby('State_B_Lambda',sort=False)['State_B_G'].apply(list))) 


    for i in range(len(State_A_Energies_df.columns)-1):
        print(i)
        State_A_Energies=State_A_Energies_df.iloc[:,[i,i+1]]
        State_A_Lambda_float=State_A_Energies.columns
        
        State_B_Energies=State_B_Energies_df.iloc[:,[i,i+1]]
        State_B_Lambda_float=State_B_Energies.columns
        
        E0=(State_A_Energies*State_A_Lambda_float).values+(State_B_Energies*State_B_Lambda_float).values
        E1=(State_A_Energies*State_A_Lambda_float[::-1]).values+(State_B_Energies*State_B_Lambda_float[::-1]).values

        dE=pd.DataFrame(E1-E0,columns=[str(State_A_Lambda_float.values[0])+"_"+str(State_B_Lambda_float.values[0])+"-"+str(State_A_Lambda_float.values[1])+"_"+str(State_B_Lambda_float.values[1]),str(State_A_Lambda_float.values[1])+"_"+str(State_B_Lambda_float.values[1])+"-"+str(State_A_Lambda_float.values[0])+"_"+str(State_B_Lambda_float.values[0])])
        dEs=dEs.append(dE)
    return dEs
#%%
def dE_Calculation3():
    dEs=pd.DataFrame()
    Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('State_A_Lambda')

    State_A_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list)),orient='index')
    State_A_Energies_df=State_A_Energies_df.transpose()
    State_B_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_B_Lambda',sort=False)['State_B_G'].apply(list)),orient="index") 
    State_B_Energies_df=State_B_Energies_df.transpose()

    for i in range(len(State_A_Energies_df.columns)-1):
        State_A_Energies=State_A_Energies_df.iloc[:,[i,i+1]]
        State_A_Energies.columns=["0","1"]
        State_A_Lambda_float=State_A_Energies_df.iloc[:,[i,i+1]].columns
        
        State_B_Energies=State_B_Energies_df.iloc[:,[i,i+1]]
        State_B_Energies.columns=["0","1"]
        State_B_Lambda_float=State_B_Energies_df.iloc[:,[i,i+1]].columns
        
        E0=State_A_Energies*State_A_Lambda_float+State_B_Energies*State_B_Lambda_float
        E1=State_A_Energies*State_A_Lambda_float[::-1]+State_B_Energies*State_B_Lambda_float[::-1]
        dE=E1-E0
        dE.columns=[str(State_A_Lambda_float.values[0])+"_"+str(State_B_Lambda_float.values[0])+"-"+str(State_A_Lambda_float.values[1])+"_"+str(State_B_Lambda_float.values[1]),str(State_A_Lambda_float.values[1])+"_"+str(State_B_Lambda_float.values[1])+"-"+str(State_A_Lambda_float.values[0])+"_"+str(State_B_Lambda_float.values[0])]
        dEs=pd.concat([dEs,dE],axis=1, sort=False)
    return dEs

#%%
def dE_ParallelCalculationPrepare():
    Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('State_A_Lambda')

    State_A_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list)),orient='index')
    State_A_Energies_df=State_A_Energies_df.transpose()
    State_B_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_B_Lambda',sort=False)['State_B_G'].apply(list)),orient="index") 
    State_B_Energies_df=State_B_Energies_df.transpose()
    return State_A_Energies_df,State_B_Energies_df

def dE_ParallelCalculation(State_A_Energies_df,State_B_Energies_df,LambdaState_Int):
    State_A_Energies=State_A_Energies_df.iloc[:,[LambdaState_Int,LambdaState_Int+1]]
    State_A_Energies.columns=["0","1"]
    State_A_Lambda_float=State_A_Energies_df.iloc[:,[LambdaState_Int,LambdaState_Int+1]].columns
    
    State_B_Energies=State_B_Energies_df.iloc[:,[LambdaState_Int,LambdaState_Int+1]]
    State_B_Energies.columns=["0","1"]
    State_B_Lambda_float=State_B_Energies_df.iloc[:,[LambdaState_Int,LambdaState_Int+1]].columns
    
    E0=State_A_Energies*State_A_Lambda_float+State_B_Energies*State_B_Lambda_float
    E1=State_A_Energies*State_A_Lambda_float[::-1]+State_B_Energies*State_B_Lambda_float[::-1]
    dE=E1-E0
    dE.columns=[str(State_A_Lambda_float.values[0])+"_"+str(State_B_Lambda_float.values[0])+"-"+str(State_A_Lambda_float.values[1])+"_"+str(State_B_Lambda_float.values[1]),str(State_A_Lambda_float.values[1])+"_"+str(State_B_Lambda_float.values[1])+"-"+str(State_A_Lambda_float.values[0])+"_"+str(State_B_Lambda_float.values[0])]
    return dE

def Run_dE_ParallelCalculation(State_A_Energies_df,State_B_Energies_df):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        dEs=pd.DataFrame()
        ###dE=[pd.concat([dEs,dE.result()],axis=1, sort=False) for dE in [executor.submit(dE_ParallelCalculation,State_A_Energies_df,State_B_Energies_df,i) for i in range(len(State_A_Energies_df.columns)-1)]]
        dE=[(dE.result()) for dE in [executor.submit(dE_ParallelCalculation,State_A_Energies_df,State_B_Energies_df,LambdaState_Int) for LambdaState_Int in range(len(State_A_Energies_df.columns)-1)]]
        for i in dE:dEs=pd.concat([dEs,i],axis=1, sort=False)
    return (dEs)

def Run_dE_ParallelCalculation2(State_A_Energies_df,State_B_Energies_df):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        dEs=pd.DataFrame()
        #dE=[(dE.result()) for dE in [executor.submit(dE_ParallelCalculation,State_A_Energies_df,State_B_Energies_df,LambdaState_Int) for LambdaState_Int in range(len(State_A_Energies_df.columns)-1)]]
        dE1=[executor.submit(dE_ParallelCalculation,State_A_Energies_df,State_B_Energies_df,LambdaState_Int) for LambdaState_Int in range(len(State_A_Energies_df.columns)-1)]
        dE=[(dE.result()) for dE in dE1]
        for i in dE:dEs=pd.concat([dEs,i],axis=1, sort=False)
    return (dEs)
#%%
def Zwanazig_Estimator(dEs_df,steps):
    dEs_df=pd.DataFrame(-0.592*np.log(np.mean(np.exp(-dEs.iloc[:steps]/0.592))))
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


def Convergence(df,Estimator,StepsChunk_Int,ReplicatiesCount_Int):
    StepsChunk_Int-=1 # the last step is not included in the reading
    Zwanzig_Final_Lst=[Estimator(df,steps_limit)[1] for steps_limit in range(StepsChunk_Int*ReplicatiesCount_Int,len(df)+1,StepsChunk_Int*ReplicatiesCount_Int)]
    StepsChunk_Lst=[steps_limit/ReplicatiesCount_Int for steps_limit in range(StepsChunk_Int*ReplicatiesCount_Int,len(df)+1,StepsChunk_Int*ReplicatiesCount_Int)]
    Convergence_df=pd.DataFrame({'Number of Steps':StepsChunk_Lst, 'dG':Zwanzig_Final_Lst })
    return Convergence_df

def Plot_Convergence(df):
    plt.plot(df['Number of Steps'],df['dG'])
    plt.title('Convergence Plot',fontsize=16)
    plt.xlabel("Number of Steps",fontsize=14)
    plt.ylabel("ΔG FEP (Kcal/mol)",fontsize=14)
    plt.savefig('Convergence.png',dpi=300)
    
def Plot_Hysteresis(df):
    p=plt.plot(df.iloc[:,2],'.',label= "ΔGf")
    p=plt.plot(df.iloc[:,4][::-1],'.',label ="ΔGr")
    plt.title('Hysteresis between ΔGf and ΔGr',fontsize=16)
    plt.xlabel("λ",fontsize=14)
    plt.ylabel("ΔG FEP (Kcal/mol)",fontsize=14)
    plt.legend()
    plt.savefig('Hysteresis.png',dpi=300)
    
def Plot_dG_by_Lambda(df):
    p=plt.plot(df.iloc[:,0],df.iloc[:,1],'.',label= "ΔGf")
    p=plt.plot(df.iloc[:,0],df.iloc[:,3][::-1],'.',label ="ΔGr")
    plt.title('dG_vs_Lambda',fontsize=16)
    plt.xlabel("λ",fontsize=14)
    plt.tick_params(axis='x', which='major', labelsize=7)
    plt.xticks(rotation=60)
    plt.ylabel("ΔG FEP (Kcal/mol)",fontsize=14)
    plt.legend()
    plt.savefig('dG_vs_Lambda.png',dpi=300)
#%%

#if '__name__' == '__main__':

#%%
#os.chdir("/Users/nour/New_qfep/") #MAC
os.chdir("Z:/jobs/Qfep_NEW/")
#os.chdir("G:/PhD/Project/En")
EnergyFiles_Lst = [filename for filename in glob.glob("FEP2*.en")]  
State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = ReadBinary(EnergyFiles_Lst)
#State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = ReadAndCollectBinariesInParallel(EnergyFiles_Lst)
State_A_df = createDataFrames(State_A_RawEnergies_Lst)
State_B_df = createDataFrames(State_B_RawEnergies_Lst)
#State_A_Energies_df,State_B_Energies_df=dE_ParallelCalculationPrepare()
dEs =  dE_Calculation(None)
dEs =  Run_dE_ParallelCalculation(State_A_Energies_df,State_B_Energies_df)
Zwanzig_df, Zwanzig_Final_dG= Zwanazig_Estimator(dEs,None)

convergenc_df= Convergence(dEs,Zwanazig_Estimator,1000,1)
print(convergenc_df)
#convergenc_df=Convergence(dEs,Z)
Plot_Convergence(convergenc_df)
#chunck=1000
#Zwanzig_Final_list=[Zwnazig_Estimator(dEs,steps)[1] for steps in range(0,len(dEs)+1,chunck)]
print(Zwanzig_Final_dG)
#Plot_Hysteresis(Zwanzig_df)
Plot_dG_by_Lambda(Zwanzig_df)

#%%
parser = argparse.ArgumentParser(description="MD/FEP Analysis")

parser.add_argument("-f","--energy_files_prefix", help = "Energy Files Prefix, ex: FEP1, FEP2")

parser.add_argument("-a","--all_replicaties",default = False, action="store_true", help = "Analyze all replicaties in the current directory")

parser.add_argument("-p","--run_in_parallel", help = "Run in parallel")

parser.add_argument("-e","--estimator",nargs='+', default='Zwanazig_Estimator',help = "Energy Estimator")

parser.add_argument("-c","--convergence_analysis",nargs='+', help = "Convergence Analysis: Estimator, by Number of Steps(fs), Number of used Replicaties")

parser.add_argument("-t","--plot", default = False, action="store_true", help = "Plot and Save")


args = parser.parse_args()
if __name__ == "__main__":
    if args.all_replicaties==True:
        EnergyFiles_Lst = [filename for filename in glob.glob('*/'+args.energy_files_prefix+'*.en')]
    else:
        EnergyFiles_Lst = [filename for filename in glob.glob(args.energy_files_prefix+'*.en')]
    
    if args.run_in_parallel==True:
        State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = ReadAndCollectBinariesInParallel(EnergyFiles_Lst)
        State_A_df = createDataFrames(State_A_RawEnergies_Lst)
        State_B_df = createDataFrames(State_B_RawEnergies_Lst)
        State_A_Energies_df,State_B_Energies_df=dE_ParallelCalculationPrepare()
        dEs =  Run_dE_ParallelCalculation(State_A_Energies_df,State_B_Energies_df)

    else:
        State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = ReadBinary(EnergyFiles_Lst)
        State_A_df = createDataFrames(State_A_RawEnergies_Lst)
        State_B_df = createDataFrames(State_B_RawEnergies_Lst)
        dEs =  dE_Calculation3()

    if args.estimator =='Zwanazig_Estimator':
        Zwanzig_df, Zwanzig_Final_dG= Zwanazig_Estimator(dEs,None)
        print(Zwanzig_Final_dG)
        
    if args.convergence_analysis is not None:
        args.convergence_analysis=args.convergence_analysis[0].split(',')
        convergenc_df=Convergence(dEs,eval(args.convergence_analysis[0]),int(args.convergence_analysis[1]),int(args.convergence_analysis[2]))
        print(convergenc_df)
        Plot_Convergence(convergenc_df)

    if args.plot ==True:
        Plot_Hysteresis(Zwanzig_df)
        Plot_dG_by_Lambda(Zwanzig_df)
else: 
    print("please use"+ " "+ "'MD/FEP Analysis.py -h'"+" ""for usege ")  



# %%
# from datetime import datetime
# start_time = datetime.now()
# dEs =  dE_Calculation(None)
# end_time = datetime.now()
# print('Duration: {}'.format(end_time - start_time))

x=[i.split('\\')[1] for i in EnergyFiles_Lst]
x.count([i.split('\\')[1] for i in EnergyFiles_Lst][0])


# %%
