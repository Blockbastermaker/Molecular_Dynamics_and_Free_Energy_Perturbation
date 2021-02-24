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
from pandas.core.frame import DataFrame
import seaborn as sns
from scipy.stats import norm
from seaborn.utils import despine
import itertools

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


def Convergence(df,Estimator,StepsChunk_Int,ReplicatiesCount_Int,EnergyOutputInterval_Int):
                                    # the last and first steps are not included in the reading
    Zwanzig_Final_Lst=[Estimator(df,steps_limit)[1] for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,len(df)+1,StepsChunk_Int*ReplicatiesCount_Int)]
    StepsChunk_Lst=[EnergyOutputInterval_Int*steps_limit/ReplicatiesCount_Int for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,len(df)+1,StepsChunk_Int*ReplicatiesCount_Int)]
    Convergence_df=pd.DataFrame({'Number of Steps':StepsChunk_Lst, 'dG':Zwanzig_Final_Lst })
    return Convergence_df


def dU_Plot():
    dEs=pd.DataFrame()
    Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('State_A_Lambda')

    State_A_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list)),orient='index')
    State_A_Energies_df=State_A_Energies_df.transpose()


def Plot_Convergence(df):
    plt.plot(df['Number of Steps'],df['dG'])
    plt.title('Convergence Plot',fontsize=16)
    plt.xlabel("Number of Steps",fontsize=14)
    plt.ylabel("ΔG FEP (Kcal/mol)",fontsize=14)
    plt.savefig('Convergence.png',dpi=300)
    plt.close()
    
# def Plot_Hysteresis(df):
#     p=plt.plot(df.iloc[1:,2],'.',label= "ΔGf")
#     p=plt.plot(df.iloc[:-1,4],'.',label ="ΔGr")
#     plt.title('Hysteresis between ΔGf and ΔGr',fontsize=16)
#     plt.xlabel("λ",fontsize=14)
#     plt.ylabel("ΔG FEP (Kcal/mol)",fontsize=14)
#     plt.legend()
#     plt.savefig('Hysteresis.png',dpi=300)
#     plt.close()

def Plot_Hysteresis(df):
    p=plt.plot(df.iloc[:,0],df.iloc[:,2],'.',label= "ΔGf")
    p=plt.plot(df.iloc[:,0],df.iloc[:,4][::-1]*-1.0,'.',label ="ΔGr")
    plt.title('Hysteresis between ΔGf and ΔGr',fontsize=16)
    plt.xlabel("λ",fontsize=14)
    plt.tick_params(axis='x', which='major', labelsize=7)
    plt.xticks(rotation=60)
    plt.ylabel("ΔG FEP (Kcal/mol)",fontsize=14)
    plt.legend()
    plt.savefig('Hysteresis.png',dpi=300)
    plt.close()


def Plot_dG_by_Lambda(df):
    p=plt.plot(df.iloc[1:,0],df.iloc[1:,1],'.',label= "ΔGf")
    p=plt.plot(df.iloc[1:,0],df.iloc[:-1,3]*-1.0,'.',label ="ΔGr")
    plt.title('dG_vs_Lambda',fontsize=16)
    plt.xlabel("λ",fontsize=14)
    plt.tick_params(axis='x', which='major', labelsize=7)
    plt.xticks(rotation=60)
    plt.ylabel("ΔG FEP (Kcal/mol)",fontsize=14)
    plt.legend()
    plt.savefig('dG_vs_Lambda.png',dpi=300)
    plt.close()

def Plot_dEs(df):
    plots=df.reset_index().plot(x='index', y=list(df.columns)[:], kind = 'line', legend=True,
        subplots = True, layout=(int(len(df.columns)/2),2),sharex = True, figsize=(16, 14)).flatten()
    for plot in range(len(plots)):
        plots[plot].legend(loc='upper right',prop={'size': 7})
        plots[plot].set_xlabel('Steps (fs)',fontsize=20)
        plt.subplots_adjust(wspace=0.2,hspace =0.5)
    plt.suptitle('ΔEs Plots', fontsize=30) # Add the text/suptitle to figure
    plt.savefig('dEs.png',dpi=300)
    plt.close()
    
    
def Generate_PDF(df,axis,window1,color1,window2,color2):
    sns.distplot(df.iloc[:,window1].values , hist = False, kde = True,color=color1,
                kde_kws = {'shade': True,'alpha':0.4},label=df.columns[window1], ax=axis[window1])
    sns.distplot( df.iloc[:,window2].values , hist = False, kde = True,color=color2,
                kde_kws = {'shade': True,'alpha':0.4},label=df.columns[window2], ax=axis[window1])
    axis[window1].legend(loc='upper right',prop={'size': 7})
    plt.subplots_adjust(wspace=0.2,hspace = 0.5)

def Plot_PDF():

    Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"],"Window":State_A_df["Lambda"].astype(str)+"_"+State_B_df["Lambda"].astype(str)})).sort_values('State_A_Lambda')
    dU_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('Window',sort=False)['E'].apply(list)),orient='index')
    df=dU_df.transpose()
    f, axis = plt.subplots(int(len(df.columns)/2), 2, figsize=(10, 10))
    plt.subplots_adjust(wspace=0.2,hspace = 0.5)
    axis = axis.flatten()
    for i in range(1,len(df.columns[:-1])-1):
        Generate_PDF(df,axis,i,'orange',i+1 ,'gray')
    Generate_PDF(df,axis,0,'blue',1,'gray')
    Generate_PDF(df,axis,-1,'red',-2,'gray')
    [axis[i].set_xlabel('U (Kcal/mol)',fontsize=18) for i in [-1,-2] ]
    plt.suptitle('Probability Density Function of U', fontsize=20)
    plt.savefig('PDF.png',dpi=300)
    #plt.close()


def Plot_PDF_Matrix():
    Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"],"Window":State_A_df["Lambda"].astype(str)+"_"+State_B_df["Lambda"].astype(str)})).sort_values('State_A_Lambda')
    dU_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('Window',sort=False)['E'].apply(list)),orient='index')
    df=dU_df.transpose()

    f, axis = plt.subplots(int(len(df.columns)), int(len(df.columns)), figsize=(15, 15))
    for window1 in range(len(df.columns)):
        for window2 in range(len(df.columns)):
            if window1==window2: color1, color2='blue','blue'
            else: color1 ,color2='gray','orange'
            sns.distplot(df.iloc[:,window1].values , hist = False, kde = True,color=color1,
                kde_kws = {'shade': True,'alpha':0.4},label=df.columns[window1], ax=axis[window1,window2])
            sns.distplot( df.iloc[:,window2].values , hist = False, kde = True,color=color2,
                kde_kws = {'shade': True,'alpha':0.4},label=df.columns[window2], ax=axis[window1,window2])
            axis[window1,window2].legend(loc='upper right',prop={'size':4})
            axis[window1,window2].tick_params(labelsize=5)
            plt.subplots_adjust(wspace=0.5,hspace = 0.5)
    [axis[-1,-i].set_xlabel('U (Kcal/mol)',fontsize=5) for i in range(int(len(df.columns))) ]
    plt.suptitle('Probability Density Function Matrix', fontsize=25)
    plt.savefig('PDF_Matrix.png',dpi=300)
    #plt.close()
#%%

#if '__name__' == '__main__':

#%%
os.chdir("/Users/nour/New_qfep/qfep_small2") #MAC
#os.chdir("Z:/jobs/Qfep_NEW/")
#os.chdir("G:/PhD/Project/En")
EnergyFiles_Lst = [filename for filename in glob.glob("FEP1*.en")]  
State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = ReadBinary(EnergyFiles_Lst)
#State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = ReadAndCollectBinariesInParallel(EnergyFiles_Lst)
State_A_df = createDataFrames(State_A_RawEnergies_Lst)
State_B_df = createDataFrames(State_B_RawEnergies_Lst)
#State_A_Energies_df,State_B_Energies_df=dE_ParallelCalculationPrepare()
dEs =  dE_Calculation3()
#dEs =  Run_dE_ParallelCalculation(State_A_Energies_df,State_B_Energies_df)
Zwanzig_df, Zwanzig_Final_dG= Zwanazig_Estimator(dEs,None)


from . estimators import Zwanzig
zz.Zwanzig()
x.Zwanzig(dEs,None)

convergenc_df= Convergence(dEs,zz.Zwanzig,1000,1,10)
#print(convergenc_df)
Plot_PDF()
#fig=dEs.plot(subplots=True,figsize=(10,8),layout=(int(len(dEs.columns)/2), 3),sharex=True,legend=True)
#fig
#plt.close("all")

#convergenc_df=Convergence(dEs,Z)
#Plot_Convergence(convergenc_df)
#chunck=1000
#Zwanzig_Final_list=[Zwnazig_Estimator(dEs,steps)[1] for steps in range(0,len(dEs)+1,chunck)]
print(Zwanzig_Final_dG)
#Plot_Hysteresis(Zwanzig_df)
Plot_dG_by_Lambda(Zwanzig_df)



#%%
def TI_Estimator(State_A_df, State_B_df):
    dU_dH_df=(pd.DataFrame({"lambda":State_A_df["Lambda"],"fep":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('lambda')
    dU_dH_df.reset_index(drop=True,inplace=True)
    dU_dH_df.index.names = ['time']
    dU_dH_df.set_index(['lambda'], append=True,inplace=True)
    return TI


#%%
##### BAR Not ............. ready
dEs_df=dEs
dEs_list_F=[]
dEs_list_R=[]
Lambdas=set([re.split('_|-',dEs_df.columns[i-1])[0] for i in range(len(dEs_df.columns))])
for i in range(1,len(dEs_df.columns),2):
    dE_R=pd.DataFrame(columns=Lambdas)
    dE_F=pd.DataFrame(columns=Lambdas)
    #x=(re.split('_|-',dEs_df.columns[i-1])[0])
    dE_F[re.split('_|-',dEs_df.columns[i])[0]]=dEs_df.iloc[:,i]
    dE_F["fep-lambda"]=float(re.split('_|-',dEs_df.columns[i])[2])
    dE_R[re.split('_|-',dEs_df.columns[i-1])[0]]=dEs_df.iloc[:,i-1]
    dE_R["fep-lambda"]=float(re.split('_|-',dEs_df.columns[i-1])[2])
    dEs_list_R.append(dE_R)
    dEs_list_F.append(dE_F)
dEs_ready_F=pd.concat(dEs_list_F,ignore_index=False,sort=False)
#del dEs_ready['0.0']
dEs_ready_R=pd.concat(dEs_list_R,ignore_index=False,sort=False)

dEs_ready= dEs_ready_F.append(dEs_ready_R,ignore_index=False,sort=False)
x=dEs_ready['fep-lambda'].mode().values[0]

dfwf = dEs_ready_F[dEs_ready_F['fep-lambda'] == x]
dfwr = dEs_ready_R[dEs_ready_R['fep-lambda'] == x]
dfwf.fillna(dfwr,inplace=True)
dEs_ready_F = dEs_ready_F[dEs_ready_F['fep-lambda'] != x]
dEs_ready_F= dEs_ready_F.append(dfwf,ignore_index=False,sort=False)
dEs_ready_R = dEs_ready_R[dEs_ready_R['fep-lambda'] != x]
dEs_ready= dEs_ready_F.append(dEs_ready_R,ignore_index=False,sort=False)
dEs_ready = dEs_ready.reindex(sorted(dEs_ready.columns), axis=1)
dEs_ready.replace(np.nan, 0, inplace=True)

#dEs_ready=pd.concat([dEs_ready_R,dEs_ready], axis=0)
#del dEs_ready_R['1.0']
#dEs_ready.fillna(dEs_ready_R,inplace=True)
# replace zeroes in initial dataframe with nan
#dEs_ready.replace(0, np.nan, inplace=True)
# replace the nan values with the reverse dataframe --
# this should not overwrite any of the fwd work values
#dEs_ready[dEs_ready.isnull()] = dEs_ready_R
#dEs_ready.append(dEs_list_R[0],ignore_index=False,sort=False)
# replace remaining nan values back to zero
# dfw = dEs_ready_R#[dEs_ready_R['fep-lambda'] != 1.0]
# #dfw = dEs_ready_R[dEs_ready_R['fep-lambda'] == min(dEs_ready_R['fep-lambda'])]
# dEs_ready=dEs_ready.append(dfw,ignore_index=False,sort=False)
# #dEs_ready=pd.concat([dEs_ready,dfw],levels =['fep-lambda'], join="inner",ignore_index=False,sort=False)
# dEs_ready = dEs_ready.reindex(sorted(dEs_ready.columns), axis=1)

#dEs_ready.replace(np.nan, 0, inplace=True)
#dEs_ready=dEs_ready.groupby(['fep-lambda'], as_index=False).count()

dEs_ready.index=dEs_ready.index.astype(float)
dEs_ready.index.names = ['time']
dEs_ready.set_index(['fep-lambda'], append=True,inplace=True)
dEs_ready.columns=dEs_ready.columns.astype(float)
dEs_ready = dEs_ready.reindex(sorted(dEs_ready.columns), axis=1)
# sort final dataframe by `fep-lambda` (as opposed to `timestep`)
u_nk = dEs_ready.sort_index(level=dEs_ready.index.names[1:])

bar_vdw = BAR().fit(u_nk)
bar_vdw.delta_f_
bar_vdw.delta_f_.loc[0.00, 1.00]
#%%
from alchemlyb.estimators import BAR
bar_vdw = BAR().fit(u_nk)

bar_vdw.delta_f_
bar_vdw.delta_f_.loc[0.00, 1.00]

df1=dEs_ready.dropna(axis=1, how="all", thresh=None, subset=None, inplace=False)
dr1=dEs_ready_R.dropna(axis=1, how="all", thresh=None, subset=None, inplace=False)
pd.merge(dEs_ready, dfw, how="outer")
dEs_ready.append( dEs_ready_R)
del df1[0.0]
#%%
#dEs matrix reday!!!!!!!
#def dE_Calculation3():
Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('State_A_Lambda')

State_A_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list)),orient='index')
State_A_Energies_df=State_A_Energies_df.transpose()
State_B_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_B_Lambda',sort=False)['State_B_G'].apply(list)),orient="index") 
State_B_Energies_df=State_B_Energies_df.transpose()
lambdas_list_A=list(State_A_Energies_df.columns)
lambdas_list_B=list(State_B_Energies_df.columns)

time= [i for i in range(len(State_A_Energies_df))]
lambdas_df=[i for i in State_A_Energies_df.columns]
States={i:[] for i in range(len(lambdas_list_A))}
for i in range(len(State_A_Energies_df.columns)):
    State_A_Energies=State_A_Energies_df.iloc[:,[i]]
    State_A_Energies.columns=["0"]
    State_A_Lambda_float=State_A_Energies_df.columns[i]    
    
    State_B_Energies=State_B_Energies_df.iloc[:,[i]]
    State_B_Energies.columns=["0"]
    State_B_Lambda_float=State_B_Energies_df.columns[i]    
    E0=State_A_Energies*State_A_Lambda_float+State_B_Energies*State_B_Lambda_float
    for x in range(len(lambdas_list_A)):
        #print(State_A_Energies,State_B_Energies)
        #print("A: ",State_A_Lambda_float,'B:',State_B_Lambda_float,'X:' ,lambdas_list_A[x], 'X+1',lambdas_list_B[x])
        E1=State_A_Energies*lambdas_list_A[x]+State_B_Energies*lambdas_list_B[x]
        dE=E1-E0
        dElam=pd.DataFrame()
        dE=dE.values.tolist()
        dE=list(itertools.chain(*dE))
        States[i].append(dE)
        #print('lambda:',State_A_Lambda_float,lambdas_list_A[x],dE.values)
        #dE.columns=[State_A_Lambda_float]
        #dicts0[str(State_A_Lambda_float)]=list(dE.values)
        #print(dicts0[str(State_A_Lambda_float)])
        #chunks.append(pd.DataFrame(dicts0))
for i in range(len(States)):
    States[i]=list(itertools.chain(*States[i]))
dEx=pd.DataFrame.from_dict(States)
dEx.columns=lambdas_list_A
lambdas_df=lambdas_df*len(State_A_Energies_df)
lambdas_df.sort()
dEx['time']=time*len(State_A_Energies_df.columns)
dEx['fep-lambda']=lambdas_df
dEx=dEx.astype('float')
dEx.set_index(['time'] ,append=False,inplace=True)
dEx.set_index(['fep-lambda'], append=True,inplace=True)
dEx.columns= dEx.columns.astype('float')
#dEs= dEs*-0.592

#%%
from alchemlyb.estimators import BAR

bar_vdw = BAR().fit(dEx)

bar_vdw.delta_f_
bar_vdw.delta_f_.loc[0, 1]

#%%

def Create_df_BAR_MBAR_2(State_A_Energies_df,States_dicts,steps):
    
    time = [i for i in range(len(State_A_Energies_df))]
    lambdas_df=[i for i in State_A_Energies_df.columns]

    for x in States_dicts.keys():
        for i in range(len(States_dicts[x])):
            States_dicts[x][i]=States_dicts[x][i][:steps]
            
    for i in range(len(States_dicts)):
        States_dicts[i]=list(itertools.chain(*States_dicts[i]))
    u_nk_df=pd.DataFrame.from_dict(States_dicts)
    u_nk_df.columns=lambdas_list_A
    lambdas_df=lambdas_df*len(State_A_Energies_df.iloc[:steps])
    lambdas_df.sort()
    u_nk_df['time']=time[:steps]*len(State_A_Energies_df.columns)
    u_nk_df['fep-lambda']=lambdas_df
    u_nk_df=u_nk_df.astype('float')
    u_nk_df.set_index(['time'] ,append=False,inplace=True)
    u_nk_df.set_index(['fep-lambda'], append=True,inplace=True)
    u_nk_df.columns= dEx.columns.astype('float')
    u_nk_df.dropna(axis=0,inplace=True)
    return u_nk_df




#%%
ds=dEx
ds.transpose()




#%%

#def BAR_Estimator(State_A_df, State_B_df):

Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('State_A_Lambda')
State_A_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list)),orient='index')
State_A_Energies_df=State_A_Energies_df.transpose()
State_B_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_B_Lambda',sort=False)['State_B_G'].apply(list)),orient="index") 
State_B_Energies_df=State_B_Energies_df.transpose()
Es=pd.DataFrame()
Es=pd.DataFrame(columns=State_A_Energies_df.columns)
State_B_Energies_df.columns=list(State_A_Energies_df.columns.values)
Es_list=[]
for i in range(len(State_A_Energies_df.columns)):
    E=State_A_Energies_df.columns[i]*State_A_Energies_df+(1-State_A_Energies_df.columns[i])*State_B_Energies_df
    E["fep-lambda"]=State_A_Energies_df.columns[i]
    Es_list.append(E)
    print(E)
Es=pd.concat(Es_list,ignore_index=False,sort=False)
Es.index=Es.index.astype(float)
Es.index.names = ['time']
Es.set_index(['fep-lambda'], append=True,inplace=True)
Es.columns=Es.columns.astype(float)
Es= Es*-0.592
Es


#%%
##2


#%%
##3
Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('State_A_Lambda')
State_A_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list)),orient='index')
State_A_Energies_df=State_A_Energies_df.transpose()
State_B_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_B_Lambda',sort=False)['State_B_G'].apply(list)),orient="index") 
State_B_Energies_df=State_B_Energies_df.transpose()
Es=pd.DataFrame()
Es=pd.DataFrame(columns=State_A_Energies_df.columns)
State_B_Energies_df=State_B_Energies_df.columns*State_B_Energies_df
State_A_Energies_df=State_A_Energies_df.columns*State_A_Energies_df
State_B_Energies_df.columns=list(State_A_Energies_df.columns.values)
Es=State_A_Energies_df+ State_B_Energies_df
dEs_list=[]
dEs_ready=pd.DataFrame(columns=State_A_Energies_df.columns.values)
for i in range(len(Es.columns)):
    print('HHHHH', i)
    for x in range(len(Es.columns)):
        dE=pd.DataFrame(columns=State_A_Energies_df.columns.values)
        dE[x]=(Es.loc[x].values-Es.loc[i].values)
        dE["fep-lambda"]=Es.columns[i]
        dEs_list.append(dE)
dEs_ready=pd.concat(dEs_list,ignore_index=False,sort=False)
dEs_ready.index=dEs_ready.index.astype(float)
dEs_ready.index.names = ['time']
dEs_ready.set_index(['fep-lambda'], append=True,inplace=True)
dEs_ready.columns=dEs_ready.columns.astype(float)
dEs


#%%
dEx=pd.read_csv("E:\\de3.csv")
dEx=dEx.astype('float')
dEx.set_index(['time'] ,append=False,inplace=True)
dEx.set_index(['fep-lambda'], append=True,inplace=True)
dEx.columns= dEx.columns.astype('float')
#dEs= dEs*-0.592
dEx.columns
#%%
from alchemlyb.estimators import BAR

bar_vdw = BAR().fit(dEx)

bar_vdw.delta_f_
bar_vdw.delta_f_.loc[0.00, 1.00]
#%%
dEs_ready.to_csv("E:\\de2.csv",index=True)

#%%
def Plot_PDF_Matrix():
    Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"],"Window":State_A_df["Lambda"].astype(str)+"_"+State_B_df["Lambda"].astype(str)})).sort_values('State_A_Lambda')
    dU_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('Window',sort=False)['E'].apply(list)),orient='index')
    df=dU_df.transpose()

    f, axis = plt.subplots(int(len(df.columns)), int(len(df.columns)), figsize=(15, 15))
    for window1 in range(len(df.columns)):
        for window2 in range(len(df.columns)):
            if window1==window2: color1, color2='blue','blue'
            else: color1 ,color2='gray','orange'
            sns.distplot(df.iloc[:,window1].values , hist = False, kde = True,color=color1,
                kde_kws = {'shade': True,'alpha':0.4},label=df.columns[window1], ax=axis[window1,window2])
            sns.distplot( df.iloc[:,window2].values , hist = False, kde = True,color=color2,
                kde_kws = {'shade': True,'alpha':0.4},label=df.columns[window2], ax=axis[window1,window2])
            axis[window1,window2].legend(loc='upper right',prop={'size':4})
            axis[window1,window2].tick_params(labelsize=5)
            plt.subplots_adjust(wspace=0.5,hspace = 0.5)
    [axis[-1,-i].set_xlabel('U (Kcal/mol)',fontsize=5) for i in range(int(len(df.columns))) ]
    plt.suptitle('Probability Density Function Matrix', fontsize=25)
    #plt.savefig('PDF_Matrix.png',dpi=300)
    #plt.close()

#%%
parser = argparse.ArgumentParser(description="MD/FEP Analysis")

parser.add_argument("-f","--energy_files_prefix", help = "Energy Files Prefix, ex: FEP1, FEP2")

parser.add_argument("-a","--all_replicaties",default = False, action="store_true", help = "Analyze all replicaties in the current directory")

parser.add_argument("-p","--run_in_parallel",action="store_true", help = "Run in parallel")

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
        convergenc_df=Convergence(dEs,eval(args.convergence_analysis[0]),int(args.convergence_analysis[1]),int(args.convergence_analysis[2]),10)
        print(convergenc_df)
        Plot_Convergence(convergenc_df)

    if args.plot ==True:
        Plot_Hysteresis(Zwanzig_df)
        Plot_dG_by_Lambda(Zwanzig_df)
        Plot_dEs(dEs)
        Plot_PDF()
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


"""Thermodynamic integration (TI).
    Parameters
    ----------
    verbose : bool, optional
        Set to True if verbose debug output is desired.
    Attributes
    ----------
    delta_f_ : DataFrame
        The estimated dimensionless free energy difference between each state.
    d_delta_f_ : DataFrame
        The estimated statistical uncertainty (one standard deviation) in 
        dimensionless free energy differences.
    states_ : list
        Lambda states for which free energy differences were obtained.
    """

def __init__( verbose=False):
        verbose = verbose

def fit(dHdl):
        """
        Compute free energy differences between each state by integrating
        dHdl across lambda values.
        Parameters
        ----------
        dHdl : DataFrame 
            dHdl[n,k] is the potential energy gradient with respect to lambda
            for each configuration n and lambda k.
        """

        # sort by state so that rows from same state are in contiguous blocks,
        # and adjacent states are next to each other
        dHdl = dHdl.sort_index(level=dHdl.index.names[1:])

        # obtain the mean and variance of the mean for each state
        # variance calculation assumes no correlation between points
        # used to calculate mean
        means = dHdl.mean(level=dHdl.index.names[1:])
        variances = np.square(dHdl.sem(level=dHdl.index.names[1:]))
        
        # get the lambda names
        l_types = dHdl.index.names[1:]

        # obtain vector of delta lambdas between each state
        dl = means.reset_index()[means.index.names[:]].diff().iloc[1:].values

        # apply trapezoid rule to obtain DF between each adjacent state
        deltas = (dl * (means.iloc[:-1].values + means.iloc[1:].values)/2).sum(axis=1)

        # build matrix of deltas between each state
        adelta = np.zeros((len(deltas)+1, len(deltas)+1))
        ad_delta = np.zeros_like(adelta)

        for j in range(len(deltas)):
            out = []
            dout = []
            for i in range(len(deltas) - j):
                out.append(deltas[i] + deltas[i+1:i+j+1].sum())

                # Define additional zero lambda
                a = [0.0] * len(l_types)

                # Define dl series' with additional zero lambda on the left and right
                dll = np.insert(dl[i:i + j + 1], 0, [a], axis=0)
                dlr = np.append(dl[i:i + j + 1], [a], axis=0)

                # Get a series of the form: x1, x1 + x2, ..., x(n-1) + x(n), x(n)
                dllr = dll + dlr

                # Append deviation of free energy difference between state i and i+j+1
                dout.append((dllr ** 2 * variances.iloc[i:i + j + 2].values / 4).sum(axis=1).sum())
            adelta += np.diagflat(np.array(out), k=j+1)
            ad_delta += np.diagflat(np.array(dout), k=j+1)

        # yield standard delta_f_ free energies between each state
        delta_f_ = pd.DataFrame(adelta - adelta.T,
                                     columns=means.index.values,
                                     index=means.index.values)

        # yield standard deviation d_delta_f_ between each state
        d_delta_f_ = pd.DataFrame(np.sqrt(ad_delta + ad_delta.T),
                                       columns=variances.index.values,
                                       index=variances.index.values)
        globals()["dw"] =delta_f_
        states_ = means.index.values.tolist()
        print(states_)
        print(delta_f_)
        print( delta_f_.loc[0.00, 1.00])
        return 
    
fit(dU_dH_df)
# %%

import numpy as np
import pandas as pd

from sklearn.base import BaseEstimator

from pymbar import MBAR as MBAR_


class MBAR(BaseEstimator):
    """Multi-state Bennett acceptance ratio (MBAR).
    Parameters
    ----------
    maximum_iterations : int, optional
        Set to limit the maximum number of iterations performed.
    relative_tolerance : float, optional
        Set to determine the relative tolerance convergence criteria.
    initial_f_k : np.ndarray, float, shape=(K), optional
        Set to the initial dimensionless free energies to use as a 
        guess (default None, which sets all f_k = 0).
    method : str, optional, default="hybr"
        The optimization routine to use.  This can be any of the methods
        available via scipy.optimize.minimize() or scipy.optimize.root().
    verbose : bool, optional
        Set to True if verbose debug output is desired.
    Attributes
    ----------
    delta_f_ : DataFrame
        The estimated dimensionless free energy difference between each state.
    d_delta_f_ : DataFrame
        The estimated statistical uncertainty (one standard deviation) in
        dimensionless free energy differences.
    theta_ : DataFrame
        The theta matrix.
    states_ : list
        Lambda states for which free energy differences were obtained.
    """

    def __init__(self, maximum_iterations=10000, relative_tolerance=0,
                 initial_f_k=None, method='hybr', verbose=False):

        self.maximum_iterations = maximum_iterations
        self.relative_tolerance = relative_tolerance
        self.initial_f_k = initial_f_k
        self.method = [dict(method=method)]
        self.verbose = verbose

        # handle for pymbar.MBAR object
        self._mbar = None

    def fit(self, u_nk):
        """
        Compute overlap matrix of reduced potentials using multi-state
        Bennett acceptance ratio.
        Parameters
        ----------
        u_nk : DataFrame 
            u_nk[n,k] is the reduced potential energy of uncorrelated
            configuration n evaluated at state k.
        """
        # sort by state so that rows from same state are in contiguous blocks
        u_nk = u_nk.sort_index(level=u_nk.index.names[1:])
        
        groups = u_nk.groupby(level=u_nk.index.names[1:])
        print(u_nk.groupby(level=u_nk.index.names[1:]))
        N_k = [(len(groups.get_group(i)) if i in groups.groups else 0) for i in u_nk.columns]        
        print([(len(groups.get_group(i)) if i in groups.groups else 0) for i in u_nk.columns] ) 
        self._mbar = MBAR_(u_nk.T, N_k,
                           maximum_iterations=self.maximum_iterations,
                           relative_tolerance=self.relative_tolerance,
                           initial_f_k=self.initial_f_k,
                           solver_protocol=self.method,
                           verbose=self.verbose)

        self.states_ = u_nk.columns.values.tolist()

        # set attributes
        out = self._mbar.getFreeEnergyDifferences(return_theta=True)
        attrs = [pd.DataFrame(i,
                              columns=self.states_,
                              index=self.states_) for i in out]

        (self.delta_f_, self.d_delta_f_, self.theta_) = attrs 
        
        return self

    def predict(self, u_ln):
        pass

# we could also just call the `fit` method
# directly, since it returns the `MBAR` object
mbar_vdw = MBAR().fit(Es)

mbar_vdw.delta_f_
mbar_vdw.delta_f_.loc[0.00, 1.00]

# %%

# %%


# %%
def dEs_matrix(State_A_Lambda,State_B_Lambda):
    dEs_matrix=pd.DataFrame()
    Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('State_A_Lambda')

    State_A_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list)),orient='index')
    State_A_Energies_df=State_A_Energies_df.transpose()
    State_B_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_B_Lambda',sort=False)['State_B_G'].apply(list)),orient="index") 
    State_B_Energies_df=State_B_Energies_df.transpose()
    lambdas_list_A=list(State_A_Energies_df.columns)
    lambdas_list_B=list(State_B_Energies_df.columns)

    time= [i for i in range(len(State_A_Energies_df))]
    lambdas_df=[i for i in State_A_Energies_df.columns]
    States={i:[] for i in range(len(lambdas_list_A))}
    for i in range(len(State_A_Energies_df.columns)):
        State_A_Energies=State_A_Energies_df.iloc[:,[i]]
        State_A_Energies.columns=["0"]
        State_A_Lambda_float=State_A_Energies_df.columns[i]    
        
        State_B_Energies=State_B_Energies_df.iloc[:,[i]]
        State_B_Energies.columns=["0"]
        State_B_Lambda_float=State_B_Energies_df.columns[i]    
        E0=State_A_Energies*State_A_Lambda_float+State_B_Energies*State_B_Lambda_float
        for x in range(len(lambdas_list_A)):
            #print(State_A_Energies,State_B_Energies)
            #print("A: ",State_A_Lambda_float,'B:',State_B_Lambda_float,'X:' ,lambdas_list_A[x], 'X+1',lambdas_list_B[x])
            E1=State_A_Energies*lambdas_list_A[x]+State_B_Energies*lambdas_list_B[x]
            dE=E1-E0
            dE.columns=[str(State_A_Lambda_float)+"_"+str(lambdas_list_A[x])]
            dEs_matrix=pd.concat([dEs_matrix,dE],axis=1, sort=False)
    dEs_matrix= dEs_matrix.astype('float')
    # df = dEs_matrix.replace(0.0, np.nan)
    # df = df.dropna(how='all', axis=1)
    dEs_matrix=dEs_matrix.transpose()
    # dEs_matrix['State A']=[(re.split('_|-',i)[0]) for i in dEs_matrix.index.values]
    # dEs_matrix['State B']=[(re.split('_|-',i)[1]) for i in dEs_matrix.index.values]
    dEs_matrix.to_csv('dEs_matrix.csv', index=True)
# %%


# %%
    ###for dE matrix AI
def Zwanzig_matrix_AI(dEs,steps)
    dEs_df=pd.DataFrame(-0.592*np.log(np.mean(np.exp(-dEs.iloc[:None]/0.592))))
    Lambdas_F=[]
    Lambdas_R=[]
    Lambdas=[]
    dGF=[]
    dGF_sum=[]
    dGR=[]
    dGR_sum=[]
    dG_Average=[]
    dGR.append(0.0)
    dG_Average.append(0.0)
    Lambdas_F.append((re.split('_|-',dEs_df.index[-1])[0])+'_'+(re.split('_|-',dEs_df.index[-1])[0]))
    for i in range(1,len(dEs_df.index),2):
        Lambdas.append(re.split('_|-',dEs_df.index[i-1])[1])
        Lambdas_R.append((re.split('_|-',dEs_df.index[i])[1])+'_'+(re.split('_|-',dEs_df.index[i])[3]))
        Lambdas_F.append((re.split('_|-',dEs_df.index[i-1])[1])+'_'+(re.split('_|-',dEs_df.index[i-1])[3]))
        dGF.append(dEs_df.iloc[i,0])
        dGR.append(dEs_df.iloc[i-1,0])
    Lambdas_R.append((re.split('_|-',dEs_df.index[-1])[1])+'_'+(re.split('_|-',dEs_df.index[-1])[1]))
    Lambdas.append(re.split('_|-',dEs_df.index[-1])[1])
    dGF.append(0.0)
    dGF=dGF[::-1]
    for i in range(len(dGF)):
        dGF_sum.append(sum(dGF[:i+1]))
        dGR_sum.append(sum(dGR[:i+1]))

    dG_average_raw=((pd.DataFrame(dGF[1:]))-pd.DataFrame(dGR[1:][::-1]))/2

    for i in range(len(list(dG_average_raw.values))):
        dG_Average.append(np.sum(dG_average_raw.values[:i+1]))

    Zwanzig_df=pd.DataFrame.from_dict({"Lambda":Lambdas,"dG_Forward":dGF,"Lambda_F":Lambdas_F,"SUM_dG_Forward":dGF_sum,"dG_Reverse":dGR[::-1],"Lambda_R":Lambdas_R,"SUM_dG_Reverse":dGR_sum[::-1],"dG_Average":dG_Average})
    Zwanzig_Final_dG = Zwanzig_df['dG_Average'].iloc[-1]
    Zwanzig_df.to_csv('Zwanzig_df_lambdas_F-R.csv')

# %%
df.index.unique()

# %%
