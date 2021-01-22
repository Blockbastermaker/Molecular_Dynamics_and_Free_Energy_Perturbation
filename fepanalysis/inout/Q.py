import glob
import pandas as pd
import concurrent
import struct

import logging
logger = logging.getLogger(__name__)


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


def dE_ParallelCalculationPrepare(State_A_df, State_B_df):

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



def createDataFrames(rawEnergy):

    Columns_name=["Lambda","Q_sum","Q_bond","Q_angle","Q_torsion","Q_improper","Q_any_ele","Q_any_vdw","Q_Q_ele","Q_Q_vdw","Q_protein_ele","Q_protein_vdw","Q_water_ele","Q_water_vdw","Restraints"]

    return pd.DataFrame(rawEnergy, columns=Columns_name)


def ReadAndCollectBinariesInParallel(EnergyFiles_Lst):

    State_A_RawEnergies = []
    State_B_RawEnergies = []
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
            ExtractedBinaries = [executor.submit(ReadBinaryParallel, file) for file in EnergyFiles_Lst]
            [State_A_RawEnergies.extend(y) for y in [ExtractedBinaries[x].result()[0] for x in range(len(ExtractedBinaries))]]
            [State_B_RawEnergies.extend(y) for y in [ExtractedBinaries[x].result()[1] for x in range(len(ExtractedBinaries))]]

    return State_A_RawEnergies, State_B_RawEnergies

def Run_dE_ParallelCalculation(State_A_Energies_df,State_B_Energies_df):

    with concurrent.futures.ThreadPoolExecutor() as executor:

        dEs=pd.DataFrame()
        dE=[(dE.result()) for dE in [executor.submit(dE_ParallelCalculation,State_A_Energies_df,State_B_Energies_df,LambdaState_Int) for LambdaState_Int in range(len(State_A_Energies_df.columns)-1)]]

        for i in dE: 
            dEs=pd.concat([dEs,i],axis=1, sort=False)

    return (dEs)

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


def dE_Calculation3(State_A_df, State_B_df):

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


def parser(args):

    if args.all_replicaties==True:
        
        EnergyFiles_Lst = [filename for filename in glob.glob('*/'+args.energy_files_prefix+'*.en')]
    
    else:
    
        EnergyFiles_Lst = [filename for filename in glob.glob(args.energy_files_prefix+'*.en')]

    if len(EnergyFiles_Lst) !=0: logger.info('Parsing Q input Energy files: ' + ' '.join(EnergyFiles_Lst))
    else: logger.error('ERROR NOT energy files found in : ' + '*/'+args.energy_files_prefix+'*.en')


    if args.run_in_parallel==True:

        State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = ReadAndCollectBinariesInParallel(EnergyFiles_Lst)

        State_A_df = createDataFrames(State_A_RawEnergies_Lst)
        State_B_df = createDataFrames(State_B_RawEnergies_Lst)

        State_A_Energies_df, State_B_Energies_df = dE_ParallelCalculationPrepare(State_A_df, State_B_df)

        dEs = Run_dE_ParallelCalculation(State_A_Energies_df,State_B_Energies_df)

    else:

        State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = ReadBinary(EnergyFiles_Lst)
        
        State_A_df = createDataFrames(State_A_RawEnergies_Lst)
        State_B_df = createDataFrames(State_B_RawEnergies_Lst)

        dEs = dE_Calculation3(State_A_df, State_B_df)

    return dEs, State_A_df, State_B_df
