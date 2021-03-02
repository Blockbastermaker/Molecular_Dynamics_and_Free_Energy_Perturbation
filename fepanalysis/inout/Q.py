import glob
import pandas as pd
import concurrent
import struct
import os
import logging
import re
logger = logging.getLogger(__name__)

class Binary():
    """
    Read and exctract the energies from Q energies files.
    
    Read and extract the needed energies form the binray energy files predeuced by
    Q softwear package, ethier using single core or multicore processing.
    
    Returns two pandas dataFramas contains the energies for state A and state B.   
    
    """
    
    
    def ReadBinary(EnergyFiles_Lst):
        """
        Read and exctract the energies from Q energies files using one cpu core.
        
        Read and extract the needed energies form the binray energy files predeuced by
        Q softwear package, using one cpu core.
        
        Parameters
        ----------
        EnergyFiles_Lst : List
                    contains the names of the names of energy file in the needed directory.
        
        Returns
        ---------
        State_A_RawEnergies : List
                            Contains state A extracted energies.
                            
        State_B_RawEnergies : List
                            Contains state B extracted energies.
        """
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
        
        """
        Read and exctract the energies from Q energies files.
        
        Read and extract the needed energies form the binray energy files predeuced by
        Q softwear package, this funcation is used by the funcation ReadAndCollectBinariesInParallel
        to read and extact binaries using multiple prosseors.
        
        Parameters
        ----------
        EnergyFiles_Lst : List
                    contains the names of the names of energy file in the needed directory.
        
        Returns
        ---------
        State_A_Lst : List
                            Contains state A extracted energies.
                            
        State_B_Lst : List
                            Contains state B extracted energies.
        """
        
        
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
        """
        Run the Reading and exctraction of the binary energies using all available proseccores.
        
        Parameters
        ----------
        EnergyFiles_Lst : List
                    contains the names of the names of energy file in the needed directory.
        
        Returns
        ---------
        State_A_RawEnergies : List
                            Contains state A extracted energies.
                            
        State_B_RawEnergies : List
                            Contains state B extracted energies.
        """
        State_A_RawEnergies = []
        State_B_RawEnergies = []
        
        with concurrent.futures.ThreadPoolExecutor() as executor:
                ExtractedBinaries = [executor.submit(dE.ReadBinaryParallel, file) for file in EnergyFiles_Lst]
                [State_A_RawEnergies.extend(y) for y in [ExtractedBinaries[x].result()[0] for x in range(len(ExtractedBinaries))]]
                [State_B_RawEnergies.extend(y) for y in [ExtractedBinaries[x].result()[1] for x in range(len(ExtractedBinaries))]]

        return State_A_RawEnergies, State_B_RawEnergies
        
    def createDataFrames(rawEnergy):

        """
        Create pandas dataframe for any state energies.
        
        Parameters
        ----------
        rawEnergy : List
                Contains the extracted energies from the ReadBinary funcation.  
        
        Returns
        ---------
        Pandas DataFrame
                            
        """


        Columns_name=["Lambda","Q_sum","Q_bond","Q_angle","Q_torsion","Q_improper","Q_any_ele","Q_any_vdw","Q_Q_ele","Q_Q_vdw","Q_protein_ele","Q_protein_vdw","Q_water_ele","Q_water_vdw","Restraints"]

        return pd.DataFrame(rawEnergy, columns=Columns_name)
    
    
class dE():
    def dE_Calculation(State_A_df, State_B_df):
        """
        Calculates the energy diffrance (dE) between lambda states.
        
        Calculates the energy diffrance (dE) between lambda states in a forwared (0 to 1)
        and backwared (1 to 0) direcation using one CPU core.
        
        Parameters
        ----------
        State_A_df : Pandas DataFrame
                Contains state A extracted energies.  
                
        State_B_df : Pandas DataFrame
                Contains state B extracted energies.  
        
        Returns
        ---------
        dEs : Pandas DataFrame
                            
        
        """

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
    
    def dE_ParallelCalculationPrepare(State_A_df, State_B_df):
        """
        Create a Pandas dataframe needed for the function dE_ParallelCalculation. 

        Prepare energies Dataframe for energy differacne (dE) calculations using multiple proscessors.
        
        Parameters
        ----------
        State_A_df : Pandas DataFrame
                Contains state A extracted energies.  
                
        State_B_df : Pandas DataFrame
                Contains state B extracted energies.  
        
        Returns
        ---------
        State_A_Energies_df : Pandas DataFrame
        
        State_B_Energies_df : Pandas DataFrame
        
        """
        Energies_df=(pd.DataFrame({"State_A_Lambda":State_A_df["Lambda"],"State_A_G":State_A_df["Q_sum"] ,"State_B_Lambda":State_B_df["Lambda"],"State_B_G":State_B_df["Q_sum"],"E":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('State_A_Lambda')

        State_A_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_A_Lambda',sort=False)['State_A_G'].apply(list)),orient='index')
        State_A_Energies_df=State_A_Energies_df.transpose()
        State_B_Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('State_B_Lambda',sort=False)['State_B_G'].apply(list)),orient="index") 
        State_B_Energies_df=State_B_Energies_df.transpose()

        return State_A_Energies_df,State_B_Energies_df

    def dE_ParallelCalculation(State_A_Energies_df,State_B_Energies_df,LambdaState_Int):
        """
        Calculates the energy diffrance (dE) between lambda states.
        
        Calculates the energy diffrance (dE) between lambda states in a forwared (0 to 1)
        and backwared (1 to 0) direcation, and its a part form run_dE_ParallelCalculation
        to energy diffrance using multiple proseccores.
        
        Parameters
        ----------
        State_A_df : Pandas DataFrame
                Contains state A extracted energies.  
                
        State_B_df : Pandas DataFrame
                Contains state B extracted energies.  
        
        LambdaState_Int : Integer
                State A Lambda 
        Returns
        ---------
        dE : Pandas DataFrame
                            
        
        """

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
        
        """
        Run the calculation of the energy differacne between the states using all available proseccores. 
        
        Parameters
        ----------
        State_A_Energies_df : Pandas DataFrame
                Contains state A extracted and prepared energies by the dE_ParallelCalculationPrepare function.  
                
        State_A_Energies_df : Pandas DataFrame
                Contains state B extracted and prepared energies by the dE_ParallelCalculationPrepare function. 
        
        Returns
        ---------
        dEs : Pandas DataFrame
        
        """
        
        
        with concurrent.futures.ThreadPoolExecutor() as executor:

            dEs=pd.DataFrame()
            dE=[(dE.result()) for dE in [executor.submit(dE.dE_ParallelCalculation,State_A_Energies_df,State_B_Energies_df,LambdaState_Int) for LambdaState_Int in range(len(State_A_Energies_df.columns)-1)]]

            for i in dE: 
                dEs=pd.concat([dEs,i],axis=1, sort=False)

        return (dEs)    
    
    def dEs_matrix(State_A_df,State_B_df):
        """
        Development in Progress !!!!!

        """
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
        #dEs_matrix.to_csv('dEs_matrix.csv', index=True)    
        return dEs_matrix

    def Get_dEs_dGs_AI(Zwanzig_df,dEs_matrix):
        """
        Development in Progress !!!!!

        """

        dgf_dict = dict(zip(Zwanzig_df.Lambda_F,Zwanzig_df.dG_Forward))
        dgr_dict = dict(zip(Zwanzig_df.Lambda_R,Zwanzig_df.dG_Reverse))

        dgf_dict_matrix={}
        for  lam, dg in dgf_dict.items():
            for x in Zwanzig_df['Lambda'].values:
                a=float(re.split('_',lam)[1])
                b=float(x)
                if a<b:
                    nf= str(b)+'_'+str(a)
                    dgf_dict_matrix[nf]=dg

        dgr_dict_matrix={}
        for  lam, dg in dgr_dict.items():
            for x in Zwanzig_df['Lambda'].values:
                a=float(re.split('_',lam)[1])
                b=float(x)
                if a>b:
                    nr= str(b)+'_'+str(a)
                    dgr_dict_matrix[nr]=dg

        dg_matrix= dgf_dict_matrix.update(dgr_dict_matrix)
        dg_matrix_df=pd.DataFrame.from_dict(dgf_dict_matrix,orient='index')
        dg_matrix_df.sort_index(ascending=True,inplace=True)
        dg_matrix_df.columns=["dG"] 
        dEs_dGs_AI=pd.concat([dg_matrix_df, dEs_matrix], axis=1)
        dEs_dGs_AI.to_csv('dEs_dGs_AI.csv', index=True)


def parser(args):

    if args.all_replicaties==True:
        
        EnergyFiles_Lst = [filename for filename in glob.glob('*/'+args.energy_files_prefix+'*.en')]
    
    else:
    
        EnergyFiles_Lst = [filename for filename in glob.glob(args.energy_files_prefix+'*.en')]

    if len(EnergyFiles_Lst) !=0: logger.info('Parsing Q input Energy files: ' + ' '.join(EnergyFiles_Lst))
    else: logger.error('ERROR NOT energy files found in : ' + '*/'+args.energy_files_prefix+'*.en')


    if args.run_in_parallel==True:

        State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = Binary.ReadAndCollectBinariesInParallel(EnergyFiles_Lst)

        State_A_df = Binary.createDataFrames(State_A_RawEnergies_Lst)
        State_B_df = Binary.createDataFrames(State_B_RawEnergies_Lst)

        State_A_Energies_df, State_B_Energies_df = dE.dE_ParallelCalculationPrepare(State_A_df, State_B_df)

        dEs = dE.Run_dE_ParallelCalculation(State_A_Energies_df,State_B_Energies_df)

    else:
        

        State_A_RawEnergies_Lst, State_B_RawEnergies_Lst = Binary.ReadBinary(EnergyFiles_Lst)
        
        State_A_df = Binary.createDataFrames(State_A_RawEnergies_Lst)
        State_B_df = Binary.createDataFrames(State_B_RawEnergies_Lst)

        dEs = dE.dE_Calculation(State_A_df, State_B_df)
        dEs2=dE.dEs_matrix(State_A_df, State_B_df)
    return dEs, State_A_df, State_B_df,dEs2

