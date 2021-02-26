import numpy as np
from numpy.core.numeric import _rollaxis_dispatcher
import pandas as pd
from pymbar import BAR as BAR_
from sklearn.base import BaseEstimator
import glob
import re
import itertools
import logging
logger = logging.getLogger(__name__)

class Estimators():
    """
    Return Estimated binding free energy (dG).
    
    Returns the dG between state A and state B using 3 differant Energy estimators
    
    Zwanzig, Thermodynamic Integration TI, or Bennett Acceptance Ratio (BAR).
    
    """
            
    def Zwanzig(dEs,steps):
        
        """
        Return the estimated binding free energy using Zwanzig estimator.
        
        Computes the binding free (dG) form molecular dynamics simulation 
        between state A and state B using Zwanzig estimator.
        
        Parameters
        ----------
        dEs : Pandas Dataframe
            contains the reduced potentail (dE) between the states.
        
        steps : interger 
            the number of the steps to be included in the calculation, set to "None" if all steps are needed.
        
        Returns
        ---------
        Zwanzig_df : Pandas Dataframe
            contains the binding free energy (dG) between the states.
            
        Examples
        --------
        >>> Zwanzig(dEs,None)
                
        >>> Zwanzig(dEs,1000)
        
        """
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
        logger.info('Final DG computed from Zwanzig estimator: ' +str(Zwanzig_Final_dG))

        return Zwanzig_df, Zwanzig_Final_dG    
    

    def Create_df_TI(State_A_df, State_B_df):
        """
        create the input dataframe needed for the Thermodynamic Integration (TI) function.
        
            Parameters
            ----------
            State_A_df : Pandas DataFrame for state A energies
            State_B_df : Pandas DataFrame for state B energies
            ----------
            Returns
            ----------
            dU_dH_df : Pandas DataFrame 
        
        """
        
        dU_dH_df=(pd.DataFrame({"lambda":State_A_df["Lambda"],"fep":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('lambda')
        dU_dH_df.reset_index(drop=True,inplace=True)
        dU_dH_df.index.names = ['time']
        dU_dH_df.set_index(['lambda'], append=True,inplace=True)
        return dU_dH_df




    def TI(State_A_df,State_B_df,steps):
        
        
        """
        Return the estimated binding free energy using Thermodynamic integration (TI) estimator.
        
        Compute free energy differences between each state by integrating
        dHdl across lambda values.
        Parameters
        ----------
        dHdl : Pandas DataFrame 
        ----------
        Returns
        ----------
        delta_f_ : DataFrame
            The estimated dimensionless free energy difference between each state.
        d_delta_f_ : DataFrame
            The estimated statistical uncertainty (one standard deviation) in 
            dimensionless free energy differences.
        states_ : list
            Lambda states for which free energy differences were obtained.
        
        TI : float
            The free energy difference between state 0 and state 1.
        """
        if steps != None:
            Energies_df=(pd.DataFrame({"lambda":State_A_df["Lambda"],"fep":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('lambda')
            Energies_df=pd.DataFrame.from_dict(dict(Energies_df.groupby('lambda',sort=False)['fep'].apply(list)),orient='index')
            Energies_df=Energies_df.transpose()
            Energies_df=Energies_df.iloc[:steps]

            dfl=pd.DataFrame(columns=['lambda','fep'])
            dU_dH_df=pd.DataFrame(columns=['lambda','fep'])
            for state in range (len(Energies_df.columns)):
                    dfl=pd.DataFrame(columns=['lambda','fep'])
                    dfl['fep']=Energies_df.iloc[:,state]
                    dfl['lambda']=Energies_df.columns.values[state]
                    dU_dH_df=dU_dH_df.append(dfl)
        else:
            dU_dH_df=(pd.DataFrame({"lambda":State_A_df["Lambda"],"fep":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('lambda')
        
        dU_dH_df.reset_index(drop=True,inplace=True)
        dU_dH_df.index.names = ['time']
        dU_dH_df.set_index(['lambda'], append=True,inplace=True)        
        # dU_dH_df=(pd.DataFrame({"lambda":State_A_df["Lambda"][:steps],"fep":State_B_df["Q_sum"][:steps] - State_A_df["Q_sum"][:steps] })).sort_values('lambda')
        # dU_dH_df.reset_index(drop=True,inplace=True)
        # dU_dH_df.index.names = ['time']
        # dU_dH_df.set_index(['lambda'], append=True,inplace=True)
        dHdl=dU_dH_df

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
        states_ = means.index.values.tolist()
        TI=( delta_f_.loc[0.00, 1.00])
        return delta_f_ , TI


    def Create_df_BAR_MBAR(State_A_df, State_B_df):
        
        """
        Create the input dataframe needed for the Bennett Acceptance Ratio (BAR) and multistate Bennett Acceptance Ratio (MBAR) estimators.
        
            Parameters
            ----------
            State_A_df : Pandas DataFrame for state A energies
            State_B_df : Pandas DataFrame for state B energies
            ----------
            Returns
            ----------
            u_nk_df : Pandas DataFrame 
        
        """
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
        States_dicts={i:[] for i in range(len(lambdas_list_A))}
        for i in range(len(State_A_Energies_df.columns)):
            State_A_Energies=State_A_Energies_df.iloc[:,[i]]
            State_A_Energies.columns=["0"]
            State_A_Lambda_float=State_A_Energies_df.columns[i]    
            
            State_B_Energies=State_B_Energies_df.iloc[:,[i]]
            State_B_Energies.columns=["0"]
            State_B_Lambda_float=State_B_Energies_df.columns[i]    
            E0=State_A_Energies*State_A_Lambda_float+State_B_Energies*State_B_Lambda_float
            for x in range(len(lambdas_list_A)):
                E1=State_A_Energies*lambdas_list_A[x]+State_B_Energies*lambdas_list_B[x]
                dE=E1-E0
                dE=dE.values.tolist()
                dE=list(itertools.chain(*dE))
                States_dicts[i].append(dE)
        for i in range(len(States_dicts)):
            States[i]=list(itertools.chain(*States_dicts[i]))
        u_nk_df=pd.DataFrame.from_dict(States)
        u_nk_df.columns=lambdas_list_A
        lambdas_df=lambdas_df*len(State_A_Energies_df)
        lambdas_df.sort()
        u_nk_df['time']=time*len(State_A_Energies_df.columns)
        u_nk_df['fep-lambda']=lambdas_df
        u_nk_df=u_nk_df.astype('float')
        u_nk_df.set_index(['time'] ,append=False,inplace=True)
        u_nk_df.set_index(['fep-lambda'], append=True,inplace=True)
        u_nk_df.columns= u_nk_df.columns.astype('float')
        u_nk_df.dropna(axis=0,inplace=True)
        return u_nk_df,States_dicts,State_A_Energies_df




    def Create_df_BAR_MBAR_2(States_dicts,State_A_Energies_df,steps):
        States_dicts_2={}
        lambdas_list_A=list(State_A_Energies_df.columns)
        time = [i for i in range(len(State_A_Energies_df))]
        lambdas_df=lambdas_list_A
        for x in States_dicts.keys():
            for i in range(len(States_dicts[x])):
                print(len(States_dicts[x]))
                States_dicts[x][i]=States_dicts[x][i][:steps]
                
        for i in range(len(States_dicts)):
            States_dicts_2[i]=list(itertools.chain(*States_dicts[i]))
        u_nk_df=pd.DataFrame.from_dict(States_dicts_2)
        u_nk_df.columns=lambdas_list_A
        lambdas_df=lambdas_df*len(State_A_Energies_df.iloc[:steps])
        lambdas_df.sort()
        u_nk_df['time']=time[:steps]*len(State_A_Energies_df.columns)
        u_nk_df['fep-lambda']=lambdas_df
        u_nk_df=u_nk_df.astype('float')
        u_nk_df.set_index(['time'] ,append=False,inplace=True)
        u_nk_df.set_index(['fep-lambda'], append=True,inplace=True)
        u_nk_df.columns= u_nk_df.columns.astype('float')
        u_nk_df.dropna(axis=0,inplace=True)
    
        BAR_df=BAR().fit(u_nk_df)
        BAR_dG = BAR_df.delta_f_.loc[0.00, 1.00]
        return BAR_dG

    def Convergence(df1,df2,Estimator,StepsChunk_Int,ReplicatiesCount_Int,EnergyOutputInterval_Int):
                                        # the last and first steps are not included in the reading
        """
        Convergence analysis
        
        Retrun a dateframe contains computed free energy dG at a differant steps intervals using 3 differant Energy estimators
    
        Zwanzig, Thermodynamic Integration TI, or Bennett Acceptance Ratio (BAR).
    
            Parameters
            ----------
            df : Pandas DataFrame 
                Contains the dEs between the states
                
            Estimator : funcation 
                    The Free energy estimating method (Zwanzig or TI or BAR)
            
            StepsChunk_Int: integer
                     The Number of Steps(fs) to be used.
                     
            ReplicatiesCount_Int: integer
                    The Number of used replicates.
            
            EnergyOutputInterval_Int: integer
                    The interval which the molecular dynamics simulation softwear
                    is writing the energies at.
            ----------
            Returns
            ----------
            Convergence_df : Pandas DataFrame 
                            Contains the computed dG at each interval.
            
            Examples
            --------
        >>> Convergence(dEs,Zwanzig,1000,1,10)
        
        >>> Convergence(dEs,TI,10000,3,10)        
        """
        if isinstance(df1, pd.DataFrame) and isinstance(df2, pd.DataFrame) :
            dGs_Lst=[Estimator(df1,df2,steps_limit)[1] for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,int((len(df1)/len(df1['Lambda'].unique())))+1,StepsChunk_Int*ReplicatiesCount_Int)]
            StepsChunk_Lst=[EnergyOutputInterval_Int*steps_limit/ReplicatiesCount_Int for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,int((len(df1)/len(df1['Lambda'].unique())))+1,StepsChunk_Int*ReplicatiesCount_Int)]

        elif isinstance(df2, pd.DataFrame) and not isinstance(df1, pd.DataFrame):
            #Estimator(df1,df2,steps_limit)
            dGs_Lst=[Estimator(df1,df2,int(steps_limit)) for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,len(df2)+1,StepsChunk_Int*ReplicatiesCount_Int)]
            StepsChunk_Lst=[EnergyOutputInterval_Int*steps_limit/ReplicatiesCount_Int for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,len(df2)+1,StepsChunk_Int*ReplicatiesCount_Int)]
        else:
            dGs_Lst=[Estimator(df1,steps_limit)[1] for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,len(df1)+1,StepsChunk_Int*ReplicatiesCount_Int)]
            StepsChunk_Lst=[EnergyOutputInterval_Int*steps_limit/ReplicatiesCount_Int for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,len(df1)+1,StepsChunk_Int*ReplicatiesCount_Int)]


            
        #StepsChunk_Lst=[EnergyOutputInterval_Int*steps_limit/ReplicatiesCount_Int for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,len(df1)+1,StepsChunk_Int*ReplicatiesCount_Int)]
        Convergence_df=pd.DataFrame({'Number of Steps':StepsChunk_Lst, 'dG':dGs_Lst })
        return Convergence_df

    def Zwanzig_matrix_AI(dEs,steps):
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
        return Zwanzig_df ,Zwanzig_Final_dG


#from pymbar import BAR as BAR_

class BAR(BaseEstimator):
    """Bennett acceptance ratio (BAR).
    Parameters
    ----------
    maximum_iterations : int, optional
        Set to limit the maximum number of iterations performed.
    relative_tolerance : float, optional
        Set to determine the relative tolerance convergence criteria.
    method : str, optional, default='false-position'
        choice of method to solve BAR nonlinear equations,
        one of 'self-consistent-iteration' or 'false-position' (default: 'false-position')
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

    def __init__(self, maximum_iterations=10000, relative_tolerance=1.0e-7, method='false-position', verbose=False):

        self.maximum_iterations = maximum_iterations
        self.relative_tolerance = relative_tolerance
        self.method = method
        self.verbose = verbose

        # handle for pymbar.BAR object
        self._bar = None

    def fit(self, u_nk):
        """
        Compute overlap matrix of reduced potentials using
        Bennett acceptance ratio.
        Parameters
        ----------
        u_nk : DataFrame
            u_nk[n,k] is the reduced potential energy of uncorrelated
            configuration n evaluated at state k.
        """
        # sort by state so that rows from same state are in contiguous blocks
        u_nk = u_nk.sort_index(level=u_nk.index.names[1:])

        # get a list of the lambda states
        self.states_ = u_nk.columns.values.tolist()

        # group u_nk by lambda states
        groups = u_nk.groupby(level=u_nk.index.names[1:])
        N_k = [(len(groups.get_group(i)) if i in groups.groups else 0) for i in u_nk.columns]

        # Now get free energy differences and their uncertainties for each step
        deltas = np.array([])
        d_deltas = np.array([])
        for k in range(len(N_k) - 1):
            # get us from lambda step k
            uk = groups.get_group(self.states_[k])
            # get w_F
            w_f = uk.iloc[:, k+1] - uk.iloc[:, k]

            # get us from lambda step k+1
            uk1 = groups.get_group(self.states_[k+1])
            # get w_R
            w_r = uk1.iloc[:, k] - uk1.iloc[:, k+1]

            # now determine df and ddf using pymbar.BAR
            df, ddf = BAR_(w_f, w_r,
                             method=self.method,
                             maximum_iterations=self.maximum_iterations,
                             relative_tolerance=self.relative_tolerance,
                             verbose=self.verbose)

            deltas = np.append(deltas, df)
            d_deltas = np.append(d_deltas, ddf**2)

        # build matrix of deltas between each state
        adelta = np.zeros((len(deltas) + 1, len(deltas) + 1))
        ad_delta = np.zeros_like(adelta)

        for j in range(len(deltas)):
            out = []
            dout = []
            for i in range(len(deltas) - j):
                out.append(deltas[i:i + j + 1].sum())

                # See https://github.com/alchemistry/alchemlyb/pull/60#issuecomment-430720742
                # Error estimate generated by BAR ARE correlated

                # Use the BAR uncertainties between two neighbour states
                if j == 0:
                    dout.append(d_deltas[i:i + j + 1].sum())
                # Other uncertainties are unknown at this point
                else:
                    dout.append(np.nan)

            adelta += np.diagflat(np.array(out), k=j + 1)
            ad_delta += np.diagflat(np.array(dout), k=j + 1)

        # yield standard delta_f_ free energies between each state
        self.delta_f_ = pd.DataFrame(adelta - adelta.T,
                                     columns=self.states_,
                                     index=self.states_)

        # yield standard deviation d_delta_f_ between each state
        self.d_delta_f_ = pd.DataFrame(np.sqrt(ad_delta + ad_delta.T),
                                       columns=self.states_,
                                       index=self.states_)

        return self


