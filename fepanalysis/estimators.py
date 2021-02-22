import numpy as np
import pandas as pd
import re
import glob
import logging
logger = logging.getLogger(__name__)

class Estimators():
            
    def Zwanzig(dEs,steps):

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

        return Zwanzig_df

    def Create_df_TI(State_A_df, State_B_df):
        
        dU_dH_df=(pd.DataFrame({"lambda":State_A_df["Lambda"],"fep":State_B_df["Q_sum"] - State_A_df["Q_sum"] })).sort_values('lambda')
        dU_dH_df.reset_index(drop=True,inplace=True)
        dU_dH_df.index.names = ['time']
        dU_dH_df.set_index(['lambda'], append=True,inplace=True)
        return dU_dH_df




    def TI(dHdl):
        
        
            """Thermodynamic integration (TI).
            
            Compute free energy differences between each state by integrating
            dHdl across lambda values.
            Parameters
            ----------
            dHdl : DataFrame 
            ----------
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
            return TI


    def Convergence(df,Estimator,StepsChunk_Int,ReplicatiesCount_Int,EnergyOutputInterval_Int):
                                        # the last and first steps are not included in the reading
        # if Estimator=='Zwanzig':
        #     Estimator=Estimators.Zwanzig
            
        Zwanzig_Final_Lst=[Estimator(df,steps_limit)[1] for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,len(df)+1,StepsChunk_Int*ReplicatiesCount_Int)]
        StepsChunk_Lst=[EnergyOutputInterval_Int*steps_limit/ReplicatiesCount_Int for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,len(df)+1,StepsChunk_Int*ReplicatiesCount_Int)]
        Convergence_df=pd.DataFrame({'Number of Steps':StepsChunk_Lst, 'dG':Zwanzig_Final_Lst })
        return Convergence_df
