import estimators
import pandas as pd
 
# def Convergence(df,Estimator,StepsChunk_Int,ReplicatiesCount_Int,EnergyOutputInterval_Int):
#                                     # the last and first steps are not included in the reading
#     Zwanzig_Final_Lst=[eval('estimators.Estimators.'+Estimator)(df,steps_limit)[1] for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,len(df)+1,StepsChunk_Int*ReplicatiesCount_Int)]
#     StepsChunk_Lst=[EnergyOutputInterval_Int*steps_limit/ReplicatiesCount_Int for steps_limit in range((StepsChunk_Int-2)*ReplicatiesCount_Int,len(df)+1,StepsChunk_Int*ReplicatiesCount_Int)]
#     Convergence_df=pd.DataFrame({'Number of Steps':StepsChunk_Lst, 'dG':Zwanzig_Final_Lst })
#     return Convergence_df
