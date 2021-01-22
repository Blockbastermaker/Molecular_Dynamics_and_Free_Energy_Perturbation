import numpy as np
import pandas as pd
import re

import logging
logger = logging.getLogger(__name__)


def Zwanzig_Estimator(dEs,steps):

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

