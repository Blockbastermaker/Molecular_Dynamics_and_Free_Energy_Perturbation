#%%

from .inout import commands
from .inout import Q
from . import estimators
from . import tools
from . import plots
import os
import logging
logger = logging.getLogger(__name__)

def generatedLogger(logFileName):

    FORMATTER = logging.Formatter("%(asctime)s — %(levelname)s — %(message)s", 
        datefmt='%Y-%m-%d %H:%M:%S')

    LOG_FILE = os.path.join(logFileName)
    file_handler = logging.FileHandler(LOG_FILE,'w')
    file_handler.setFormatter(FORMATTER)

    logger = logging.getLogger(__name__)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.ERROR)

    logger.addHandler(console_handler)

    logger.setLevel(logging.INFO)


def main():

    args = commands.parseArguments()

    generatedLogger('feptool.log')

    logger.info('Starting FEP analysis')

    dEs, State_A_df, State_B_df = Q.parser(args)
    DG_df = None

    if args.estimator =='Zwanzig_Estimator':

        DG_df = estimators.Estimators.Zwanzig(dEs,None)
        dU_dH_df=estimators.Estimators.Create_df_TI(State_A_df, State_B_df)
        Zwanzig_dG = DG_df['dG_Average'].iloc[-1]
        TI_dG = estimators.Estimators.TI(dU_dH_df)
        print("ZW: ",Zwanzig_dG )
        print("TI: ", TI_dG)
        
    if args.convergence_analysis is not None:
        print(args.convergence_analysis[0].split(','))
        args.convergence_analysis=args.convergence_analysis[0].split(',')
        convergenc_df = estimators.Estimators.Convergence(dEs,args.convergence_analysis[0], int(args.convergence_analysis[1]), int(args.convergence_analysis[2]), int(args.convergence_analysis[3]))
        print(convergenc_df)
        
        plots.Plot_Convergence(convergenc_df)

    if args.plot ==True:

        plots.plotting.Plot_Hysteresis(DG_df)
        plots.plotting.Plot_dG_by_Lambda(DG_df)
        plots.plotting.Plot_dEs(dEs)
        plots.plotting.Plot_PDF(State_A_df, State_B_df) 
        plots.plotting.Plot_PDF_Matrix(State_A_df, State_B_df)    

if __name__ == "__main__":

    main()


# %%
# from datetime import datetime
# start_time = datetime.now()
# dEs =  dE_Calculation(None)
# end_time = datetime.now()
# print('Duration: {}'.format(end_time - start_time))

# x=[i.split('\\')[1] for i in EnergyFiles_Lst]
# x.count([i.split('\\')[1] for i in EnergyFiles_Lst][0])


# %%
