#%%
from .inout import commands
from .inout import Q
from . import estimators
from . import tools
from . import plots
from pymbar import MBAR as MBAR_ 
from alchemlyb.estimators import MBAR
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

    dEs, State_A_df, State_B_df,dEs2 = Q.parser(args)
    DG_df = None

    if args.estimator =='Zwanzig_Estimator':

        # DG_df, Zwanzig_dG= estimators.Estimators.Zwanzig(dEs,None)
        # DG_df2, Zwanzig_dG= estimators.Estimators.Zwanzig_matrix_AI(dEs,None)
        #Q.dE.Get_dEs_dGs_AI(DG_df2,dEs2)
        # #estimators.Estimators.Zwanzig_matrix_AI(dEs,None)
        Zwanzig_dG= estimators.Estimators.Zwanzig_matrix_AI2(dEs2,None,args.dEs_dGs_file_name)
        # dU_dH_df=estimators.Estimators.Create_df_TI(State_A_df, State_B_df)
        # Zwanzig_dG = DG_df['dG_Average'].iloc[-1]
        # TI_dG = estimators.Estimators.TI(State_A_df, State_B_df,None)[1]
        # # print("ZW: ",Zwanzig_dG)
        # print("TI: ", TI_dG)
        # u_nk_df,states_dicts,State_A_Energies_df= estimators.Estimators.Create_df_BAR_MBAR(State_A_df, State_B_df)
        # BAR_dG= estimators.BAR().fit(u_nk_df)
        # print('BAR: ',BAR_dG.delta_f_.loc[0.00, 1.00])
        # # MBAR_df=estimators.Estimators.Create_df_MBAR(states_dicts,State_A_Energies_df,None)
        # # print('MBAR',MBAR_df)

        # MBAR_df= MBAR().fit(u_nk_df)
        # print('MBAR: ',MBAR_df.delta_f_.loc[0.00, 1.00])
        
    if args.convergence_analysis is not None:
        args.convergence_analysis=args.convergence_analysis[0].split(',')

        convergenc_Zwanzig = estimators.Estimators.Convergence(dEs,'BLANK',estimators.Estimators.Zwanzig, int(args.convergence_analysis[1]), int(args.convergence_analysis[2]), int(args.convergence_analysis[3]))
        
        convergenc_TI = estimators.Estimators.Convergence(State_A_df, State_B_df,estimators.Estimators.TI, int(args.convergence_analysis[1]), int(args.convergence_analysis[2]), int(args.convergence_analysis[3]))

        convergenc_BAR = estimators.Estimators.Convergence(states_dicts, State_A_Energies_df,estimators.Estimators.Create_df_dG_BAR, int(args.convergence_analysis[1]), int(args.convergence_analysis[2]), int(args.convergence_analysis[3]))
    
        convergenc_MBAR = estimators.Estimators.Convergence(states_dicts, State_A_Energies_df,estimators.Estimators.Create_df_dG_MBAR, int(args.convergence_analysis[1]), int(args.convergence_analysis[2]), int(args.convergence_analysis[3]))

        plots.plotting.Plot_Convergence(convergenc_Zwanzig,convergenc_TI,convergenc_BAR,convergenc_MBAR)
    if args.plot ==True:

        plots.plotting.Plot_Hysteresis(DG_df)
        plots.plotting.Plot_dG_by_Lambda(DG_df)
        plots.plotting.Plot_dEs(dEs)
        plots.plotting.Plot_PDF(State_A_df, State_B_df) 
        plots.plotting.Plot_PDF_Matrix(State_A_df, State_B_df)    
        plots.plotting.Plot_PDF2(dEs)
        plots.plotting.Plot_PDF_Matrix2(dEs)  

if __name__ == "__main__":

    main()
