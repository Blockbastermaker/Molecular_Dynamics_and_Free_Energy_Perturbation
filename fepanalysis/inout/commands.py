import argparse

def parseArguments():

    parser = argparse.ArgumentParser(description="MD/FEP Analysis")

    parser.add_argument("-f","--energy_files_prefix", help = "Energy Files Prefix, ex: FEP1, FEP2")
    parser.add_argument("-a","--all_replicaties",default = False, action="store_true", help = "Analyze all replicaties in the current directory")
    parser.add_argument("-p","--run_in_parallel",action="store_true", help = "Run in parallel")
    parser.add_argument("-e","--estimator",nargs='+', default='Zwanzig_Estimator',help = "Energy Estimator")
    parser.add_argument("-c","--convergence_analysis",nargs='+', help = "Convergence Analysis: Estimator, by Number of Steps(fs), Number of used Replicaties")
    parser.add_argument("-t","--plot", default = False, action="store_true", help = "Plot and Save")
    parser.add_argument("-n","--dEs_dGs_file_name", help = "dEs_dGs File name, ex: Adrenaline_FEP1")
    args = parser.parse_args()

    return args