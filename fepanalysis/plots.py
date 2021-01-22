import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def Plot_Convergence(df):
    plt.plot(df['Number of Steps'],df['dG'])
    plt.title('Convergence Plot',fontsize=16)
    plt.xlabel("Number of Steps",fontsize=14)
    plt.ylabel("ΔG FEP (Kcal/mol)",fontsize=14)
    plt.savefig('Convergence.png',dpi=300)
    plt.close()


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

def Plot_PDF(State_A_df, State_B_df):

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
    plt.close()


def Plot_PDF_Matrix(State_A_df, State_B_df):
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
    plt.close()