import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import chi2_contingency
import scipy.stats as stats
from math import log10, log2
import matplotlib_venn as mpv

def ICD10_code_to_chapter(let):
    if let == 'nan':
        return 'NaN';
    elif let[0] == 'A' or let[0] == 'B':
        return 'A00–B99';
    elif let[0] == 'C' or (let[0] == 'D' and int(let[1])>=0 and int(let[1])<5):
        return 'C00–D48';
    elif let[0] == 'D' and int(let[1])>=5 and int(let[1])<9:
        return 'D50–D89';
    elif let[0] == 'E':
        return 'E00–E90';
    elif let[0] == 'H' and int(let[1])>=0 and int(let[1])<6:
        return 'H00–H59';
    elif let[0] == 'H' and int(let[1])>=6 and int(let[1])<=9:
        return 'H60–H95';
    elif let[0] == 'K':
        return 'K00–K93';
    elif let[0] == 'P':
        return 'P00–P96';
    elif let[0] == 'S' or let[0] == 'T':
        return 'S00–T98';
    elif let[0] in ['V','W','X','Y']:
        return 'V01–Y98';
    elif let[0] in ['F', 'G','I', 'J', 'L', 'M', 'N', 'O','Q','R','Z','U']:
        return '{}00–{}99'.format(let[0], let[0]);
    else:
        return let;
    
def ICDchapter_to_name(chp):
    if chp == 'nan': return 'NaN';
    elif chp == 'A00–B99': return 'Certain infectious and parasitic diseases';
    elif chp == 'C00–D48': return 'Neoplasms';
    elif chp == 'D50–D89': return 'Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism';
    elif chp == 'E00–E90': return 'Endocrine, nutritional and metabolic diseases';
    elif chp == 'F00–F99': return 'Mental and behavioural disorders';
    elif chp == 'G00–G99': return 'Diseases of the nervous system';
    elif chp == 'H00–H59': return 'Diseases of the eye and adnexa';
    elif chp == 'H60–H95': return 'Diseases of the ear and mastoid process';
    elif chp == 'I00–I99': return 'Diseases of the circulatory system';
    elif chp == 'J00–J99': return 'Diseases of the respiratory system';
    elif chp == 'K00–K93': return 'Diseases of the digestive system';
    elif chp == 'L00–L99': return 'Diseases of the skin and subcutaneous tissue';
    elif chp == 'M00–M99': return 'Diseases of the musculoskeletal system and connective tissue';
    elif chp == 'N00–N99': return 'Diseases of the genitourinary system';
    elif chp == 'O00–O99': return 'Pregnancy, childbirth and the puerperium';
    elif chp == 'P00–P96': return 'Certain conditions originating in the perinatal period';
    elif chp == 'Q00–Q99': return 'Congenital malformations, deformations and chromosomal abnormalities';
    elif chp == 'R00–R99': return 'Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified';
    elif chp == 'S00–T98': return 'Injury, poisoning and certain other consequences of external causes';
    elif chp == 'V01–Y98': return 'External causes of morbidity and mortality';
    elif chp == 'Z00–Z99': return 'Factors influencing health status and contact with health services';
    elif chp == 'U00–U99': return 'Codes for special purposes';
    else: return ' ';

    
    import matplotlib_venn as mpv

def countPtsDiagnosis_Dict(csvfileordf, totalpts, diagkeys = diagkeys):
    '''
        diagnosisdict = countPtsWithWithoutDiagnosis(csvfile, totalpts)
        input: csvfileordf - string that leads to the csv file with diagnosis, or dataframe of diagnosis.
                   columns of csv/df should include: PatientID, ICD10_code, diagkeys
                totalpts - total number of patients in data
                diagkeys - diagnostic categories. default: ['DiagnosisName','Level3_Category','Level2_Category']
        outputs: diagnosisdict - diagnosis dictionary of dataframes. 
                  diagnosisdict[x] will give you the dataframe for the a diagkey category.
                  each dataframe includes columns: diagkey_category, Count, Count_r
                  where Count is the number of patient with a diagnosis, and Count_r = totalpts - Count
    '''
    if type(csvfileordf)==str: # if csv file...
        ptDiag = pd.read_csv(csvfileordf) # Read in file
        ptDiag['ValueL'] = ptDiag['Value'].apply(lambda x: ICD10_code_to_chapter(str(x)[0:3])) # ICD 10 Category
        if not(ptDiag['PatientID'].unique().shape[0] == totalpts):
            raise Exception('Patient Cohort Unique ID number doesn''t match up to Pt Number')
    elif type(csvfileordf)==pd.DataFrame: # if dataframe...
        ptDiag = csvfileordf
    else:
        raise Exception(str(type(csvfileordf)) + 'not supported. Please give csv file name as string, or dataframe.')
       
    # loop through diagkeys and save counts 
    ptDiagCount = dict()
    for n in diagkeys:
        diagtemp = ptDiag[['PatientID',n]].drop_duplicates() # drop duplicate diagnosis for each patient
        ptDiagCount[n]= pd.DataFrame(diagtemp[n].value_counts()).reset_index()
        ptDiagCount[n].columns = [n,'Count']
        ptDiagCount[n]['Count_r'] = totalpts - ptDiagCount[n]['Count']
        
    return ptDiagCount

def vennDiagramFromDiagnosisDict(ptdict1, ptdict2, diagkeys = diagkeys, outfile = None,
                                 set_labels = ('Group 1','Group 2'), title = 'Overlap between two groups',
                                 figsize = (10,4)):
    ''' 
        vennDiagramFromDiagnosisDict(ptdict1, ptdict2, outfile, diagkeys = diagkeys, outfile = None,
                                 set_labels = ('Group 1','Group 2'), title = 'Overlap between two groups',
                                 figsize = (10,4))
        input: ptdict1 and ptdict2 are dictionary of dataframes. 
                Both dictionary should have entries to every category in diagkeys to generate the venn diagram.
                outfile  - file for output figure
                set_labels - labels for the two groups
                title - title of figure
                figsize - figure size
        outputs: none. Creates a venn diagram. Saves it if outfile is given.
    '''
    fig = plt.figure(figsize = figsize)
    
    for i, n in enumerate(diagkeys):
        plt.subplot(1,len(diagkeys),i+1)
        pset1 = set(ptdict1[n][n])
        pset2 = set(ptdict2[n][n])
        mpv.venn2([pset1, pset2], set_labels = set_labels) # Use matplot_venn package
        plt.title(n + ' Overlap')
    fig.suptitle(title, fontsize = 16);
    if outfile: plt.savefig(outfile, bbox_inches = 'tight');
    else: plt.show();

def sigTestCounts(allcounts, n = None, verbose = False, diag = False): 
    ''' 
    combined = sigTestCounts(allcounts, n = None, verbose = False, diag = False)
    Inputs:
        allcounts - dataframe or dictionary.
            Dataframe is of format index with feature name, and 4 columns 
            that make up the contingency table of interest in the order of
            case positive, case negative, control positive, control negative. 
            Each row will be reshaped into a 2x2 contingency table. 
            If dictionary of dataframes, will extract dataframe with key n
        n - default: None. 
            If allcounts is dictionary, n is the key to extract the dataframe
            of interest.
        verbose - default: False
        diag - default: False. If true, will append ICD10 category to the output
            dataframes.
    -----
    Outputs:
        combined - dataframe with fisher or chi square stats and odds ratios 
        appended. '''
    # First, for fischer choose any row with less than 5 patients in a category
    if type(allcounts) is dict:
        allcounts = allcounts[n]
    print(n,': Amount: ', allcounts.shape[0])
    
    temp_less5 = allcounts[allcounts.min(axis=1)<5] # take all with counts less than 5
    fisher1 = pd.DataFrame()
    if temp_less5.shape[0]>0:
        print('\t Fisher Exact for <5 pts in a category, num:', temp_less5.shape[0])
        fisher =  temp_less5 \
            .apply(lambda x: stats.fisher_exact(np.array(x).reshape(2,2)), axis = 1) \
            .apply(pd.Series)
        fisher.columns = ['OddsRatio', 'pvalue']
        if verbose: print('\t\t fisher:',fisher.shape)

        maxratio = fisher['OddsRatio'][fisher['OddsRatio']<np.inf].max();
        minratio = fisher['OddsRatio'][fisher['OddsRatio']>0].min();
        fisher = fisher.replace(np.inf, maxratio+1) 
        fisher['log2_oddsratio'] = fisher['OddsRatio']\
            .apply(lambda x: log2(minratio/2) if (x==0) else log2(x))

        minpvalue = fisher['pvalue'][fisher['pvalue']>0].min();
        fisher['pvalue']=fisher['pvalue'].replace(0,minpvalue/2)
        fisher['-log_pvalue']=fisher['pvalue'].apply(lambda x: -log10(x))

        fisher1 = fisher.merge(temp_less5, how = 'right', left_index=True, right_index = True)
        if verbose: print('\t\t fisher1',fisher1.shape)

    # now take the rest of the patients
    temp_more5 = allcounts[allcounts.min(axis=1)>=5]
    print('\t Chi square for >=5 pts in a category, num:', temp_more5.shape[0])

    fisher =  temp_more5 \
        .apply(lambda x: stats.fisher_exact(np.array(x).reshape(2,2)), axis = 1) \
        .apply(pd.Series)
    fisher.columns = ['OddsRatio', 'fpvalue']
    
    maxratio = fisher['OddsRatio'][fisher['OddsRatio']<np.inf].max();
    minratio = fisher['OddsRatio'][fisher['OddsRatio']>0].min();
    fisher = fisher.replace(np.inf, maxratio+1) 
    fisher['log2_oddsratio'] = fisher['OddsRatio']\
        .apply(lambda x: log2(minratio) if (x==0) else log2(x))
    minpvalue = fisher['fpvalue'][fisher['fpvalue']>0].min();
    fisher['fpvalue']=fisher['fpvalue'].replace(0,minpvalue/2)
    fisher['-log_fpvalue']=fisher['fpvalue'].apply(lambda x: -log10(x))
    if verbose: print('\t\t fisher',fisher.shape)

    chisquare = temp_more5.apply(lambda x: \
                           chi2_contingency(np.array(x).reshape(2,2)), axis=1) \
                            .apply(pd.Series)
    chisquare.columns = ['chistat','pvalue','dof','expected']
    chisquare = chisquare.merge(temp_more5, how = 'right',left_index=True, right_index = True)
    minpvalue = chisquare['pvalue'][chisquare['pvalue']>0].min();
    chisquare['pvalue']=chisquare['pvalue'].replace(0,minpvalue/2)
    chisquare['-log_pvalue']=chisquare['pvalue'].apply(lambda x: -log10(x))
    if verbose: print('\t\t chisquare:',chisquare.shape)

    combined = chisquare.merge(fisher, left_index=True, right_index = True, how = 'left')
    combined = combined.append(fisher1)
    if verbose: print('\t\t combined 1:', combined.shape)
    
    if diag: 
        temp = alzdiag[[n,'ValueL']].append(condiag[[n, 'ValueL']]).drop_duplicates() # create mapping
        temp = temp[temp['ValueL'] != 'NaN'].groupby(n)['ValueL'].apply(list)
        combined = combined.merge(temp, how = 'left', left_index=True, right_index=True, suffixes=(False, False))
        if verbose: print('\t\t combined 2:', combined.shape)

    print('\t Final num: ', combined.shape[0])
    
    return combined;

def sigTestCountsDict(allcountsdict, dictkeys, verbose = False, diag = False):
    ''' combined = sigTestCountsDict(allcountsdict, verbose = False, diag = False)
    Inputs:
    allcountsdict - dictionary of dataframes.
        Each dataframe is of format index with feature name, and 4 columns 
        that make up the contingency table of interest in the order of
        case positive, case negative, control positive, control negative. 
        Each row will be reshaped into a 2x2 contingency table. 
    dictkeys - iterable object with the keys of allcountsdict.
    verbose - default: False
    diag - default: False. If true, will append ICD10 category to the output
        dataframes.
    ----
    Outputs: 
    combined - dictionary of dataframes with fisher or chi square 
    stats and odds ratios appended to each dataframe. '''
    
    combineddict = dict()
    for n in dictkeys:
        print('Significance testing on ', n)
        combineddict[n] = sigTestCounts(allcountsdict, n,verbose= verbose,diag = diag)
    return combineddict

def getXCounts(in_df, totalpts, col):
    ''' counts = getXCounts(in_df, totalpts, col)
        inputs: in_df: input dictionary. should have columns 'PatientDurableKey' and 'col'
                totalpts: total patients for counting patients without 
                col: column to count entries of
        outputs: counts = dataframe with counts and counts_r
    '''
    temp = in_df[['PatientDurableKey',col]].drop_duplicates()
    temp = temp.groupby(col)['PatientDurableKey'].nunique().reset_index()
    temp.columns = [col, 'Count']
    temp['Count_r'] = totalpts - temp['Count']
    return temp

def getXCountsStratify(stratvarname, stratvars, colCount, in_df = None, in_dict=None, numptsvar = None, 
                               equalize_num = False, random_state = 40):
    ''' 
    stratXcount = getXCountsStratify(stratvarname, stratvars, colCount, in_df = None, in_dict=None, numptsvar = None, 
                               equalize_num = False, random_state = 40)
    Inputs: stratvarname: string. name of variable, like 'Sex'
            stratvars: list of stratvarname group, like ['Male','Female']
            colCount: column to count entried over
            in_df or in_dict = either dataframe of unstratified data, or dictionary of dataframes pre-stratified
            numptsvar = dictionary, # patients in each category. Leave empty if passed in dataframe
            equalize_num = output counts have equal pts in each category
            random_state = random state number.
    Outputs: stratXcount: dictionary of the counts in each stratified variable.
    '''
    Xdictvar = dict()
    if in_dict: #check if dictionary of dataframes, or a full dataframe
        Xdictvar = in_dict;
        if numptsvar is None:
            raise Exception('variable numptsvar is empty, pass in dictionary with number of patients in each category.')
    elif in_df:
        numptsvar = dict()
        for var in stratvars:
            Xdictvar[var] = in_df[in_df[stratvarname] == var]
            numptsvar[var] = Xdictvar[var][['PatientDurableKey',stratvarname]].drop_duplicates().shape[0]
    else:
        raise Exception('did not pass in full dataframe or dictionary of stratified dataframes.')
     
    # subsample patients
    kNumPatientsSampled = min(numptsvar.values())
    varmin = list(numptsvar.keys())[np.array(list(numptsvar.values())).argmin()]
    if equalize_num:
        print('Number sampled:', kNumPatientsSampled)
        for var in stratvars:
            if var != varmin:
                subsampledPatientKeys = Xdictvar[var]['PatientDurableKey'].drop_duplicates()\
                    .sample(kNumPatientsSampled, random_state = randomstate)
    
    stratXcount = dict()
    for var in stratvars:
        stratXcount[var] = dict()
        if numptsvar[var]>kNumPatientsSampled and equalize_num:
            subsampledPatientKeys = Xdictvar[var]['PatientDurableKey'].drop_duplicates()\
                        .sample(kNumPatientsSampled, random_state = random_state)
            Xdictvar_s = Xdictvar[var][Xdictvar[var]['PatientDurableKey'].isin(subsampledPatientKeys)]
            stratXcount[var] = getXCounts(Xdictvar_s, kNumPatientsSampled, colCount)
        else: 
            stratXcount[var] = getXCounts(Xdictvar[var], numptsvar[var], colCount)
    return stratXcount;

# from bioinfokit
class general:
    def __init__(self):
        pass

    rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                   '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                   '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                   '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')

    def get_figure(show, r, figtype, fig_name):
        if show:
            plt.show()
        else:
            plt.savefig(fig_name+'.'+figtype, format=figtype, bbox_inches='tight', dpi=r)

    def axis_labels(x, y, axlabelfontsize=None, axlabelfontname=None):
        plt.xlabel(x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.ylabel(y, fontsize=axlabelfontsize, fontname=axlabelfontname)
        # plt.xticks(fontsize=9, fontname="sans-serif")
        # plt.yticks(fontsize=9, fontname="sans-serif")

    def axis_ticks(xlm=None, ylm=None, axtickfontsize=None, axtickfontname=None, ar=None):
        if xlm:
            plt.xlim(left=xlm[0], right=xlm[1])
            plt.xticks(np.arange(xlm[0], xlm[1], xlm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.xticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.yticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

    def check_for_nonnumeric(pd_series=None):
        if pd.to_numeric(pd_series, errors='coerce').isna().sum() == 0:
            return 0
        else:
            return 1
        
# Manhattan plot modified bioinfokit package source code
# https://github.com/reneshbedre/bioinfokit
def mhat(df="dataframe", chr=None, pv=None, color=None, dim=(6,4), r=300, ar=90, gwas_sign_line=False,
             gwasp=5E-8, dotsize=8, markeridcol=None, markernames=None, gfont=8, valpha=1, show=False, figtype='png',
             axxlabel=None, axylabel=None, axlabelfontsize=9, axlabelfontname="Arial", axtickfontsize=9, figtitle = 'manhattan plot',
             axtickfontname="Arial", ylm=None, gstyle=1, yskip = 1, plotlabelrotation = 0, figname='manhattan', 
             invert = False, fig = None, ax = None, xtickname=False):

        _x, _y = 'Chromosomes', r'$ -log_{10}(P)$'
        rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                       '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                       '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                       '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')
        
        # minus log10 of P-value
        if invert:
            df['tpval'] = np.log10(df[pv])
        else: 
            df['tpval'] = -np.log10(df[pv])
        df = df.sort_values(chr)
        
        # add indices
        df['ind'] = range(len(df))
        df_group = df.groupby(chr)
        if color is not None and len(color) == 2:
            color_1 = int(df[chr].nunique() / 2) * [color[0]]
            color_2 = int(df[chr].nunique() / 2) * [color[1]]
            if df[chr].nunique() % 2 == 0:
                color_list = list(reduce(lambda x, y: x+y, zip(color_1, color_2)))
            elif df[chr].nunique() % 2 == 1:
                color_list = list(reduce(lambda x, y: x+y, zip(color_1, color_2)))
                color_list.append(color[0])
        elif color is not None and len(color) == df[chr].nunique():
            color_list = color
        elif color is None:
            # select colors randomly from the list based in number of chr
            # color_list = sample(rand_colors, df[chr].nunique())
            color_list = rand_colors[:df[chr].nunique()]
        else:
            print("Error: in color argument")
            sys.exit(1)

        # label x-axis
        xlabels = []
        xticks = []
        if fig is None:
            fig, ax = plt.subplots(figsize=dim)
        i = 0
        for label, df1 in df.groupby(chr):
            df1.plot(kind='scatter', x='ind', y='tpval', color=color_list[i], s=dotsize, alpha=valpha, ax=ax)
            df1_max_ind = df1['ind'].iloc[-1]
            df1_min_ind = df1['ind'].iloc[0]
            xlabels.append(label)
            xticks.append((df1_max_ind - (df1_max_ind - df1_min_ind) / 2))
            i += 1

        # add GWAS significant line
        if gwas_sign_line is True:
            ax.axhline(y=-np.log10(gwasp), linestyle='--', color='#7d7d7d', linewidth=1)
        if markernames is not None:
            geneplot_mhat(df, markeridcol, chr, pv, gwasp, markernames, gfont, gstyle, ax=ax, plotlabelrotation =  plotlabelrotation)
        ax.margins(x=0)
        ax.margins(y=0)
        ax.set_xticks(xticks)
        ax.set_ylim([0, max(df['tpval'] + 1)+10])
        if ylm:
            ylm = np.arange(ylm[0], ylm[1], ylm[2])
        else:
            ylm = np.arange(0, max(df['tpval']+1), yskip)
        ax.set_yticks(ylm)
        if xtickname:
            ax.set_xticklabels(map(ICDchapter_to_name,xlabels), rotation = ar, va = 'top',ha='right')
        else: 
            ax.set_xticklabels(xlabels, rotation=ar)
        ax.set_yticklabels(ylm, fontsize=axtickfontsize, fontname=axtickfontname, rotation=ar)
        
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        ax.set_xlabel(_x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        ax.set_ylabel(_y, fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.title(figtitle)
        return fig, ax

def geneplot_mhat_logp(df, markeridcol, chr, tpval, gwasp, markernames, gfont, gstyle, ax, plotlabelrotation, vertalign = 'bottom'):
        loggwasp = np.log10(gwasp)
        if markeridcol is not None:
            if markernames is not None and markernames is True:
                for i in df[markeridcol].unique():
                    if abs(df.loc[df[markeridcol] == i, tpval].iloc[0]) >= abs(loggwasp):
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0],
                                    str(i), fontsize=gfont, rotation = plotlabelrotation, va = vertalign)
                        elif gstyle == 2:
                            plt.annotate(i, xy=(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
            elif markernames is not None and isinstance(markernames, (tuple, list)):
                for i in df[markeridcol].unique():
                    if i in markernames:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0],
                                str(i), fontsize=gfont, rotation = plotlabelrotation)
                        elif gstyle == 2:
                            plt.annotate(i, xy=(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
            elif markernames is not None and isinstance(markernames, dict):
                for i in df[markeridcol].unique():
                    if i in markernames:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0],
                                 markernames[i], fontsize=gfont, rotation = plotlabelrotation)
                        elif gstyle == 2:
                            plt.annotate(markernames[i], xy=(
                            df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, tpval].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
        else:
            raise Exception("provide 'markeridcol' parameter")
            
# miami plot
def miami(df="dataframe", chromo=None, logp1=None, logp2=None, color=None, dim=(10,10), r=300, ar=90, gwas_sign_line=False,
             gwasp=5E-8, dotsize=8, markeridcol=None, markernames=None, gfont=8, valpha=1, show=False, figtype='png',
             axxlabel=None, axylabel=None, axlabelfontsize=9, axlabelfontname="Arial", axtickfontsize=9, figtitle = 'miami plot',
             label1='firstgroup', label2 = 'secondgroup',
             axtickfontname="Arial", ylm=None, gstyle=1, yskip = 1, plotlabelrotation = 0, figname='miami', invert = False, fig = None, ax = None):
        
        _x, _y = 'Chromosomes', r'$ -log_{10}(P)$'

        df['tpval'] = df[logp1]
        df['tpval2'] = -df[logp2]
        df = df.sort_values(chromo)

        df['ind'] = range(len(df))
        df_group = df.groupby(chromo)

        rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                               '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                               '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                               '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')
        color_list = rand_colors[:df[chromo].nunique()]

        xlabels = []
        xticks = []
        
        if fig is None:
            fig, (ax0, ax) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1,20]}, figsize = dim)
            ax0.axis('off')
            fig.tight_layout()

        i = 0
        for label, df1 in df.groupby(chromo):
            df1.plot(kind='scatter', x='ind', y='tpval', color=color_list[i], s=dotsize, alpha=valpha, ax=ax)
            df1.plot(kind='scatter', x='ind', y='tpval2', color=color_list[i], s=dotsize, alpha=valpha, ax=ax)
            df1_max_ind = df1['ind'].iloc[-1]
            df1_min_ind = df1['ind'].iloc[0]
            xlabels.append(label)
            xticks.append((df1_max_ind - (df1_max_ind - df1_min_ind) / 2))
            i += 1

        ax.axhline(y=0, color='#7d7d7d', linewidth=.5, zorder = 0)

        # add GWAS significant line
        if gwas_sign_line is True:
            ax.axhline(y=np.log10(gwasp), linestyle='--', color='#7d7d7d', linewidth=1)
            ax.axhline(y=-np.log10(gwasp), linestyle='--', color='#7d7d7d', linewidth=1)
        if markernames is not None:
            marker.geneplot_mhat_logp(df, markeridcol, chromo, 'tpval', gwasp, markernames, gfont, gstyle, ax=ax, plotlabelrotation = plotlabelrotation)
            marker.geneplot_mhat_logp(df, markeridcol, chromo, 'tpval2', gwasp, markernames, gfont, gstyle, ax=ax, plotlabelrotation = -plotlabelrotation, vertalign = 'top')

        ax.margins(x=0)
        ax.margins(y=0)
        ax.set_xticks(xticks);
        (ymin, ymax) = (min(df['tpval2']-1)-10, max(df['tpval']+1)+10)
        ax.set_ylim([ymin, ymax])
        ax0.set_ylim([ymin, ymax])
        
        ax0.text(0,ymin/2,label2, fontsize=axlabelfontsize, fontname=axlabelfontname, rotation = 90, va = 'center')
        ax0.text(0,ymax/2,label1, fontsize=axlabelfontsize, fontname=axlabelfontname, rotation = 90, va = 'center')
        
        if ylm:
            ylm = np.arange(ylm[0], ylm[1], ylm[2])
        else:
            ylm = np.concatenate((np.arange(0,min(df['tpval2']-10),-yskip), np.arange(0, max(df['tpval']+10), yskip)))
            ax.set_yticks(ylm);
        ax.set_xticklabels(xlabels, rotation=ar)
        ax.set_yticklabels(ylm.astype(int), fontsize=axtickfontsize, fontname=axtickfontname);
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        ax.set_xlabel(_x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        ax.get_yaxis().get_label().set_visible(False)
        
        ax0.text(.5,0,_y, fontsize=axlabelfontsize, fontname=axlabelfontname, rotation = 90, va = 'center')
        
        plt.title(figtitle)
        general.get_figure(show, r, figtype, figname)
        return fig, ax