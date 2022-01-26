import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
from random import sample
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import seaborn as sns
#from adjustText import adjust_text

# Code is adapted from bioinfokit package by Renesh Bendre
# https://github.com/reneshbedre/bioinfokit

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
        #plt.close()

    def axis_labels(x, y, axlabelfontsize=None, axlabelfontname=None):
        plt.xlabel(x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.ylabel(y, fontsize=axlabelfontsize, fontname=axlabelfontname)

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

        
class marker:
    def __init__(self):
        pass

    def geneplot_mhat(df, markeridcol, chr, pv, gwasp, markernames, gfont, gstyle, ax, plotlabelrotation):
        if markeridcol is not None:
            if markernames is not None and markernames is True:
                for i in df[markeridcol].unique():
                    if df.loc[df[markeridcol] == i, pv].iloc[0] <= gwasp:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0], str(i), fontsize=gfont, rotation = plotlabelrotation)
                        elif gstyle == 2:
                            plt.annotate(i, xy=(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
            elif markernames is not None and isinstance(markernames, (tuple, list)):
                for i in df[markeridcol].unique():
                    if i in markernames:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0], str(i), fontsize=gfont, rotation = plotlabelrotation)
                        elif gstyle == 2:
                            plt.annotate(i, xy=(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
            elif markernames is not None and isinstance(markernames, dict):
                for i in df[markeridcol].unique():
                    if i in markernames:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0], markernames[i], fontsize=gfont, rotation = plotlabelrotation)
                        elif gstyle == 2:
                            plt.annotate(markernames[i], xy=(
                            df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
        else:
            raise Exception("provide 'markeridcol' parameter")
            
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
    
    def mhat(df="dataframe", chr=None, pv=None, color=None, dim=(6,4), r=300, ar=90, gwas_sign_line=False,
             gwasp=5E-08, dotsize=8, markeridcol=None, markernames=None, gfont=8, valpha=1, show=False, figtype='png',
             axxlabel=None, axylabel=None, axlabelfontsize=9, axlabelfontname="Arial", axtickfontsize=9, figtitle = 'manhattan plot',
             axtickfontname="Arial", ylm=None, gstyle=1, yskip = 1, plotlabelrotation = 0, figname='manhattan', 
             invert = False, fig = None, ax = None, xtickname=False, rand_colors = None):

        _x, _y = 'Chromosomes', r'$ -log_{10}(P)$'
        if rand_colors is None:
            rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                           '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                           '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                           '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')
        '''
         rand_colors = ('#f67280', '#00a8cc', '#ffd082', '#fb8d62', '#6e5773', '#21bf73', '#d5c455', '#c9753d',
                       '#ad62aa','#d77fa1', '#a0855b', '#ffd800', '#da2d2d', '#6f9a8d', '#a8ff3e', '#b2fcff',
                       '#a0c334', '#b5525c', '#c06c84', '#3a3535', '#9b45e4', '#f6da63', '#9dab86', '#0c093c',
                       '#f6f078', '#64c4ed', '#da4302', '#5edfff', '#08ffc8', '#ca3e47', '#f7ff56', '#6c5ce7')
        '''
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
            marker.geneplot_mhat(df, markeridcol, chr, pv, gwasp, markernames, gfont, gstyle, ax=ax, plotlabelrotation =  plotlabelrotation)
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
        #general.get_figure(show, r, figtype, figname)
        return fig, ax
    
    
    def miami(df="dataframe", chromo=None, logp1=None, logp2=None, color=None, dim=(10,10), r=300, ar=90, gwas_sign_line=False,
             gwasp=5E-08, dotsize=8, markeridcol=None, markernames=None, gfont=8, valpha=1, show=False, figtype='png',
             axxlabel=None, axylabel=None, axlabelfontsize=9, axlabelfontname="Arial", axtickfontsize=9, figtitle = 'miami plot',
             label1='firstgroup', label2 = 'secondgroup', rand_colors = None,
             axtickfontname="Arial", ylm=None, gstyle=1, yskip = 1, plotlabelrotation = 0, figname='miami', invert = False, fig = None, ax = None):

        _x, _y = 'Chromosomes', r'$ -log_{10}(P)$'

        df['tpval'] = df[logp1]
        df['tpval2'] = -df[logp2]
        df = df.sort_values(chromo)

        df['ind'] = range(len(df))
        df_group = df.groupby(chromo)

        if rand_colors is None:
            rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                               '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                               '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                               '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')
        #color_list = sample(rand_colors, df[chromo].nunique())
        color_list = rand_colors[:df[chromo].nunique()]

        xlabels = []
        xticks = []
        
        if fig is None:
            #fig, ax = plt.subplots(figsize = dim)
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
        #ax.set_ylabel(_y, fontsize=axlabelfontsize, fontname=axlabelfontname)
        
        ax0.text(.5,0,_y, fontsize=axlabelfontsize, fontname=axlabelfontname, rotation = 90, va = 'center')
        
        plt.title(figtitle)
        #general.get_figure(show, r, figtype, figname)
        return fig, ax












