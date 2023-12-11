#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 01:00:04 2019

@author: peter
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import math

sns.set(context='paper')

def GFID(myDict):
    """Get First in Dictionary
    
    A function to strip the first and assumed only eliment out of a
    dictionary.
    
    Args:
       myDict (Dict):  Dictioary wraped around single member.
       
    Returns:
       what ever was in the dictionary.
    """
    return (myDict[next(iter(myDict))])

def unevenSubplots(number, col, sharex=True, sharey=False, 
                   figsize=(12,10)):
    fig = plt.figure(figsize=figsize)
    rows = number//col-(number%col==0)+1
    axs = []
    for i in range(number):
        if sharex and i>0 and (not sharey):
            axs.append(fig.add_subplot(rows, col, i+1, sharex=axs[0]))
        elif (not sharex) and i>0 and sharey:
            axs.append(fig.add_subplot(rows, col, i+1, sharey=axs[0]))
        elif sharex and i>0 and sharey:
            axs.append(fig.add_subplot(rows, col, i+1, sharex=axs[0],
                                       sharey=axs[0]))
        else:
            axs.append(fig.add_subplot(rows, col, i+1))
    if sharey:
        for i, ax in enumerate(axs):
            if i%col>0:
                for label in ax.get_yticklabels():
                    label.set_visible(False)
                ax.yaxis.offsetText.set_visible(False)
                ax.yaxis.label.set_visible(False)
    if sharex:
        for i, ax in enumerate(axs):
            if i+col < number:
                for label in ax.get_xticklabels():
                    label.set_visible(False)
                ax.xaxis.offsetText.set_visible(False)
                ax.xaxis.label.set_visible(False)
    return fig, axs

def plotPanels(df, xVal, yVal, panelVal, hueVal = None, save=None):
    panVals = pd.unique(df[panelVal])
    if hueVal is not None:
        hueVals = pd.unique(df[hueVal])
    else:
        hueVals = [None]
    rows = math.ceil(math.sqrt(3*len(panVals)/5))
    cols = len(panVals)//rows+(len(panVals)%rows>0)
    fig, axs = unevenSubplots(len(panVals), cols)
    for pan, ax in zip(panVals, axs):
        ax.set_title(pan)
        for hueIndex, theHueVal in enumerate(hueVals):
            subDF = df[df[panelVal]==pan]
            if hueVal is not None: 
                subDF = subDF[subDF[hueVal]==theHueVal]
            for _, g in subDF.groupby((subDF[xVal].diff() < 0).cumsum()):
                ax.plot(g[xVal],g[yVal],"C"+str(hueIndex))
    fig.tight_layout()
    if save is not None:
        fig.savefig(save)

def trim_axs(axs, N):
    """removes unwanted subplots from subplots function
    
    Convenience function used to remove unnneeded axises generated
    by the subplots function and return them in a flat structure
    
    Args:
       axs (array of Axes objects):  output axises from subplots
       N (int): number of subplots you actualy want.
       
    Returns:
       array of Axes objects: your remaining axises
    """
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]

def getTCSelectionMaxes(TCList,selectionList=None,varSelection=None,
                       valRemove=[]):
    """Get TimeCourse maximum values by variable
    
    Gets the maximum value for each variable over all TimeCouses and
    all their time points. Useful in ordering and grouping multi plots.
    
    Args:
       TCList (list):  A list of TimeCourse Dataframes
       
    Kwargs:
       selectionList (list of ints): an optional list of the TimeCourses
           (by their indexing in the TCList) to be included in the
           calculation.
       varSelection (list of str): an optional list of variables from the
           TimeCourses we want returned in the results.
       valRemove (list of float): an optional list of values that are to be
           excluded from the list. A quick and dirty way to exclude
           variables that have all been set to some pre determined value.
       
    Returns:
       serise: the maximum value found for each variable
    """
    if selectionList is None:
        maxes = [course.max() for course in TCList]
    else:
        maxes = [TCList[selection].max() for selection in selectionList]
    maxes = pd.DataFrame(maxes).max()
    maxes = maxes.drop(labels=['Time'])
    if varSelection is not None:
        maxes = maxes.filter(items=varSelection)
    maxes = maxes.sort_values(ascending=False)
    for value in valRemove:
        maxes = maxes[maxes!=value]
    return maxes

def breakSeriesByScale(mySerise,relativeScaleBreakPoint=0.1,
                       maxRunLength=6):
    """breaks Series into subparts of similar size
    
    Takes a serise and breaks it into parts of comparable value size.
    Sub parts are returned sorted internaly and relative to each other in
    decending order. Useful for ordering and grouping plots with very
    difrent value ranges.
    
    Args:
       mySerise (serise):  A serise of floats to break into parts
       
    Kwargs:
       relativeScaleBreakPoint (float): Specifies the smallest (by scale)
           drop in value that should force an entry into a new subpart.
       maxRunLength (int): Max size of any subpart.
       
    Returns:
       list of list of str: the indexes of the serise are returned broken
       into subparts (the inner list) acording the the size of their
       asociated values.
    """
    outerList=[]
    innerList=[]
    lastVal=None
    for index, value in mySerise.sort_values(ascending=False).items():
        if lastVal is None:
            innerList.append(index)
            lastVal=value
            continue
        if (lastVal*relativeScaleBreakPoint>value or 
            len(innerList)==maxRunLength):
            outerList.append(innerList)
            innerList=[]
        else:
            innerList.append(index)
        lastVal=value
    if len(innerList)>0:
        outerList.append(innerList)
    outerList = [myList for myList in outerList if len(myList)>=1]
    return outerList

class profileLikelyhoodVisualisor:
    def __init__(self,profLike):
        idenDict = profLike[0]
        refDict = profLike[1]
        self.scaledDict = {}
        self.fitDict = {}
        self.aggDF = []
        for idenPage in idenDict.keys():
            temp = idenDict[idenPage].copy()
            temp["variable"] = idenPage
            temp = temp.sort_values(by=["RSS"])
            temp["LEVEL"] = temp.groupby(idenPage).cumcount()
            temp = temp.sort_values(by=["LEVEL",idenPage])
            temp = temp.reset_index(drop=True)
            self.aggDF.append(temp)            
        self.aggDF = pd.concat(self.aggDF,ignore_index=True)
        for variable, origVal in refDict.items():
            if origVal!=0 and variable in self.aggDF.columns:
                self.aggDF[variable] = self.aggDF[variable]/origVal
        for idenPage in idenDict.keys():
            self.scaledDict[idenPage] = idenDict[idenPage].copy()
            temp = self.scaledDict[idenPage][idenPage]==refDict[idenPage]
            refRow = self.scaledDict[idenPage][temp].squeeze()
            if not isinstance(refRow,pd.Series):
                try:
                    refRow = refRow.iloc[0]
                except:
                    print(idenPage)
                    print(refRow)
            for col in refRow.index:
                if refRow[col]!=0:
                    self.scaledDict[idenPage][col] = (
                            self.scaledDict[idenPage][col]/refRow[col])
            self.fitDict[idenPage] = []
            for col in [i for i in self.scaledDict[idenPage].columns
                        if i in idenDict.keys() and i!=idenPage]:
                df = self.scaledDict[idenPage][[idenPage]].copy()
                df["myOnes"] = 1
                m, c = np.linalg.lstsq(df[[idenPage,"myOnes"]],
                                       self.scaledDict[idenPage][col])[0]
                Sq = ((self.scaledDict[idenPage][col]-1)**2).sum()
                self.fitDict[idenPage].append({"m":m, "c":c, "col":col,
                            "sq":Sq})
            #self.fitDict[idenPage].sort(key=(lambda x:np.abs(x["m"])))
            self.fitDict[idenPage].sort(key=(lambda x:np.abs(x["m"])),
                        reverse=True)
        self.displayOrder = [{"k":k, "o":np.abs(v[0]["m"])} for k,v
                     in self.fitDict.items()]
        self.displayOrder.sort(key=(lambda x:x["o"]),reverse=True)
        self.displayOrder = [i["k"] for i in self.displayOrder]
        self.displayOrder = [self.displayOrder[4*i:4*(i+1)] for i
                             in range(len(self.displayOrder))]
        self.displayOrder = [i for i in self.displayOrder if len(i)>0]
    
    def plotProfiles(self,showRows,showLimit=5, save = None,
                     style = None):
        if style is None:
            myStyle = "darkgrid"
        else:
            myStyle = style
        with sns.axes_style(style):
            fig, axs = plt.subplots(len(showRows), showLimit,
                                    sharex=True, figsize=(12,10))
            for row, i in zip(showRows,range(len(showRows))):
                for col, j in zip([i["col"] for i in self.fitDict[row]],
                                  range(showLimit)):
                    axs[i][j].scatter(self.scaledDict[row][row],
                       self.scaledDict[row][col])
                    axs[i][j].set_title(col)
                axs[i][0].set_ylabel(row, size='large')
            fig.tight_layout()
            if save is not None:
                fig.savefig(save)
            
    def plotRSS(self,showVars,showLimit=5, save = None, style = None):
        if style is None:
            myStyle = "darkgrid"
        else:
            myStyle = style
        with sns.axes_style(style):
            fig, axs = plt.subplots((len(showVars)-1)//showLimit+1,
                                    showLimit, sharex=True,
                                    figsize=(12,10))
            axs = trim_axs(axs, len(showVars))
            for myVar, i, ax in zip(showVars,range(len(showVars)),axs):
                temp = self.aggDF[self.aggDF["variable"]==myVar]
                temp = temp[[myVar,"RSS","LEVEL"]]
                avalableLevels = temp["LEVEL"].unique()
                for myLevel in avalableLevels:
                    ax.plot(temp[temp["LEVEL"]==myLevel][myVar],
                            temp[temp["LEVEL"]==myLevel]["RSS"])
                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_title(myVar)
            fig.tight_layout()
            if save is not None:
                fig.savefig(save)

class timeCourseVisualiser:
    def __init__(self,data):
        if not isinstance(data, list):
            data = [data]
        runningList=[]
        for i in range(len(data)):
            valueColumns=list(data[i].columns)
            if "Time" in valueColumns:
                valueColumns.remove("Time")
            if "" in valueColumns:
                valueColumns.remove("")
            dataFrame=pd.melt(data[i],id_vars=["Time"],
                    value_vars=valueColumns)
            dataFrame['index']=i
            runningList.append(dataFrame)
        self.longData = pd.concat(runningList,ignore_index=True)
    
    def multiPlot2(self,indexSelect=None,varSelect=None,wrapNumber=5,
                  compLines=None, save = None):
        df = self.longData
        if varSelect is not None:
            if not isinstance(varSelect, list):
                varSelect = [varSelect]
            df = df[df['variable'].isin(varSelect)]
        if indexSelect is not None:
            if not isinstance(indexSelect, list):
                indexSelect = [indexSelect]
            df=df.loc[df['index'].isin(indexSelect)]
        indexes = list(df['index'].unique())
        colors = sns.husl_palette(len(indexes)).as_hex()
        if compLines is not None:
            compVars=list(compLines.columns)
            dfB = compLines.copy()
            dfB["Time"]=dfB.index
            dfB = pd.melt(dfB,id_vars=["Time"],
                          value_vars=compVars)
            dfB["index"]=-1
            df = pd.concat([dfB,df],ignore_index=True)
            colors.append("#000000")
            indexes.append(-1)
        #colors = dict(zip(indexes, colors))
        grid = sns.FacetGrid(df, col="variable", col_wrap=wrapNumber, 
                             palette=colors)
        grid.map(sns.lineplot,"Time","value","index")
        if save is not None:
            grid.savefig(save)
            
    def multiPlot(self,indexSelect=None,varSelect=None,wrapNumber=5,
                  compLines=None, save = None, xlim = None, ylim = None,
                  forceYAxisZero = True, colourOverride = None,
                  style = None, legend = None, varAsAxis = False,
                  xAxisLabel = None, yAxisLabel = None, figsize = (12,10),
                  legendLoc = 'lower right', compErrors=None):
        """Plots grid if time course variables
        
        Creats grid of rainbow coloured plots for each variable in
        colection of TimeCourses. Optionaly plots refrence data against
        it as black dots.
        
        Kwargs:
           indexSelect (list of ints): Indexs of time courses that should
               be included in the plot.
           varSelect (list of str): Variables from the Time course that
               should be graphed.
           wrapNumber (int): max number of columns to alow in grid.
           compLines (dataFrame): refrence data to graph against time
               courses. Sould be in wide format and have the times as the
               index.
           save (str): A path to where to save the image. if omited image
               not saved
        """
        if isinstance(compLines,list):
            compVars = [list(i.columns) for i in compLines]
            dfB = [i.copy() for i in compLines]
            for i in range(len(compLines)):
                if "Time" not in compVars[i]:
                    dfB[i]["Time"]=dfB[i].index
                else:
                    compVars[i].remove("Time")
                dfB[i] = pd.melt(dfB[i], id_vars=["Time"],
                   value_vars=compVars[i])
        elif compLines is not None:
            compVars=list(compLines.columns)
            dfB = compLines.copy()
            if "Time" not in compVars:
                dfB["Time"]=dfB.index
            else:
                compVars.remove("Time")
            dfB = pd.melt(dfB,id_vars=["Time"],
                          value_vars=compVars)
        if isinstance(compErrors,list):
            dfC = [i.copy() for i in compErrors]
            for i in range(len(dfC)):
                if "Time" not in dfC[i].columns:
                    dfC[i]["Time"]=dfC[i].index
                dfC[i] = pd.melt(dfC[i], id_vars=["Time"],
                   value_vars=[j for j in dfC[i].columns if j!="Time"])
        elif compErrors is not None:
            dfC = compErrors.copy()
            dfC = pd.melt(dfC, id_vars=["Time"],
                value_vars=[j for j in dfC.columns if j!="Time"])
        if varSelect is None:
            varSelect=list(self.longData['variable'].unique())
        if indexSelect is None:
            indexSelect=list(self.longData['index'].unique())
        if not isinstance(indexSelect,list):
            indexSelect = [indexSelect]
        if len(varSelect)<wrapNumber:
            #cols = math.floor(math.sqrt(len(varSelect)))
            cols = math.ceil(math.sqrt(len(varSelect)))
        else:
            cols = wrapNumber
        rows = math.ceil(len(varSelect)/cols)
        if style is None:
            myStyle = "darkgrid"
        else:
            myStyle = style
        with sns.axes_style(style):
            fig, axs = plt.subplots(rows, cols, sharex=True,
                                        figsize=figsize)
            if (rows>1):
                axs = trim_axs(axs, len(varSelect))
            elif (cols==1):
                axs = [axs]
            if colourOverride is not None:
                myColorMap = plt.get_cmap(name="cool")
            else:
                myColorMap = plt.get_cmap(name="hsv",
                                          lut=len(indexSelect)+1)
            for ax, theVar, j in zip(axs, varSelect, range(len(varSelect))):
                if varAsAxis:
                    if isinstance(yAxisLabel,list):
                        ax.set_ylabel(theVar+" "+yAxisLabel[j])
                    elif yAxisLabel is not None:
                        ax.set_ylabel(theVar+" "+yAxisLabel)
                    else:
                        ax.set_ylabel(theVar)
                else:
                    ax.set_title(theVar)
                    if isinstance(yAxisLabel,list):
                        ax.set_ylabel(yAxisLabel[j])
                    elif yAxisLabel is not None:
                        ax.set_ylabel(yAxisLabel)
                if xAxisLabel is not None:
                    ax.set_xlabel(xAxisLabel)
                df = self.longData
                df = df[df['variable']==theVar]
                if indexSelect is not None:
                    for theIndex, i in zip(indexSelect,
                                           range(len(indexSelect))):
                        df2 = df[df['index']==theIndex]
                        if colourOverride is not None:
                            ax.plot(df2["Time"], df2["value"],
                                    linestyle='solid',
                                    color=myColorMap(colourOverride[i]))
                        else:
                            ax.plot(df2["Time"], df2["value"],
                                    linestyle='solid',
                                    color=myColorMap(i))
                if isinstance(compLines,list):
                    for i, theIndex in enumerate(indexSelect):
                        dfB2 = dfB[theIndex][
                                dfB[theIndex]['variable']==theVar]
                        if isinstance(compErrors,list):
                            dfC2 = dfC[theIndex][
                                dfC[theIndex]['variable']==theVar]
                            dfC2 = dfC2["value"]
                        else:
                            dfC2 = None
                        if colourOverride is not None:
                            #ax.plot(dfB2["Time"], dfB2["value"],"o",
                            #        color=myColorMap(colourOverride[i]))
                            ax.errorbar(dfB2["Time"], dfB2["value"],
                                        yerr = dfC2, xerr = dfC2, ls='none',
                                        color=myColorMap(colourOverride[i]),
                                        marker='o')
                        else:
                            #ax.plot(dfB2["Time"], dfB2["value"],"o",
                            #        color=myColorMap(i))
                            ax.errorbar(dfB2["Time"], dfB2["value"],
                                        yerr = dfC2, xerr = dfC2, ls='none',
                                        color=myColorMap(i), marker='o')
                elif compLines is not None:
                    dfB2 = dfB[dfB['variable']==theVar]
                    if compErrors is not None:
                        dfC2 = dfC[dfC['variable']==theVar]
                        dfC2 = dfC2["value"]
                    else:
                        dfC2 = None
                    #ax.plot(dfB2["Time"], dfB2["value"],"ko")
                    ax.errorbar(dfB2["Time"], dfB2["value"],
                                yerr = dfC2, xerr = dfC2, ls='none',
                                color="k", marker='o')
                if xlim is not None:
                    ax.set_xlim(xlim)
                if isinstance(ylim, dict):
                    if theVar in ylim:
                        ax.set_ylim(ylim[theVar])
                elif ylim is not None:
                    ax.set_ylim(ylim)
                elif forceYAxisZero:
                    ax.set_ylim([0, None])
            if legend is not None:
                if colourOverride is not None:
                    custom_lines = [Line2D([0], [0], color=myColorMap(
                            colourOverride[i]), lw=4)
                            for i in range(len(indexSelect))]
                else:
                    custom_lines = [Line2D([0], [0], color=myColorMap(i),
                                           lw=4)
                            for i in range(len(indexSelect))]
                if ((not isinstance(compLines,list)) and
                    (compLines is not None)):
                    custom_lines.append(Line2D([0], [0], 
                                               color="k", lw=4))
                fig.legend(custom_lines, legend,
                           loc = legendLoc)
            fig.tight_layout()
            if save is not None:
                fig.savefig(save)
            
    def barChart(self, time, indexSelect=None, varSelect=None,
                 wrapNumber=5, compLines=None, save = None,
                 style = None, colourOveride = None, varOnAxis = False,
                 spacing = None, figsize=(12,10)):
        if compLines is not None:
            compVars=list(compLines.columns)
            dfB = compLines.copy()
            if "Time" not in compVars:
                dfB["Time"]=dfB.index
            else:
                compVars.remove("Time")
            dfB = pd.melt(dfB,id_vars=["Time"],
                          value_vars=compVars)
        if varSelect is None:
            varSelect=list(self.longData['variable'].unique())
        if indexSelect is None:
            indexSelect=list(self.longData['index'].unique())
        if isinstance(indexSelect,dict):
            renameIndex = [v for _,v in indexSelect.items()]
            indexSelect = [k for k,_ in indexSelect.items()]
        elif not isinstance(indexSelect,list):
            indexSelect = [indexSelect]
            renameIndex = None
        if len(varSelect)<wrapNumber:
            cols = math.floor(math.sqrt(len(varSelect)))
        else:
            cols = wrapNumber
        rows = math.ceil(len(varSelect)/cols)
        if style is None:
            myStyle = "darkgrid"
        else:
            myStyle = style
        with sns.axes_style(myStyle):
            fig, axs = plt.subplots(rows, cols, sharex=True, 
                                    figsize=figsize)
            if (rows>1):
                axs = trim_axs(axs, len(varSelect))
            elif (cols==1):
                axs = [axs]
            # may need to add black
            myColorMap = plt.get_cmap(name="hsv", lut=len(indexSelect)+1)
            
            for ax, theVar in zip(axs, varSelect):
                df = self.longData
                df = df[df['variable']==theVar]
                df = df[df['Time']==time]
                df = df[df['index'].isin(indexSelect)]
                if spacing is not None:
                    bar_pos = np.arange(len(df['index'])+sum(spacing))
                    bar_pos = [i for i,j in zip(bar_pos,spacing) if not j]
                else:
                    bar_pos = np.arange(len(df['index']))
                if colourOveride is None:
                    colorList = [[j for j,x in enumerate(indexSelect) 
                                  if x == i][0] for i in df['index']]
                    colorList = [myColorMap(i) for i in colorList]
                else:
                    colorList = colourOveride
                ax.bar(bar_pos, df["value"], color=colorList)
                if varOnAxis:
                    ax.set_ylabel(theVar)
                else:
                    ax.set_title(theVar)
                ax.set_xticks(bar_pos)
                ax.set_xticklabels(renameIndex)
                if compLines is not None:
                    df = dfB[dfB["Time"] == time]
                    df = df[df['variable'] == theVar]
                    if len(df)==1:
                        ax.axhline(y=df.squeeze()["value"])
            fig.tight_layout()
            if save is not None:
                fig.savefig(save)
        
class parameterEstimationVisualiser:
    def __init__(self,data):
        if not isinstance(data, list):
            data = [data]
        BPList=[]
        RSSList=[]
        wideList=[]
        for i in range(len(data)):
            df = data[i][list(data[i])[0]].copy()
            theColumns=list(df.columns)
            if "" in theColumns:
                df=df.drop(columns="")
            theColumns=list(df.columns).remove('RSS')
            df['subIndex'] = df.index
            df=pd.melt(df,id_vars=['subIndex','RSS'],
                       value_vars=theColumns)
            df['index']=i
            BPList.append(df)
            df = data[i][list(data[i])[0]].copy()
            df['subIndex'] = df.index
            df=df[['subIndex',"RSS"]]
            df['index']=i
            RSSList.append(df)
            df = data[i][list(data[i])[0]].copy()
            df['index']=i
            wideList.append(df)
        self.BPData = pd.concat(BPList,ignore_index=True)
        self.RSSData = pd.concat(RSSList,ignore_index=True)
        self.wideData = pd.concat(wideList,ignore_index=True)
        
    def violinPlot(self,indexSelect=None,paramSelect=None,RSSSelect=None,
                   save = None, indexNames=None):
        """Creates Violin plot showing variation in paramiter space
        
        Kwargs:
           indexSelect (list of ints): Indexs of sets of paramiter
               estimation results that should be ploted. Takes a maximum of
               2 to allow comparision of two sets. optionaly a single
               interger may be used. Omision leads to all sets of paramiter
               estimations being processed in a combined way.
           varSelect (list of str): Variables from the paramiter
               estimations that should be graphed.
           RSSSelect (list of floats): Specifies a maximum RSS value for
               data in the paramiter estimations to be included in the plot.
               Ordering should match indexSelect.
           save (str): A path to where to save the image. if omited image
               not saved
           indexNames (list of str): names to use in the plot legend. If
               Omited numerical index is used. Ordering should match
               indexSelect.
        """
        df=self.BPData
        if RSSSelect is not None:
            if not isinstance(RSSSelect, list):
                RSSSelect = [RSSSelect]
            if indexSelect is None:
                df = df.loc[df['RSS']<=RSSSelect[0]]
        if indexSelect is not None:
            if not isinstance(indexSelect, list):
                indexSelect = [indexSelect]
            if RSSSelect is not None:
                for i in range(len(indexSelect)):
                    df = df.loc[(df['RSS']<=RSSSelect[i]) |
                            (df['index']!=indexSelect[i])]
            df = df.loc[df['index'].isin(indexSelect)]
        if paramSelect is not None:
            if not isinstance(paramSelect, list):
                paramSelect = [paramSelect]
            df = df.loc[df['variable'].isin(paramSelect)]
        plt.figure()
        if indexNames is not None:
            if indexSelect is not None:
                if len(indexSelect)!=len(indexNames):
                    return None
            df = df.copy()
            df["index"] = df["index"].replace(
                    dict(zip(indexSelect,indexNames)))
        if indexSelect is None:
            if paramSelect is None:
                vp = sns.violinplot(x="variable", y="value", data=df, cut=0)
            else:
                vp = sns.violinplot(x="variable", y="value", data=df,
                                    order=paramSelect, cut=0)
        else:
            if paramSelect is None:
                vp = sns.violinplot(x="variable", y="value", data=df,
                                    hue="index", split=True, cut=0)
            else:
                vp = sns.violinplot(x="variable", y="value", data=df,
                                    hue="index", split=True,
                                    order=paramSelect, cut=0)
        if save is not None:
            vp.get_figure().savefig(save)
            
    def waterFall(self,save = None, indexNames=None, style = None):
        """creates watrfall plot of RSS values in paramiter estimations
        
        Kwargs:
           save (str): A path to where to save the image. if omited image
               not saved
           indexNames (list of str): names to use in the plot legend.
        """
        if indexNames is not None:
            if not isinstance(indexNames,list):
                indexNames = [indexNames]
            df = self.RSSData.copy()
            df["index"] = [indexNames[x] for x in df["index"]]
        else:
            df = self.RSSData
        if style is None:
            myStyle = "darkgrid"
        else:
            myStyle = style
        with sns.axes_style(style):
            plt.figure()
            lp = sns.lineplot(data=df,x="subIndex",y="RSS",hue="index")
            if save is not None:
                lp.get_figure().savefig(save)
            
    def refPointVsRSS(self, myRef, save=None, indexNames=None):
        names=list(myRef.index)
        refPoint = myRef.copy()
        refPoint["index"] = 1
        refPoint["RSS"] = 1
        diffTable = self.wideData.copy()
        diffTable = (diffTable/refPoint)-1
        diffTable["index"] = diffTable["index"]+1
        diffTable["index"] = diffTable["index"].astype(int)
        diffTable["RSS"] = diffTable["RSS"]+1
        diffTable["difrence"] = 0
        for name in names:
            diffTable["difrence"] = (diffTable["difrence"] +
                     diffTable[name]*diffTable[name])
        diffTable["difrence"] = diffTable["difrence"] / len(names)
        if indexNames is not None:
            if not isinstance(indexNames,list):
                indexNames = [indexNames]   
            diffTable["index"] = [indexNames[x] for x in diffTable["index"]]
        plt.figure()
        sp = sns.scatterplot(data=diffTable, x="RSS", y="difrence",
                             hue="index")
        if save is not None:
            sp.get_figure().savefig(save)