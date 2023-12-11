#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 11:33:49 2022

@author: peter
"""

from python.utilityTools import *
from python.visualisationTools import *
from GPUtilityFunctions import importData, scale, tableToDict, normalize
import json

cmdDefaults = {"dataName":"InnsBrokenPaper",
               "maxPerGraph":15,
               "scaleToMatch":"yes",
               "outScaleName":"scales.csv",
               #"inScale":"scales.csv",
               "cmdRec":"comandRecord3.p",
               #"cmdRec":"comandRecordZR75.p",
               "elements":"elements3.p",
               #"elements":"elements_ZR75.p",
               "id":"",
               #"id":"ZR75",
               #"PL":"profile_ZR75.csv",
               "style":"white"
               }

def TCtoPoints(tc,df,pTimeVar):
    t = [i for i in range(len(tc)) if tc.iloc[i][pTimeVar]==0][-1]
    t = tc.iloc[t:].copy()
    t = t.drop(columns=["Time"]).rename(columns={pTimeVar:"Time"})
    t = [t[t["Time"]>=i].iloc[0] for i in df.index]
    t = pd.DataFrame(t)
    t.index = df.index
    return t

def getBestScale(TCdict,dfDict):
    myKeys = set(TCdict.keys()).intersection(set(dfDict.keys()))
    outScale = {}
    for col in set.union(*[set(i.columns) for _,i in dfDict.items()]):
        num = 0
        denom = 0
        for k in myKeys:
            if (col in dfDict[k].columns) and (col in TCdict[k].columns):
                num += (dfDict[k][col]*TCdict[k][col]).sum()
                denom += (TCdict[k][col]*TCdict[k][col]).sum()
        if denom!=0:
            outScale[col] = num/denom
    # sum_i|TC_i*s-df_i|=m
    return outScale

def unpackTables(df,indexKey="index",keyVar="dataKey"):
    tempList = []
    i = 0
    while True:
        temp = tableToDict(df,
                           indexKey = indexKey,
                           indexVal = i,
                           keyVar = keyVar)
        if temp != {}:
            tempList.append(temp)
            i += 1
        else:
            break
    return tempList

def normByPoints(TC,df):
    if "time" in df.columns:
        TP = list(df["time"])
    elif "Time" in df.columns:
        TP = list(df["Time"])
    else:
        TP = list(df.index)
    TP = [[j for j,r in TC.iterrows() if r["Time"]>=i][0] for i in TP]
    means = TC.loc[TP].mean().to_dict()
    means = {k:v for k,v in means.items() if k in df.columns 
             and k.lower()!="time"}
    outTC = TC.copy()
    for k,v in means.items():
        outTC[k] = outTC[k]/v
    return outTC, means

def graphMeansRatios(dfMeans,meanData):
    md = meanData.copy()
    md["df"] = dfMeans
    plotData = {}
    for dk,df in md.items(): 
        obs = set.intersection(*[set(v.keys()) for k,v in df.items()])
        dfRatios = {}
        for i in obs:
            t = {k:max(v[i],0) for k,v in df.items()}
            u = sum([j for _,j in t.items()])
            if u>0:
                dfRatios[i] = {k:v/u for k,v in t.items()}
            else:
                dfRatios[i] = {k:0 for k,_ in t.items()}
        plotData[dk] = dfRatios
    obs = set.intersection(*[set(v.keys()) for k,v in plotData.items()])
    order = list(meanData.keys())
    order.sort(key=float)
    order = ["df"] + order
    myKeys = set.union(*[set.union(*[set(j.keys()) for _,j in i.items()]) 
                         for _,i in plotData.items()])
    plotData = {k:{j:[plotData[i][k][j] for i in order] for j in myKeys} 
                for k in obs}
    
    fig, axs = plt.subplots(nrows=1+(len(obs)-(len(obs)-1)%5)//5, ncols=5)
    axs = [j for i in axs for j in i]
    order = [str(i) for i in order]
    for myObs, ax in zip(obs,axs):
        bottom = [0 for _ in range(len(order))]
        for myKey in myKeys:
            ax.bar(order, plotData[myObs][myKey], 1, bottom=bottom, 
                   label=myKey)
            bottom = [bottom[i]+j for i,j in enumerate(plotData[myObs][myKey])]
        ax.set_ylabel(myObs[:-3]+" (AU)")
    plt.tight_layout()
    
def plotMyData(dictOfTC,dictOfPoints,varSelect,index,save=None, style = None):
    frameVars = [varSelect[15*i:15*(i+1)] for i in range(len(varSelect))]
    frameVars = [i for i in frameVars if len(i)>0]
    myColorMap = plt.get_cmap(name="hsv", lut=len(dictOfTC)+1)
    myColorMap = {k:myColorMap(i) for i,k in enumerate(dictOfTC.keys())}
    if style is None:
            myStyle = "darkgrid"
    else:
        myStyle = style
    with sns.axes_style(myStyle):
        for j, myVars in enumerate(frameVars):
            fig, axs = plt.subplots(math.ceil(len(myVars)/5), 
                                    min(5,len(myVars)), 
                                    sharex=True, figsize=(12,10))
            if not isinstance(axs,np.ndarray):
                axs = [axs]
            else:
                axs = axs.flat
            for var, ax in zip(myVars, axs):
                ax.set_ylabel(var)
                for k, TC in dictOfTC.items():
                    if var in TC.columns:
                        ax.plot(TC["Time"], TC[var], color=myColorMap[k],
                                linestyle='solid')
                    if k in dictOfPoints and var in dictOfPoints[k].columns:
                        ax.errorbar(dictOfPoints[k].index,
                                    dictOfPoints[k][var],
                                    ls='none', marker='o', color=myColorMap[k])
            fig.suptitle("index: " + str(index), fontsize=16)
            custom_lines = {k:Line2D([0], [0], color=v, lw=4) 
                            for k, v in myColorMap.items()}
            legend = list(custom_lines.keys())
            custom_lines = [custom_lines[i] for i in custom_lines]
            fig.legend(custom_lines, legend)
            fig.tight_layout()
            if save is not None:
                fig.savefig(save+"-i"+str(index)+"f"+str(j)+".png")
            
def homeBrewRSS(sim,ex):
    t = set.intersection(set(sim.keys()),set(ex.keys()))
    t = {k:(sim[k],ex[k]) for k in t}
    RSS = 0
    for _, (i,j) in t.items():
        for myTime, exRow in j.iterrows():
            for _, simRow in i.iterrows():
                if simRow["Time"]>=myTime:
                    break
            for k,v in exRow.items():
                if k in simRow:
                    RSS += (v-simRow[k])**2
    return RSS
    

if __name__ == "__main__":
    cmdDict, cmdFlags = getCmdLineArgs()
    cmdDict = parseCmdDict(cmdDict,cmdDefaults)
    oldCmd = loadPick(["data",cmdDict["dataName"],cmdDict["cmdRec"]], 
                      relative=True)
    elements = loadPick(["data",cmdDict["dataName"],cmdDict["elements"]], 
                        relative=True)
    if "timeMul" in oldCmd["cmdDict"]:
        df = importData([oldCmd["cmdDict"]["data"]],
                        timeMultiplyer=oldCmd["cmdDict"]["timeMul"])
    else:
        df = importData([oldCmd["cmdDict"]["data"]])
    
    if "inScale" in cmdDict:
        refScales = pd.read_csv(resolvePath(["data",
                                             oldCmd["cmdDict"]["dataName"],
                                             cmdDict["inScale"]],
                                            relative=True))
    else:
        refScales = None
    
    dataKey = None
    if "configFile" in oldCmd["cmdDict"]:
        f = open(oldCmd["cmdDict"]["configFile"],'r')
        myJ = f.read()
        f.close()
        myJ = json.loads(myJ)
        if "dataKey" in myJ:
            dataKey = myJ["dataKey"]
    if dataKey is None:     
        dataKey = pd.read_csv(resolvePath([oldCmd["cmdDict"]["dataKey"]],
                                          relative=True),
                              header=None, index_col = 0).squeeze().to_dict()
    
    for k in df:
        df[k] = df[k].drop(columns=[i for i in df[k].columns
                                    if not i in dataKey])
        df[k] = df[k].rename(columns=dataKey)

    if ((not "parallelAntStr" in oldCmd["cmdDict"]) or 
        (oldCmd["cmdDict"]["parallelAntStr"]=="")):
        scales = pd.read_csv(resolvePath(["data",oldCmd["cmdDict"]["dataName"], 
                                          oldCmd["cmdDict"]["scaleOut"]], 
                                         relative=True))
        scales = scales.drop(columns=['Unnamed: 0']).squeeze().to_dict()
        df = scale(df,scales)
    elif cmdDict["scaleToMatch"].lower() in ["yes","true"]:
        df = {k.replace("-","").replace("_",""):v for k,v in df.items()}
    else:
        df = {k.replace("-","").replace("_",""):v for k,v in df.items()}
        df, dfMeans = normalize(df)
    
    TCPath = resolvePath(["data", oldCmd["cmdDict"]["dataName"],
                          oldCmd["cmdDict"]["TCOutputs"]], relative=True)
    
    timeCourse = pd.read_csv(TCPath)
    timeCourse = unpackTables(timeCourse,
                              indexKey = oldCmd["cmdDict"]["TCIndex"],
                              keyVar = oldCmd["cmdDict"]["TCKeyVar"])
        
    PPath = resolvePath(["data",oldCmd["cmdDict"]["dataName"],
                         oldCmd["cmdDict"]["paramOut"]], relative=True)
    params = pd.read_csv(PPath)
    params = params.drop(columns=["Unnamed: 0"])
    PEVis = parameterEstimationVisualiser({"notag":params})
    forceDir(["figures",oldCmd["cmdDict"]["dataName"]],relative=True)
    PEVis.waterFall(save=resolvePath(["figures",oldCmd["cmdDict"]["dataName"],
                                      "waterfall"+cmdDict["id"]+".png"], 
                                     relative=True),
                    style = cmdDict["style"])
    
    maxTime = max([j for _,i in df.items() for j in i.index])
    #if (("parallelAntStr" in oldCmd["cmdDict"]) and 
    #    (oldCmd["cmdDict"]["parallelAntStr"]!="")):
    #    maxTime = maxTime+2
    validData = {}
    meanData = {}
    scaleOut = []
    for i,d in enumerate(timeCourse):
        if i>oldCmd["cmdDict"]["TCToComp"]:
            break
        if all([v.iloc[-1][oldCmd["cmdDict"]["pTimeVar"]]>=maxTime 
                for _,v in d.items()]):
            print([v.iloc[-1][oldCmd["cmdDict"]["pTimeVar"]] for _,v in d.items()])
            validData[i]={}
            meanData[i]={}
            myInsert = {k:TCtoPoints(v,df[k],oldCmd["cmdDict"]["pTimeVar"]) 
                        for k,v in d.items() if k in df}
            scales = getBestScale(myInsert,df)
            if refScales is not None:
                t = params.iloc[i]["paramSet"]
                t = refScales[refScales["index"]==t].squeeze().to_dict()
                print(t)
                scales.update(t)
            for k,v in d.items():
                idx = [idx for idx,cond 
                       in enumerate(list(v[oldCmd["cmdDict"]["pTimeVar"]]==0)) 
                       if cond][-1]
                TC = v.iloc[idx:]
                TC = TC.drop(columns=["Time"]).rename(
                        columns={oldCmd["cmdDict"]["pTimeVar"]:"Time"})
                if ((not "parallelAntStr" in oldCmd["cmdDict"]) or 
                    (oldCmd["cmdDict"]["parallelAntStr"]=="")):
                    validData[i][k]=TC
                elif cmdDict["scaleToMatch"].lower() in ["yes","true"]:
                    validData[i][k]=TC
                    for col,scale in scales.items():
                        if col in validData[i][k].columns:
                            validData[i][k][col] = validData[i][k][col]*scale
                else:
                    validData[i][k], meanData[i][k] = normByPoints(TC,df[k])
            scales["index"]=i
            scaleOut.append(scales)
    scaleOut = pd.DataFrame(scaleOut)
    if "outScaleName" in cmdDict:
        scaleOut.to_csv(resolvePath(["data", oldCmd["cmdDict"]["dataName"],
                                     cmdDict["outScaleName"]], relative=True))    
    t = set.union(*[set(i.keys()) for _,i in validData.items()])
    validData = {k:{i:d[k] for i,d in validData.items() if k in d} for k in t}
    if (("parallelAntStr" in oldCmd["cmdDict"]) and 
        (oldCmd["cmdDict"]["parallelAntStr"]!="")):
        yLim = [0,2]
    else:
        yLim = None
    emptyRSS = {}
    for k, v in validData.items():
        TC = [i for _,i in v.items()]
        TCV = timeCourseVisualiser(TC)
        varSel = [i for i in df[k].columns 
                  if i in elements["pre elements"]["assignments"]]
        j = [i for i in df[k].columns 
             if i in elements["pre elements"]["assignments"]]
        df[k] = df[k][j]
        TCV.multiPlot(varSelect=varSel, compLines=df[k], ylim = yLim,
                      style = cmdDict["style"],
                      save=resolvePath(["figures",
                                        oldCmd["cmdDict"]["dataName"],
                                        k+"-TimeCourse"+cmdDict["id"]+
                                        "-wb.png"], 
                                       relative=True))
        if ((not "parallelAntStr" in oldCmd["cmdDict"]) or 
            (oldCmd["cmdDict"]["parallelAntStr"]=="")):
            varSel = elements["pre elements"]["metabolites"]
            varSel = [varSel[i*cmdDict["maxPerGraph"]:
                      (i+1)*cmdDict["maxPerGraph"]] 
                      for i in range(1+len(varSel)//cmdDict["maxPerGraph"])]
            varSel = [i for i in varSel if len(i)>0]
            for i, v in enumerate(varSel):
                pass
                TCV.multiPlot(varSelect=v, style = cmdDict["style"],
                              save = resolvePath(
                                 ["figures", oldCmd["cmdDict"]["dataName"],
                                  k+"-TimeCourse"+cmdDict["id"]+"-"+str(i)+
                                  ".png"], relative=True))
    if ((not "parallelAntStr" in oldCmd["cmdDict"]) or 
        (oldCmd["cmdDict"]["parallelAntStr"]=="") or 
        (cmdDict["scaleToMatch"].lower() in ["yes","true"])):
        pass
    else:
        graphMeansRatios(dfMeans,meanData)
    t = set.union(*[set(i.keys()) for _,i in validData.items()])
    validData = {i:{k:d[i] for k,d in validData.items() if i in d} for i in t}
    varSel = list(set.union(*[set(i.columns) for _,i in df.items()]))
    for i, d in validData.items():
        emptyRSS[i] = homeBrewRSS(d,df)
        plotMyData(d,df,varSel,i,save =
                   resolvePath(["figures", oldCmd["cmdDict"]["dataName"],
                                "CaseTimeCourse"+cmdDict["id"]],relative=True),
                   style = cmdDict["style"])
    emptyRSS = [[k,v] for k,v in emptyRSS.items()]
    emptyRSS.sort(key=lambda x:x[1])
    t = resolvePath(["figures", oldCmd["cmdDict"]["dataName"],"bestFit"+
                     cmdDict["id"]+".csv"],
                    relative=True)
    f = open(t,"w")
    f.write("index, RSS\n")
    for i in emptyRSS:
        f.write(", ".join([str(j) for j in i])+"\n")
    f.close()
    
    profile = resolvePath(["data", oldCmd["cmdDict"]["dataName"],
                           cmdDict["PL"]], relative=True)
    if os.path.exists(profile):
        profile = pd.read_csv(profile)
        profile = profile[['adjKey','adjVal','RSS']].groupby(
                ['adjKey','adjVal']).min()
        i = profile.groupby(by="adjKey").ngroups
        nRows = 4
        nCols = 5
        nFrames = math.ceil(i/(nRows*nCols))
        frames = [(j,min(i-j*nRows*nCols,nRows*nCols)) for j in range(nFrames)]
        frames = {j:(math.ceil(k/nCols),min(nCols,k)) for j,k in frames}
        if cmdDict["style"] is None:
            myStyle = "darkgrid"
        else:
            myStyle = cmdDict["style"]
        with sns.axes_style(myStyle):
            frames = {k:list(plt.subplots(v[0], v[1], sharex=True, 
                                          figsize=(12,10))) 
                      for k,v in frames.items()}
            for k in frames:
                if not isinstance(frames[k][1],np.ndarray):
                    frames[k][1] = [frames[k][1]]
                else:
                    frames[k][1] = list(frames[k][1].flat)
            axs = [i for _,v in frames.items() for i in v[1]]
            for (k,v), ax in zip(profile.groupby(by="adjKey"),axs):
                t=v.droplevel(level="adjKey")
                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_title(k)
                ax.plot(t.index, t["RSS"])
            for k,v in frames.items():
                v[0].tight_layout()
                i = resolvePath(["figures", oldCmd["cmdDict"]["dataName"],
                                            "profile"+cmdDict["id"]+"-"+str(k)+
                                            ".png"], relative=True)
                v[0].savefig(i)
