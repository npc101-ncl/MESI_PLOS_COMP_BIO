#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 16:54:43 2022

@author: peter
"""

from python.utilityTools import *
from python.pycotoolsHelpers import *
from GPUtilityFunctions import *
import pandas as pd
import json

def calculateScales(rawDfList,refDfList,toScale):
    """
        must have time as index and must be floats
    """
    myScale = {k:None for k in toScale} 
    joinDfList = [rawDf.join(refDf, lsuffix='_raw', rsuffix='_ref', 
                             how='inner') 
                  for rawDf, refDf in zip(rawDfList, refDfList)]
    for k in myScale:
        numerator = 0
        denominator = 0
        for dfInList in joinDfList:
            if not k+'_raw' in dfInList.columns:
                continue
            if not k+'_ref' in dfInList.columns:
                continue
            numerator += (dfInList[k+'_raw']*dfInList[k+'_ref']).sum()
            denominator += (dfInList[k+'_raw']*dfInList[k+'_raw']).sum()
        if denominator!=0:
            myScale[k]=numerator/denominator
    return myScale
            
def guessScales(rawDfList,dataRelations,forceAverage,sf):
    scale = {k:None for k in dataRelations}
    for ts,ci in dataRelations.items():
        df = [df[[ts,ci]] for df in rawDfList 
              if (ts in df.columns) and (ci in df.columns)]
        if len(df)==0:
            print("insufisent data to guess scale for "+ts+" relative to "+
                  ci+".")
            continue
        if forceAverage:
            for i in range(len(df)):
                df[i][ci] = df[i][ci].mean()
        df = pd.concat(df, ignore_index=True)
        scale[ts] = df[ci].mean()/(sf*df[ts].std()+df[ts].mean()) 
    return scale

def doScales(rawDfList,refDfList,dataRelations,forceAverage=True,sf = 2):
    if refDfList is None:
        return guessScales(rawDfList,dataRelations,forceAverage,sf)
    else:
        return calculateScales(rawDfList,refDfList,dataRelations.keys())
            
def convDfCannonIndex(df):
    if "Time" in df.columns:
        df2 = df.set_index("Time")
    elif "time" in df.columns:
        df2 = df.set_index("time")
    else:
        df2 = df
    df2.index.name = "Time"
    return df2

def clipScales(scales, dfList, dataRelations, scaleMax = None, scaleMin = 0, 
               forceAverage=True):
    newScales = scales    
    if isinstance(scaleMin,dict):
        for k, v in scaleMin.items():
            if k in newScales:
                if newScales[k]<v:
                    newScales[k]=v
    elif scaleMin is not None:
        newScales = {k:max(v,scaleMin) for k,v in newScales.items()}   
    for ts, ci in dataRelations.items():
        for df in dfList:
            if not (ts in df.columns and ci in df.columns and ts in newScales):
                continue
            if forceAverage:
                a = min(df[ci].mean()/df[ts])
            else:
                a = min(df[ci]/df[ts])
            if newScales[ts]>a:
                newScales[ts]=a
    if isinstance(scaleMax,dict):
        for k, v in scaleMax.items():
            if k in newScales:
                if newScales[k]>v:
                    newScales[k]=v
    elif scaleMax is not None:
        newScales = {k:min(v,scaleMax) for k,v in newScales.items()}
    newScales = {k:max(v,0) for k,v in newScales.items()} 
    if any([i==0 for _,i in newScales.items()]):
        print("warning scales for "+
              ", ".join([k for k,i in newScales.items() if i==0])+
              " are 0.")
    return newScales

def getTimesFromAntStr(antStr,pTimeVar,observables):
    myObs = {i[0]+"_"+i[1]:i[0] for i in observables}
    mylines = antStr.splitlines()
    mylines = [i.split("//",1)[0] for i in mylines]
    mylines = [i.split("#",1)[0] for i in mylines]
    mylines = [i.strip() for i in mylines]
    mylines = [i[2:].split(":",1) for i in mylines if i.startswith("at ")]
    mylines = [[i[0],i[1].split(",")] for i in mylines if len(i)==2]
    mylines = [[i[0],[j.split("=") for j in i[1]]] for i in mylines]
    mylines = [[i[0].strip(),[[k.strip() for k in j] for j in i[1] 
                              if len(j)==2]]
               for i in mylines]
    mylines = [[i[0].split(")",1)[0],i[1]] for i in mylines 
               if len(i[0].split(")",1))==2]
    mylines = [[i[0][1:],i[1]] for i in mylines 
               if i[0].startswith("(")]
    mylines = [[i[0],[j for j in i[1] if j[0] in myObs]] for i in mylines]
    mylines = [[i[0],[j[0] for j in i[1] if myObs[j[0]]==j[1]]]
               for i in mylines]
    t = [j for i in mylines for j in i[1]]
    t = [i for i in set(t) if t.count(i)>1]
    if len(t)>0:
        print("Duplicate asinments of observables "+", ".join(t)+
              " in antimony string.")
    mylines = {j:i[0] for i in mylines for j in i[1]}
    for k in mylines:
        if "&&" in mylines[k]:
            mylines[k] = 0.0
            continue
        t = [i for i in ["==",">=","<=",">","<"] if i in mylines[k]]
        if len(t)==0:
            mylines[k] = None
            continue
        t = t[0]
        s = [i.strip() for i in mylines[k].split(t)]
        if len(s)!=2:
            mylines[k] = None
            continue
        if t[0]=="<" or (t[0]=="=" and s[1]==pTimeVar):
            s.reverse()
        if s[0]==pTimeVar:
            mylines[k] = float(s[1])
        else:
            mylines[k] = None
    t = [i for i in myObs if not i in mylines or mylines[i] is None]
    if len(t)>0:
        print("Missing times for observables "+", ".join(t)+".")
    return mylines

def TCToTable(TC,dictKey,indexKey):
    df = TC
    for k in df:
        for i in range(len(df[k])):
            df[k][i][indexKey] = i 
        df[k] = pd.concat(df[k], ignore_index=True)
        df[k][dictKey] = k
    df = pd.concat([v for _,v in df.items()], ignore_index=True)
    return df

def extractFinalValues(timeCourse,observables,times,pTimeVar):
    temp = {k:[i.iloc[-1].to_dict() for i in v] for k,v in timeCourse.items()}
    myKeys = list(temp.keys())
    maxTime = times[max(times)]
    temp = [temp[k] for k in myKeys]
    temp = list(map(list,zip(*temp)))
    temp = [i for i in temp if all([j[pTimeVar]>=maxTime for j in i])]
    temp = list(map(list,zip(*temp)))
    temp = {k:v for k,v in zip(myKeys,temp)}
    temp = {k:[[[observables[obs],val,times[obs]] for obs,val in i.items() 
                 if obs in times]
                for i in v]
               for k,v in temp.items()}
    temp = {k:[pd.DataFrame(i,columns=["cols","vals","Time"]) for i in v]
            for k,v in temp.items()}
    temp = {k:[i.pivot(index='Time', columns='cols', values='vals') for i in v]
            for k,v in temp.items()}
    temp = {k:[i.reset_index() for i in v] for k,v in temp.items()}
    return temp

def dictFromCSV(localPath):
    t = localPath
    if not isinstance(t,list):
        t = [t]
    fn = resolvePath(t, relative=True)
    f = open(fn,"r")
    d = f.read()
    f.close
    if "," in d:
        out = pd.read_csv(resolvePath(t, relative=True), header=None, 
                          index_col = 0).squeeze().to_dict()
    else:
        out = {}
    return out

"""cmdFlags:
    "averageTotals"
"""

cmdDefaults = {"copys":200,
               "preAntStr":"antString2.txt",
               "antStr":"antString2M.txt",
               "parallelAntStr":"",
               "id":"test",
               "data":"zr75_data.xlsx",
               "dataKey":"DTCon.csv",
               "dataRel":"dataRelations.csv",
               "totRel":"totalsRelations.csv",
               "scaleRefs":"PEScaleRefs.csv",
               "paramOut":"parameters.csv",
               "paramRecycle":"",
               "recycleCase":0,
               "recycleVars":"",
               "scaleOut":"scales.csv",
               "scaleIn":"",
               "TCIndex":"index",
               "SRIndexVal":0,
               "TCKeyVar":"dataKey",
               "dataName":"testGPP",
               "TCOutputs":"PETimeCourses.csv",
               "DBOutputs":"DBTimeCourses.csv",
               "TCToComp":10,
               "pTimeVar":"PT",
               "sdForScale":2,
               "scaleMin":0.1,
               "scaleMax":1000,
               "stabAlow":120.0,
               "parameterUB":100.0,
               "parameterLB":0.01,
               "methP":"method:particle_swarm,swarm_size:5,"+
                       "iteration_limit:10",
               "methS":"method:hooke_jeeves",
               "maxRunTime":60*60*47,
               "secondRunTime":60*60*5,
               "twoStage":"no",
               "ratioWeight":1,
               "mergeTotals":False,
               "doNotEstimate":"Insulin:AA",
               "addICforTC":"",
               "jobLimit":1000,
               "PLCase":0,
               "profileOut":"profile.csv",
               "paramTrans":"",
               "paramKey":"paramSet",
               "useForTrans":10,
               "cmdRec":"comandRecord.p",
               "elemOut":"elements.p",
               "timeMul":60.0,
               "copasiPath":"/Applications/copasi"}
alternate = True

if __name__ == "__main__":
    # takes comand line arguments and proceeses them into a usable form
    cmdDict, cmdFlags = getCmdLineArgs()
    cmdDict = parseCmdDict(cmdDict,cmdDefaults)
    for k in ["methP","methS"]:
        if ":" in cmdDict[k]:
            cmdDict[k] = extractCompositArgument(cmdDict[k],sep=",",
                   dictSep=":")
            if not "weight_method" in cmdDict[k]:
                # adds a default weight methiod for the methiod dictionary
                # this setting should end up with all messures having a weight
                # of 1
                cmdDict[k]["weight_method"] = "standard_deviation"
    if "scaleModelKey" in cmdDict:
        cmdDict["scaleModelKey"] = extractCompositArgument(
                cmdDict["scaleModelKey"],sep=",",dictSep=":")
    else:
        cmdDict["scaleModelKey"] = {}
    cmdDict["doNotEstimate"] = extractCompositArgument(
            cmdDict["doNotEstimate"])
    cmdDict["addICforTC"] = extractCompositArgument(
            cmdDict["addICforTC"])
    cmdDict["recycleVars"] = extractCompositArgument(
            cmdDict["recycleVars"])
    if "PLRange" in cmdDict:
        cmdDict["PLRange"] = extractCompositArgument(
                cmdDict["PLRange"])
        cmdDict["PLRange"] = [float(i) for i in cmdDict["PLRange"]]
        
    if "configFile" in cmdDict:
        f = open(cmdDict["configFile"],'r')
        myJ = f.read()
        f.close()
        myJ = json.loads(myJ)
        if "cmdDict" in myJ:
            cmdDict.update(myJ["cmdDict"])
        if "cmdFlags" in myJ:
            cmdFlags = set([k for k,v in myJ["cmdFlags"].items() 
                            if v]).union(cmdFlags)
            cmdFlags = cmdFlags - set([k for k,v in myJ["cmdFlags"].items() 
                                       if not v])
    
    mySuperComputer = "slurm" in cmdFlags
    twoStage = cmdDict["twoStage"].lower() == "yes"
    
    if len(cmdDict["paramRecycle"])>0:
        recycleParams = pd.read_csv(resolvePath(["data",cmdDict["dataName"],
                                                 cmdDict["paramRecycle"]],
                                                relative=True))
        if len([i for i in cmdDict["recycleVars"]
                if not i in recycleParams.columns])>0:
            print("not in recycled params",
                  [i for i in cmdDict["recycleVars"]
                   if not i in recycleParams.columns])
        cmdDict["recycleVars"] = [i for i in cmdDict["recycleVars"]
                                  if i in recycleParams.columns]
        recycleParams = recycleParams[cmdDict["recycleVars"]]
        recycleParams = recycleParams.iloc[cmdDict["recycleCase"]].to_dict()
    else:
        recycleParams = None
    
    # set a time to end the parameter estimation early if on clustor
    myEndTime = getEndTime(seconds = cmdDict["maxRunTime"])
    if twoStage:
        myMidTime = getEndTime(
                seconds = cmdDict["maxRunTime"]-cmdDict["secondRunTime"])
    else:
        myMidTime = myEndTime
    
    # if on local machine add path to copasiSE and force 2 copies for testing
    if not mySuperComputer:
        cmdDict["copys"]=2
        addCopasiPath(cmdDict["copasiPath"])
        
    # save comands used for run to file for refrence
    if not (("PLRange" in cmdDict) and len(cmdDict["PLRange"])>0):
        savePick(["data",cmdDict["dataName"],cmdDict["cmdRec"]], 
                 {"cmdDict":cmdDict, "cmdFlags":cmdFlags}, relative=True)
    
    # extract information about from the antimony string for the base model   
    preModel = modelRunner(antString=loadTxt([cmdDict["preAntStr"]],
                                             relative=True), 
                           run_dir=resolvePath(["copasiRuns",cmdDict["id"]],
                                               relative=True))
    preModel.clearRunDirectory()
        
    origEliments = preModel.getModelEliments()
    toEstimate = [i for i in origEliments["kineticParams"] 
                  if (not i in cmdDict["doNotEstimate"])
                  and (not i in cmdDict["recycleVars"])]
    observables = origEliments["assignments"]+origEliments["metabolites"]
    
    # extract data form the antimony string that has been modified with psudo
    # time and analyticaly set initial conditions
    myModel = modelRunner(antString=loadTxt([cmdDict["antStr"]],
                                            relative=True), 
                          run_dir=resolvePath(["copasiRuns",cmdDict["id"]],
                                              relative=True))
    
    myModel.clearRunDirectory()
    
    eliments = myModel.getModelEliments()
    t = [i.split("_") for i in eliments["kineticParamsGlobal"]]
    t = [["_".join(i[:-1]),i[-1]] for i in t
         if len(i)>1 and i[-1].startswith("DVT")]
    observables = [i for i in t if i[0] in observables]
    
    t = {"pre elements":origEliments, "post elements":eliments}
    
    if len(cmdDict["parallelAntStr"])>0:
        myModel = modelRunner(antString=loadTxt([cmdDict["parallelAntStr"]],
                                                relative=True), 
                              run_dir=resolvePath(["copasiRuns",cmdDict["id"]],
                                                  relative=True))
        t["parallel elements"] = myModel.getModelEliments()
        myModel.clearRunDirectory()
    if not (("PLRange" in cmdDict) and len(cmdDict["PLRange"])>0):
        savePick(["data",cmdDict["dataName"],cmdDict["elemOut"]], t, 
                 relative=True)
    
    # import data and the specifications for the data's significance
    # eg what is a component of what, what relates to a total condition in the
    # antimony string
    df = importData([cmdDict["data"]],timeMultiplyer=cmdDict["timeMul"])
    if "configFile" in cmdDict and "dataKey" in myJ:
        dataKey = myJ["dataKey"]
    else:
        dataKey = dictFromCSV(cmdDict["dataKey"])
    if "configFile" in cmdDict and "dataRel" in myJ:
        dataRelations = myJ["dataRel"]
    else:
        dataRelations = dictFromCSV(cmdDict["dataRel"])
    if "configFile" in cmdDict and "totRel" in myJ:
        totalRelations = myJ["totRel"]
    else:
        totalRelations = dictFromCSV(cmdDict["totRel"])
    if len(cmdDict["scaleIn"])>0:
        t = resolvePath(["data",cmdDict["dataName"],cmdDict["scaleIn"]],
                        relative=True)
        scaleOveride = pd.read_csv(t)
    else:
        scaleOveride = None
    
    """
    dataKey = pd.read_csv(resolvePath([cmdDict["dataKey"]], relative=True),
                          header=None, index_col = 0).squeeze().to_dict()
    dataRelations = pd.read_csv(resolvePath([cmdDict["dataRel"]],
                                            relative=True), 
                                header=None, index_col = 0).squeeze().to_dict()
    totalRelations = pd.read_csv(resolvePath([cmdDict["totRel"]],
                                             relative=True), 
                                header=None, index_col = 0).squeeze().to_dict()
    """
    
    for k in df:
        df[k] = df[k].drop(columns=[i for i in df[k].columns
                                    if not i in dataKey])
        df[k] = df[k].rename(columns=dataKey)
    
    if len(cmdDict["parallelAntStr"])==0:
        scaleRefPath = resolvePath(["data",cmdDict["dataName"],
                                    cmdDict["scaleRefs"]], relative=True)
        try:
            scaleRef = pd.read_csv(scaleRefPath, index_col = False)
        except FileNotFoundError:
            scaleRef = None
        if scaleRef is not None:
            scaleRef = tableToDict(scaleRef, indexKey = cmdDict["TCIndex"],
                                   indexVal = cmdDict["SRIndexVal"],
                                   keyVar = cmdDict["TCKeyVar"])
            scaleRef = {k:convDfCannonIndex(df) 
                        for k, df in scaleRef.items()}
            scaleRef = [[v,scaleRef[k]] for k,v in df.items()]
            scaleRef = list(map(list, zip(*scaleRef)))
            dfList = scaleRef[0]
            scaleRef = scaleRef[1]
        else:
            dfList = [i for _,i in df.items()]
            
        scales = doScales(dfList, scaleRef, dataRelations,
                          forceAverage = "averageTotals" in cmdFlags,
                          sf = cmdDict["sdForScale"])
        
        scales = clipScales(scales, dfList, dataRelations,
                            scaleMax = cmdDict["scaleMax"], 
                            scaleMin = cmdDict["scaleMin"], 
                            forceAverage = "averageTotals" in cmdFlags)
        
        df = scale(df,scales)
        totals = {k:{t:v[col].mean() for col,t in totalRelations.items() 
                 if col in v.columns} for k,v in df.items()}
    else:
        totals = {k:{t:v[col].mean() for col,t in totalRelations.items() 
                 if col in v.columns} for k,v in df.items()}
        if scaleOveride is not None:
            noScales = [i for i in scaleOveride.columns 
                        if i in cmdDict["scaleModelKey"]]
            scaleOveride = scaleOveride.set_index("index")
            scaleOveride = scaleOveride[list(cmdDict["scaleModelKey"])]
            scaleOveride = scaleOveride.rename(
                    columns=cmdDict["scaleModelKey"])
        else:
            noScales = [] 
        df, dfMeans = normalize(df,acrossLines = alternate, noScale=noScales)
        
    times = getTimesFromAntStr(loadTxt([cmdDict["antStr"]], relative=True),
                               cmdDict["pTimeVar"],observables)
    
    cellLines = {i.replace("-","").replace("_",""):i for i in df.keys()}
    
    # code to get paralel observables mapping paralel obsevables to post
    # ant string observables and cell line
    
    # will need a new format for observables to incoperate cell line as a 
    # factor
    maxTime = times[max(times)]
    
    if len(cmdDict["parallelAntStr"])==0:
        observables = {i[0]+"_"+i[1]:i[0] for i in observables}
        
        # need new code here for new paralised data format
        PECSV = {lk:{k:lv.at[v,observables[k]] for k,v in times.items()}
                 for lk,lv in df.items()}
        for k in PECSV:
            PECSV[k]["Time"]=maxTime+cmdDict["stabAlow"]
        PECSV = {k:pd.DataFrame([v]) for k,v in PECSV.items()}
        
        PECSV, indepV = tuple(map(list, 
                                  zip(*[[v,totals[k]] for k,v 
                                        in PECSV.items()])))
        
        CSVPath = []
        for i, myCSV in enumerate(PECSV):
            CSVPath.append(resolvePath(["data",cmdDict["dataName"],
                                        cmdDict["id"]+"PEInput"+str(i)+".csv"],
                                       relative=True))
            myCSV.to_csv(CSVPath[-1])
    else:
        times = {i+"_"+j:v for i,v in times.items() for j in cellLines.keys()}
        observables = {i[0]+"_"+i[1]+"_"+j:[i[0],v] for i in observables 
                       for j,v in cellLines.items()}
        PECSV = {k:df[observables[k][1]].at[v,observables[k][0]] 
                 for k,v in times.items()}
        if not alternate:
            t = set.union(*[set(d.keys()) for cl,d in dfMeans.items()])
            t = {obs:sum([d[obs] for _,d in dfMeans.items() if obs in d])
                 for obs in t}
            PECSV.update({obs+"_"+cl.replace("_","").replace("-","")+"_ratio":
                          cmdDict["ratioWeight"]*v/t[obs] 
                          for cl,d in dfMeans.items() for obs,v in d.items()})
        PECSV["Time"]=maxTime+cmdDict["stabAlow"]+2
        PECSV = pd.DataFrame([PECSV])
        indepV = {i+"_"+k:j for k,v in cellLines.items() 
                  for i,j in totals[v].items()}
        if not alternate:
            indepV.update({"rWeight":cmdDict["ratioWeight"]})
        tNames = set.union(*[set(v.keys()) for _,v in totals.items()])
        tNames = [i for i in eliments["kineticParams"] if i.endswith("_T") 
                  and not i in tNames]
        indepV = [indepV]
        if cmdDict["mergeTotals"]:
            toEstimate = toEstimate + [i for i in tNames 
                                       if i not in cmdDict["recycleVars"]]
        else:
            toEstimate = toEstimate + [i+"_"+j for i in tNames 
                                       for j in cellLines.keys()
                                       if i+"_"+j not in 
                                       cmdDict["recycleVars"]]
        CSVPath = resolvePath(["data",cmdDict["dataName"],
                               cmdDict["id"]+"PEInput.csv"], relative=True)
        print(PECSV)
        PECSV.to_csv(CSVPath)
        CSVPath = [CSVPath]
     
    print("CSVPath:",CSVPath)
    print("copyNum:",cmdDict["copys"])
    print("rocket:",mySuperComputer)
    print("estimatedVar:",toEstimate)
    print("upperParamBound:",cmdDict["parameterUB"])
    print("lowerParamBound:",cmdDict["parameterLB"])
    print("method:",cmdDict["methP"])
    print("endTime:",myMidTime)
    print("indepToAdd:",indepV)
    if recycleParams is not None:
        print("recycleParams:",recycleParams)
    
if (__name__ == "__main__" and ("PLRange" in cmdDict) and 
   len(cmdDict["PLRange"])>0):
    params = pd.read_csv(resolvePath(["data",cmdDict["dataName"],
                                      cmdDict["paramOut"]],
                                     relative=True),
                         index_col=0)
    params = params.iloc[cmdDict["PLCase"]].to_dict()
    params = {k:v for k,v in params.items() if k!="RSS"}
    toEstimate = [i for i in toEstimate if i in params]
    if "userName" in cmdDict:
        myThrottle={"user":cmdDict["userName"], "jobLimit":cmdDict["jobLimit"]}
    else:
        myThrottle = None
    params = myModel.runProfileLikelyhood(CSVPath,cmdDict["PLRange"],
                                          toEstimate,
                                          rocket=mySuperComputer,
                                          upperParamBound=
                                          cmdDict["parameterUB"],
                                          lowerParamBound=
                                          cmdDict["parameterLB"],
                                          method=cmdDict["methS"],
                                          overrideParam=params,
                                          indepToAdd=indepV,
                                          depth=cmdDict["copys"],
                                          throttle=myThrottle)
    i,j,_ = params
    for k in set.intersection(set(i.keys()),set(j.keys())):
        i[k]["adjKey"]=k
        i[k]["adjVal"]=i[k][k]/j[k]
    params = pd.concat([j for _,j in i.items()],ignore_index=True)
    params.to_csv(resolvePath(["data",cmdDict["dataName"],
                               cmdDict["profileOut"]], relative=True))
    myModel.clearRunDirectory()
elif __name__ == "__main__" and len(cmdDict["paramTrans"])>0:
    paramsOR = pd.read_csv(resolvePath(["data",cmdDict["dataName"],
                                        cmdDict["paramTrans"]],
                                     relative=True),
                         index_col=0)
    if "useForTrans" in cmdDict:
        if isinstance(cmdDict["useForTrans"],int):
            cmdDict["useForTrans"] = [i for i in range(cmdDict["useForTrans"])]
        elif not isinstance(cmdDict["useForTrans"],list):
            cmdDict["useForTrans"] = [i for i in range(10)]
    else:
        cmdDict["useForTrans"] = [i for i in range(10)]
    paramsOR = paramsOR.iloc[cmdDict["useForTrans"]]
    i = [i for i in paramsOR.columns if not i in toEstimate]
    paramsOR = paramsOR.drop(columns=i)
    toEstimate = [i for i in toEstimate if not i in paramsOR.columns]
    # insert code for scales
    if (scaleOveride is not None) and len(scaleOveride.columns)>0:
        params = pd.concat([paramsOR,scaleOveride],axis=1)
    else:
        params = paramsOR
    if len(toEstimate)>0:
        passCols = toEstimate
        if params is not None:
            passCols = list(set(passCols).union(set(params.columns)))
        params = myModel.runParamiterEstimation(CSVPath,
                                                copyNum=cmdDict["copys"],
                                                rocket=mySuperComputer,
                                                estimatedVar=toEstimate,
                                                upperParamBound=
                                                cmdDict["parameterUB"],
                                                lowerParamBound=
                                                cmdDict["parameterLB"],
                                                method=cmdDict["methP"],
                                                overrideParam=params,
                                                endTime=myMidTime,
                                                indepToAdd=indepV)
        if twoStage:
            print("stage 2 start")
            params = GFID(params)[passCols]
            params = myModel.runParamiterEstimation(CSVPath,copyNum=1,
                                                    rocket=mySuperComputer,
                                                    estimatedVar=toEstimate,
                                                    upperParamBound=
                                                    cmdDict["parameterUB"],
                                                    lowerParamBound=
                                                    cmdDict["parameterLB"],
                                                    method=cmdDict["methS"],
                                                    overrideParam=params,
                                                    endTime=myEndTime,
                                                    indepToAdd=indepV,
                                                    randStartVal=False)
    else:
        params = [r for _,r in params.iterrows() 
                  for _ in range(cmdDict["copys"])]
        params = {"noParam":pd.DataFrame(params)}
        params["noParam"]["RSS"] = None
    for k in params:
        params[k][cmdDict["paramKey"]] = [i for i in cmdDict["useForTrans"]
                                          for _ in range(cmdDict["copys"])]
        for i in paramsOR.columns:
            if not i in params[k].columns and i!="RSS":
                params[k][i] = [j for _,j in paramsOR[i].items()
                                for _ in range(cmdDict["copys"])]
    params = {k:v.sort_values(by="RSS").reset_index(drop=True)
              for k,v in params.items()}
elif __name__ == "__main__":
    passCols = toEstimate
    if recycleParams is not None:
        passCols = list(set(passCols).union(set(recycleParams.keys())))
    params = myModel.runParamiterEstimation(CSVPath,copyNum=cmdDict["copys"],
                                            rocket=mySuperComputer,
                                            estimatedVar=toEstimate,
                                            upperParamBound=
                                            cmdDict["parameterUB"],
                                            lowerParamBound=
                                            cmdDict["parameterLB"],
                                            method=cmdDict["methP"],
                                            overrideParam=recycleParams,
                                            endTime=myMidTime,
                                            indepToAdd=indepV)
    
    if twoStage:
        print("stage 2 start")
        params = GFID(params)[toEstimate]
        if recycleParams is not None:
            for k,v in recycleParams.items():
                params[k] = v
        params = myModel.runParamiterEstimation(CSVPath,copyNum=1,
                                                rocket=mySuperComputer,
                                                estimatedVar=toEstimate,
                                                upperParamBound=
                                                cmdDict["parameterUB"],
                                                lowerParamBound=
                                                cmdDict["parameterLB"],
                                                method=cmdDict["methS"],
                                                overrideParam=params,
                                                endTime=myEndTime,
                                                indepToAdd=indepV,
                                                randStartVal=False)
        if recycleParams is not None:
            for j in params:
                for k,v in recycleParams.items():
                    if k in params[j].columns:
                        params[j][k]=v
        
    params = {k:v.sort_values(by="RSS").reset_index(drop=True)
              for k,v in params.items()}
if __name__ == "__main__" and not (("PLRange" in cmdDict) and 
                                   len(cmdDict["PLRange"])>0):
    GFID(params).to_csv(resolvePath(["data",cmdDict["dataName"],
                                     cmdDict["paramOut"]],
                                    relative=True))
    if len(cmdDict["paramTrans"])>0:
        for k in params:
            params[k] = params[k].drop(columns=[cmdDict["paramKey"]])
    if len(cmdDict["parallelAntStr"])==0:
        pd.DataFrame([scales]).to_csv(resolvePath(["data",cmdDict["dataName"],
                                                   cmdDict["scaleOut"]],
                                                  relative=True))
    
    if len(cmdDict["parallelAntStr"])!=0:
        maxTime = maxTime+2
        myModel.clearRunDirectory()
        myCond = GFID(params)
        passCols = [i for i in passCols if i in myCond.columns]
        myCond = myCond.iloc[:cmdDict["TCToComp"]][passCols].copy()
        for k,v in indepV[0].items():
            myCond[k]=v
        print("DB timecourse Start")
        timeCourse = myModel.runTimeCourse(maxTime+cmdDict["stabAlow"],
                                              adjustParams=myCond,
                                              stepSize=10)
        DBPath = resolvePath(["data",cmdDict["dataName"],cmdDict["DBOutputs"]],
                             relative=True)
        timeCourse = TCToTable({"psudo":timeCourse},cmdDict["TCKeyVar"],
                               cmdDict["TCIndex"])
        timeCourse.to_csv(DBPath)
        icDict = [i for i,j in timeCourse.nunique().items() if j==1]
        icDict = [i for i in icDict if (("_" in i) and 
                  (i.split("_",1)[-1] in origEliments["kineticParams"]))]
        icDict = timeCourse[icDict].iloc[0].to_dict()
    
    timeCourse = {}      
    myModel.clearRunDirectory()
    myModel = modelRunner(antString=loadTxt([cmdDict["antStr"]],
                                            relative=True),
                          run_dir=resolvePath(["copasiRuns",cmdDict["id"]],
                                              relative=True))
    if len(cmdDict["parallelAntStr"])==0:
        for k,myDict in totals.items():
            myCond = GFID(params)
            myCond = myCond.iloc[:cmdDict["TCToComp"]][toEstimate].copy()
            for t,v in myDict.items():
                myCond[t] = v
            print(k,"timecourse Start")
            timeCourse[k] = myModel.runTimeCourse(maxTime+cmdDict["stabAlow"],
                                                  adjustParams=myCond,
                                                  stepSize=10)
    else:
        timeCourse = {} 
        for k,v in cellLines.items():
            myCond = GFID(params).iloc[:cmdDict["TCToComp"]].copy()
            t = ["_"+i for i in (set(cellLines.keys())-set([k]))]
            t = [i for i in myCond.columns 
                 if not any([i.endswith(j) for j in t])]
            myCond = myCond[t]
            myCond = myCond.rename(columns={i:i[:-(len(k)+1)] for i in t 
                                            if i.endswith("_"+k)})
            t = {i.split("_",1)[1]:j for i,j in icDict.items() 
                 if k==i.split("_",1)[0]}
            for i,j in t.items():
                if not i in myCond.columns and i in cmdDict["addICforTC"]:
                    myCond[i]=j
            myCond = myCond.drop(columns=["RSS","rWeight"], errors="ignore")
            for i,j in totals[v].items():
                myCond[i]=j
            print(myCond.iloc[0].to_dict())
            print(k,"timecourse Start")
            timeCourse[k] = myModel.runTimeCourse(maxTime+cmdDict["stabAlow"],
                                                  adjustParams=myCond,
                                                  stepSize=10)
            
    if len(cmdDict["parallelAntStr"])==0:
        refTC = extractFinalValues(timeCourse,observables,times,
                                   cmdDict["pTimeVar"])
        refTC = TCToTable(refTC,cmdDict["TCKeyVar"],cmdDict["TCIndex"])
        refTC.to_csv(scaleRefPath)
    
    TCPath = resolvePath(["data",cmdDict["dataName"],cmdDict["TCOutputs"]],
                         relative=True)
    timeCourse = TCToTable(timeCourse,cmdDict["TCKeyVar"],cmdDict["TCIndex"])
    timeCourse.to_csv(TCPath)


