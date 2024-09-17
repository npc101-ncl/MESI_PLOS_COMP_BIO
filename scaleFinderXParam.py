#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 16:47:06 2020

@author: peter
"""

import site, os, re, pickle, time, sys
import pandas as pd
import numpy as np
from python.pycotoolsHelpers import *
from python.analysisTools import *
from python.visualisationTools import *
from python.utilityTools import *
from scipy.interpolate import interp1d
from scaleFinderXPrimer import makeInterpolatedTable, applyScale
from scaleFinderXPrimer import importData, genScale, genConstShift
from scaleFinderXPrimer import genRunIn, duplicator

cmdLineArg = sys.argv[1:]
CLDict, CLSet = getCmdLineArgs()

# import switches data and pramiters from primer / comand line arguments
working_directory = os.path.dirname(os.path.abspath(__file__))
if len(cmdLineArg)>0:
    data_dir = os.path.join(working_directory,'data', cmdLineArg[0])
else:
    data_dir = os.path.join(working_directory,'data', 'noninterpolatedX')

if "meth" in CLDict.keys():
    set_methiod = CLDict["meth"]
    set_methiod = set_methiod.split(sep=",")
    set_methiod = [ent.split(sep=":") for ent in set_methiod]
    set_methiod = {ent[0]:ent[1] for ent in set_methiod if len(ent)==2}
    for key in set_methiod.keys():
        if key!="method":
            set_methiod[key]=int(set_methiod[key])
else:
    allowable_methiods = ["particle_swarm_test", "particle_swarm_default", 
                          "particle_swarm_heroic"]
    set_methiod = [meth for meth in allowable_methiods if meth
                   in cmdLineArg]
    if len(set_methiod)==0:
        set_methiod = allowable_methiods[-1]
    else:
        set_methiod = set_methiod[-1]

isZr75 = "ZR75" in CLSet
isTransfer = "trans" in CLSet
isCombo = "combo" in CLSet
isS6K = "S6K" in CLSet
stabilizeIC = "ICS" in CLSet
noSSFlag = "noSS" in CLSet

if "scaleMin" in CLDict.keys():
    scaleMin=CLDict["scaleMin"]
else:
    scaleMin = 0.01

if "mTorConst" in CLSet:
    optMTORC = {"MTORC1_wb":1,"MTORC2_wb":1}
else:
    optMTORC = None

if "pi3kConst" in CLSet:
    optMTORC["PI3K_sum"] = 1

file = open(os.path.join(data_dir,'runSwitches.p'),'rb')
RS = pickle.load(file)
file.close()

if isZr75:
    data_filename = RS["data_filenames"]["zr75"]
else:
    data_filename = RS["data_filenames"]["t47d"]

secondsToRun = 60*60*47
secondsToRun2 = 60*60*45
endTime = time.time()+secondsToRun 
endTime2 = time.time()+secondsToRun2 

myUpperBound=10
mySuperComputer = "slurm" in CLSet

prequilFlag = "pre" in CLSet
    
if mySuperComputer: 
    myCopyNum=200
else:
    myCopyNum=2
    
if not mySuperComputer:
    addCopasiPath("/Applications/copasi")
    
kineticParams = RS["kineticParams"]

kineticBGRParams = RS["kineticBGRParams"]

kineticParamsOpt = RS["kineticParamsOpt"]

overrideDict = RS["overrideDict"]
    
if __name__ == "__main__":
    
    # set up paramiterisation based on switches
    if isS6K:
        kineticParams.extend(kineticParamsOpt)
        runName=RS["run_name_prefix"]+"-S6K"
    else:
        runName=RS["run_name_prefix"]
    if isZr75:
        runName=runName+"-ZR75"
    else:
        runName=runName+"-T47D"
    if "IRS1SynthParam" in RS and "IRSSynth" in CLSet:
        kineticParams.extend(RS["IRS1SynthParam"])
    if "isFullBGR" in CLSet:
        kineticParams.extend(kineticBGRParams)
    
    # create model object to run tasks on model
    if isCombo and (not isZr75):
        myCase = 'MCF7T47D'
    elif isCombo and isZr75:
        myCase = 'MCF7ZR75'
    elif isTransfer and (not isZr75):
        myCase = 'T47D'
    elif isTransfer and isZr75:
        myCase = 'ZR75'
    else:
        myCase = 'MCF7'
    cpfilep = os.path.join(working_directory,'copasiRuns',
                           runName+'-'+myCase)
    myModel = modelRunner(RS["antimony_string"], cpfilep)
    
    RS["modelEliments"] = myModel.getModelEliments()
    columns_from_data = {k:v for k,v
                         in RS["columns_from_data"].items() 
                         if (v in RS["modelEliments"]["metabolites"]) 
                         or (v in RS["modelEliments"]["assignments"])}
    RS["columns_from_data_in_model"] = columns_from_data
    data_relations = {k:v for k,v in RS["data_relations"].items()
                      if (k in RS["modelEliments"]["metabolites"])
                      or (k in RS["modelEliments"]["assignments"])}
    data_relations = {k:[i for i in v
                         if (i in RS["modelEliments"]["metabolites"]) 
                         or (i in RS["modelEliments"]["assignments"])]
                      for k,v in data_relations.items()}
    RS["data_relations_in_model"] = data_relations
    
    metabolites = RS["metabolites"]
    
    if RS["autoEstVar"]:
        kineticParams = RS["modelEliments"]["kineticParams"]        
        metabolites = RS["modelEliments"]["metabolites"]
        
    if isCombo:
        metabolites = (["ModA_"+i for i in metabolites]+
                       ["ModB_"+i for i in metabolites])
        RS["comb_metabolites"] = metabolites
    
    if "noMet" in CLSet:
        pass
    else:
        kineticParams.extend(metabolites)
    if "estimationBlockPattern" in RS:
        tempA = [i for i in kineticParams
                 if i.endswith(RS["estimationBlockPattern"])]
    else:
        tempA = None
    if "estimationForcePattern" in RS:
        tempB = [i for i in kineticParams
                 if i.endswith(RS["estimationForcePattern"])]
    else:
        tempB = None
    if "estimationRegExPattern" in RS:
        tempC = [i for i in kineticParams for j in RS["estimationRegExPattern"]
                 if re.search(j, i) is not None]
    else:
        tempC = None
    if tempC is not None:
        kineticParams = list(set(tempC))
    if tempA is not None:
        kineticParams = list(set(kineticParams)-set(tempA))
    if tempB is not None:
        kineticParams = list(set(kineticParams)+set(tempB))
    if "estimationBlockList" in RS:
        kineticParams = list(set(kineticParams)-
                             set(RS["estimationBlockList"]))
    if "estimationForceList" in RS:
        kineticParams = list(set(kineticParams)+
                             set(RS["estimationForceList"]))
    
    # importing the MCF7 data set
    df = importData('MCF-7',data_filename,columns_from_data,
                    addConstColls=optMTORC)
    dfSD = importData('MCF-7',data_filename,columns_from_data,
                      usecols="X:AS")
    dfSD = dfSD*2 # get standard deviation of data
    
    # importing the T47D (or ZR75) data set to ensure scalings remain
    # physicaly plausable
    if isZr75:
        df2 = importData('ZR-75-1',data_filename,columns_from_data,
                         addConstColls=optMTORC)
    else:
        df2 = importData('T47D',data_filename,columns_from_data,
                         addConstColls=optMTORC)
    if not (isTransfer or isCombo):
        df2 = pd.concat([df, df2],ignore_index=True)
    
    # optional removal of TSC2 data from the data set
    if RS["suppressTSC2Data"]:
        df=df.drop(['TSC2_wb','TSC2pT1462_wb'], axis=1, errors='ignore')
        df2=df2.drop(['TSC2_wb','TSC2pT1462_wb'], axis=1, errors='ignore')
    
    if isCombo:
        df = [df,df2]
        df2 = None
    
    # optional interpolation of experamental data
    if RS["is_interpolated"]:
        df = makeInterpolatedTable(df)
    
    # opional shifting of MCF7 data so initial values start at 0.1 and 1   
    if RS["is_data_constant_shifted"] and (not isTransfer):
        df, df2 = genConstShift(df, df2, data_relations)
    
    # generate scalings if posable from previous run data
    if isCombo:
        prefix = "Combo-"
    else:
        prefix = ""
    if prefix+"recycle-"+runName in RS.keys():
        targets = RS[prefix+"recycle-"+runName]
        if prefix+"recycleSS-"+runName in RS.keys():
            targetsSS = RS[prefix+"recycleSS-"+runName]
        else:
            targetsSS = {}
        myScale = genScale(df2,df,data_relations, scaleTarget=targets,
                           initialTarget=targetsSS)
    else:
        myScale = genScale(df2,df,data_relations)
    if isTransfer:
        myScale = RS["scale-"+runName]
    myScale = {k:max(scaleMin,v) for k,v in myScale.items()}
    df = applyScale(df, myScale)
    if isTransfer:
        df = applyScale(df2, myScale)
    dfSD = applyScale(dfSD, myScale) # get scaled standard deviation of data
    
    # generate list of initial conditions to assume from data
    overrideFilter = RS["modelEliments"]["kineticParams"]+metabolites
    if isTransfer and (not isCombo):
        overrideDictB = loadPick(["data",cmdLineArg[0],
                                  runName+'-MCF7-params.p'],
            relative=True)
        overrideDictB = GFID(overrideDictB)
        overrideDictB = overrideDictB[[col for col in overrideDictB.columns
                                       if (col!="RSS")
                                       and (not col in [met+"_F" for met
                                                        in metabolites])
                                       and (not col in [met+"_z" for met
                                                        in metabolites])]]
        overrideDictB = overrideDictB.iloc[:3]
    elif isCombo:
        if "initialsByMean" in CLSet:
            dfORT = [i.mean(axis=0) for i in df]
        else:
            dfORT = [i.iloc[0] for i in df]
        overrideDictB = ((dfORT[0])[list(overrideDict)].rename(
                index={k:"ModA_"+v for k,v 
                       in overrideDict.items()})).to_dict()
        overrideDictB.update(((dfORT[1])[list(overrideDict)].rename(
                index={k:"ModB_"+v for k,v 
                       in overrideDict.items()})).to_dict())
        overrideFilter = (overrideFilter + 
                          ["ModA_"+i for i in overrideFilter 
                           if i.endswith("_T")]+
                          ["ModB_"+i for i in overrideFilter 
                           if i.endswith("_T")])
    else:
        if "initialsByMean" in CLSet:
            dfORT = df.mean(axis=0)
        else:
            dfORT = df.iloc[0]
        overrideDictB = (dfORT[list(overrideDict)].rename(
                index=overrideDict)).to_dict()
    if isinstance(overrideDictB,dict):
        print(overrideDictB)
        overrideDictB = {k:v for k,v in overrideDictB.items()
                         if (k in overrideFilter)}
        print(overrideDictB)
        print(overrideFilter)
    elif isinstance(overrideDictB,pd.DataFrame):
        overrideDictB = overrideDictB[[k for k in overrideDictB.columns
                        if (k in overrideFilter)]]
    
    # remove given initial conditions from list of paramiters to estimate
    kineticParams = [i for i in kineticParams if i not in list(overrideDictB)]
    
    
    # clear the run directory of junk (don't keep anything permenent in
    # there)
    myModel.clearRunDirectory()
    
    # generate file of scaled dynamics data
    if isCombo and (not isZr75):
        fnameI = "MCF7T47D"
    elif isCombo and isZr75:
        fnameI = "MCF7ZR75"
    elif isTransfer and (not isZr75):
        fnameI = "T47D"
    elif isTransfer and isZr75:
        fnameI = "ZR75"
    else:
        fnameI = "MCF7"
    fname = os.path.join(working_directory,
                         runName+'_'+fnameI+'_exp_data.csv')
    if isCombo:
        df[0].rename(columns = {i:"ModA_"+i for i in df[0].columns 
          if not i in ["time","Time"]}, inplace=True)
        df[1].rename(columns = {i:"ModB_"+i for i in df[1].columns 
          if not i in ["time","Time"]}, inplace=True)
        df = df[0].join(df[1], on='Time', how="inner")
    df.to_csv(fname)
    
    if prequilFlag:
        if not isinstance(df,list):
            df = [df]
        fname4 = []
        timeOffset = []
        for i, df2 in enumerate(df):
            timeOffsetT, fnameT = genRunIn(df2,working_directory,
                                           runName,i)
            timeOffset.append(timeOffsetT)
            fname4.append(fnameT)
           
    
    # generate file of scaled initial unpeterbed steady state data
    fname2 = os.path.join(working_directory,
                              runName+'_'+fnameI+'_exp_S_data.csv')
    df[df.index == 0].to_csv(fname2, index=False)

    if stabilizeIC:
        df=df[[]]
        for col in [met+"_F" for met in metabolites]:
            df[col] = 0
        fname3 = os.path.join(working_directory,
                              runName+'_'+fnameI+'_exp_ICs_data.csv')
        df.to_csv(fname3)    
    
    # set independent variables associated with dynamic and steady state
    # data
    if prequilFlag:
        indepVar = [{"zeroPeriod":TOff} for TOff in timeOffset]
        fname = fname4
    elif stabilizeIC:
        indepVar = [{},{"AA":0, "Insulin":0},{"AA":0, "Insulin":0}]
        fname = [fname,fname2,fname3]
    elif noSSFlag:
        indepVar = [{}]
        fname = [fname]
    else:
        indepVar = [{},{"AA":0, "Insulin":0}]
        fname = [fname,fname2]
        
    if isCombo:
        myModel.clearRunDirectory()
        newAntStr = duplicator(RS["antimony_string"],cpfilep,
                               paramToJoin = list(set([i for i 
                                                       in kineticParams
                                                       if i.startswith("k")]+
                                                      ["AA","Insulin"])))
        print("\n"+newAntStr+"\n")
        myModel.clearRunDirectory()
        myModel = modelRunner(newAntStr, cpfilep)
    
    if isTransfer:
        kineticParams = metabolites
    for k,v in ({"fname":fname, "copyNum":myCopyNum,
                 "rocket":mySuperComputer, "estimatedVar":kineticParams, 
                 "upperParamBound":myUpperBound, "method":set_methiod, 
                 "overrideParam":overrideDictB, "endTime":endTime2, 
                 "indepToAdd":indepVar}).items():
        # gets sent to slurm file can be useful
        print(k,":",v)
    # do paramiter estimation
    params = myModel.runParamiterEstimation(fname,copyNum=myCopyNum,
                                            rocket=mySuperComputer,
                                            estimatedVar=kineticParams,
                                            upperParamBound=myUpperBound,
                                            method=set_methiod,
                                            overrideParam=overrideDictB,
                                            endTime=endTime2,
                                            indepToAdd=indepVar)
    
    if "globalChase" in CLSet:
        params = GFID(params)
        params.to_csv(os.path.join(data_dir,runName+'-MCF7-test1.csv'))
        params = params[[col for col in params.columns
                         if (col!="RSS")
                         and (not col in [met+"_F" for met in metabolites])
                         and (not col in [met+"_z" for met in metabolites])]]
        params.to_csv(os.path.join(data_dir,runName+'-MCF7-test2.csv'))

        params = myModel.runParamiterEstimation(fname,copyNum=1,
                                                rocket=mySuperComputer,
                                                estimatedVar=kineticParams,
                                                upperParamBound=
                                                myUpperBound,
                                                method=
                                                {'method':'hooke_jeeves'},
                                                overrideParam=params,
                                                endTime=endTime,
                                                indepToAdd=indepVar,
                                                randStartVal=False)
        for key in params.keys():
            params[key].sort_values(by=['RSS'], inplace=True)
            params[key].reset_index(drop=True, inplace=True)
            print(params[key])
        GFID(params).to_csv(os.path.join(data_dir,runName+'-MCF7-test3.csv'))
        print("2nd PE ends")
    
    # save paramiter estimation
    if isCombo and (not isZr75):
        fnamep = os.path.join(data_dir, runName+'-MCF7T47D-params.p')
    elif isCombo:
        fnamep = os.path.join(data_dir, runName+'-MCF7ZR75-params.p')
    elif isTransfer and (not isZr75):
        fnamep = os.path.join(data_dir, runName+'-T47D-params.p')
    elif isTransfer and isZr75:
        fnamep = os.path.join(data_dir, runName+'-ZR75-params.p')
    else:
        fnamep = os.path.join(data_dir, runName+'-MCF7-params.p')
    file = open(fnamep,'wb')
    pickle.dump(params, file)
    file.close()
    
    # clear the run directory of junk (don't keep anything permenent in
    # there)
    myModel.clearRunDirectory()
    
    # run time courses for each paramiterisation found for the stimulated
    # conditions
    tempParam = params[next(iter(params))]
    simTime = 120*60
    if prequilFlag:
        timeOffset = timeOffset[0]
        tempParam["zeroPeriod"] = timeOffset
        simTime = simTime + timeOffset
    if "noMet" in CLSet:
        tempParam = tempParam[[i for i in tempParam.columns 
                               if i.startswith("k") or i.endswith("_T")]]
    timeCourse = myModel.runTimeCourse(simTime,
                                       adjustParams=tempParam,
                                       stepSize=10,
                                       genReactions=True)
    
    # save time course for stimulated condition
    if isCombo and (not isZr75):
        fnametc = os.path.join(data_dir, runName+'-MCF7T47D-timeCourse.p')
    elif isCombo:
        fnametc = os.path.join(data_dir, runName+'-MCF7ZR75-timeCourse.p')
    elif isTransfer and (not isZr75):
        fnametc = os.path.join(data_dir, runName+'-T47D-timeCourse.p')
    elif isTransfer and isZr75:
        fnametc = os.path.join(data_dir, runName+'-ZR75-timeCourse.p')
    else:
        fnametc = os.path.join(data_dir, runName+'-MCF7-timeCourse.p')
    file = open(fnametc,'wb')
    if prequilFlag:
        timeCourseTemp = [tdf[tdf["Time"]>=timeOffset] for tdf
                          in timeCourse]
        for i in range(len(timeCourseTemp)):
            timeCourseTemp[i]["Time"] = (timeCourseTemp[i]["Time"] -
                          timeOffset)
        pickle.dump(timeCourseTemp, file)
    else:
        pickle.dump(timeCourse, file)
    file.close()
    
    # clear the run directory of junk (don't keep anything permenent in
    # there)
    myModel.clearRunDirectory()
    
    # run time courses for each paramiterisation found for the unstimulated
    # (steady state) conditions
    tempParam = params[next(iter(params))]
    tempParam["AA"]=0
    tempParam["Insulin"]=0
    tempParam.to_csv(os.path.join(data_dir,runName+'-MCF7-test4.csv'))
    if "noMet" in CLSet:
        tempParam = tempParam[[i for i in tempParam.columns 
                               if i.startswith("k") or i.endswith("_T")
                               or i in ["AA","Insulin"]]]
    if not prequilFlag:
        timeCourseSS = myModel.runTimeCourse(120*60,
                                             adjustParams=tempParam,
                                             stepSize=10,
                                             genReactions=True)
    
    # save time course for unstimulated condition
    if isCombo and (not isZr75):
        fnamess = os.path.join(data_dir, runName+'-MCF7T47D-timeCourseSS.p')
    elif isCombo:
        fnamess = os.path.join(data_dir, runName+'-MCF7ZR75-timeCourseSS.p')
    elif isTransfer and (not isZr75):
        fnamess = os.path.join(data_dir, runName+'-T47D-timeCourseSS.p')
    elif isTransfer and isZr75:
        fnamess = os.path.join(data_dir, runName+'-ZR75-timeCourseSS.p')
    else:
        fnamess = os.path.join(data_dir, runName+'-MCF7-timeCourseSS.p')
    file = open(fnamess,'wb')
    if prequilFlag:
        timeCourseSS = [tdf[tdf["Time"]<=timeOffset] for tdf in timeCourse]
        timeCourse = timeCourseTemp
    pickle.dump(timeCourseSS, file)
    file.close()
    
    myModel.clearRunDirectory()
    
    # target generating code 
    tempParam={key:val[val["RSS"]!=np.inf] for key, val in params.items()}
    candidates=RSSClusterEstimation(GFID(tempParam))[0]["size"]
    
    if isinstance(df,list):
        df = df[0]
    
    # generate dictionary of 'average' time courses for each phospho
    # protein 
    def genTCRefDict(data_relations,timeCourse,candidates,df,prefix=""):
        myDict={}
        for _, colList in data_relations.items():
            for col in colList:
                myDict[col] = pd.concat([timeCourse[i][
                        timeCourse[0].Time.isin(df.index)][prefix+col]
                        for i in range(candidates)],axis=1).mean(axis=1)
                myDict[col].index = list(df.index)
        for col in ["MTORC1_wb","MTORC2_wb"]:
            if col in timeCourse[0]:
                myDict[col] = pd.concat([timeCourse[i][
                        timeCourse[0].Time.isin(df.index)][prefix+col]
                        for i in range(candidates)],axis=1).mean(axis=1)
                myDict[col].index = list(df.index)
        return myDict
    
    if isCombo:
        myDict = [genTCRefDict(data_relations,timeCourse,candidates,df,
                               prefix) for prefix in ["ModA_","ModB_"]]
    else:
        myDict = genTCRefDict(data_relations,timeCourse,candidates,df)        
    
    # generate dictionary of 'average' last point for each phospho
    # protein in unstimulated 'steady state' time courses
    def genSSRefDict(data_relations,timeCourseSS,candidates,df,prefix=""):
        myDictSS={}
        for _, colList in data_relations.items():
            for col in colList:
                myDictSS[col] = pd.concat([timeCourseSS[i].iloc[[-1]][
                        prefix+col]
                    for i in range(candidates)],axis=1).mean(axis=1)
                myDictSS[col].index = [list(df.index)[0]]
        for col in ["MTORC1_wb","MTORC2_wb"]:
            if col in timeCourseSS[0]:
                myDictSS[col] = pd.concat([timeCourseSS[i].iloc[[-1]][
                        prefix+col]
                    for i in range(candidates)],axis=1).mean(axis=1)
                myDictSS[col].index = [list(df.index)[0]]
        return myDictSS
    
    if isCombo:
        myDictSS = [genSSRefDict(data_relations,timeCourse,candidates,df,
                                 prefix) for prefix in ["ModA_","ModB_"]]
    else:
        myDictSS = genSSRefDict(data_relations,timeCourseSS,candidates,df)
    
    # reopen RS file it might have changed during run time
    file = open(os.path.join(data_dir,'runSwitches.p'),'rb')
    RS = pickle.load(file)
    file.close()
    if isCombo:
        prefix = "Combo-"
    else:
        prefix = ""
    if not isTransfer:
        # save Targets for next rounds scalings        
        RS[prefix+"recycle-"+runName] = myDict
        RS[prefix+"recycleSS-"+runName] = myDictSS
        
        print("dynamic\n",myDict)
        print("steady\n",myDictSS)
        
        # save this rounds scalings for later use
        RS[prefix+"scale-"+runName] = myScale
    
        if prefix+"cycle-"+runName in RS.keys():
            RS[prefix+"cycle-"+runName] = RS[prefix+"cycle-"+runName] + 1
        else:
            RS[prefix+"cycle-"+runName] = 1
        
    RS["prequilFlag-"+runName] = prequilFlag
    
    # save RS file
    file = open(os.path.join(data_dir,'runSwitches.p'),'wb')
    pickle.dump(RS,file)
    file.close()