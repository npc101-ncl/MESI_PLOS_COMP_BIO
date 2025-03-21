#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 16:01:01 2025

@author: peter
"""

from python.utilityTools import *
import pandas as pd
import tellurium as te

myDir = 'redPaper'

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

def extractBlock (keyword,myPath):
    data = pd.read_excel(myPath,engine='openpyxl',header=None) 
    t = [[(i,j) for j,_ in enumerate(data.columns)] 
         for i,_ in enumerate(data.index)]
    t = pd.DataFrame(t)
    t = t.to_numpy().flatten()
    t = t[(data==keyword).to_numpy().flatten()]
    if len(t)==0:
        return None
    data = data.iloc[t[0][0]:,t[0][1]:]
    data = data.set_index(0)
    t = data.index.isnull()
    t = [i for i in range(1,len(t)) if (t[i] and (not t[i-1]) and i>1)]
    if len(t)>0:
        data = data.iloc[:t[0]]
    t = list(data.isnull().all())
    t = [i for i in range(1,len(t)) if (t[i] and (not t[i-1]))]
    if (len(t)>0):
        data = data.iloc[:,:t[0]]
    t = data.index.isnull()
    t = [i for i in range(1,len(t)) if ((not t[i]) and t[i-1])]
    if (len(t)>0):
        u = data.iloc[:t[0]].transpose()
        for i in range(1,len(u)):
            for j in range(len(u.columns)):
                if pd.isna(u.iloc[i,j]):
                    u.iloc[i,j]=u.iloc[i-1,j]      
        u = pd.MultiIndex.from_frame(u)
        data = data.iloc[t[0]:]
        data.columns = u
    else:
        data.columns = data.iloc[0]
        data = data.iloc[1:]
    return data

def parsIndex (df):
    t = [i.split() for i in df.index]
    t = [(i[0].replace("-",""),int(i[1])) for i in t]
    t = pd.MultiIndex.from_tuples(t)
    df.index = t
    return df

CR = loadPick(['data',myDir,'comandRecord.p'],relative=True)

cmdDict = CR['cmdDict']

cmdDict = parseCmdDict(cmdDict,cmdDefaults)

if ('parallelAntStr' in cmdDict and cmdDict['parallelAntStr']!=0):
    antStr = cmdDict['parallelAntStr']
elif ('antStr' in cmdDict and cmdDict['antStr']!=0):
    antStr = cmdDict['antStr']
else:
    print('bad file')
    exit()
    
antStr = loadTxt(antStr, relative=True)

t = loadTxt(cmdDict['preAntStr'], relative=True)
sbmlStr = te.antimonyToSBML(t)
   
params = pd.read_csv(os.path.join(work_dir_U,'data',myDir,cmdDict['paramOut']))
params = [i for i in params.columns if not i in ['Unnamed: 0', 'RSS']]

PEInput = pd.read_csv(os.path.join(work_dir_U,'data',myDir,
                                   cmdDict["id"]+"PEInput.csv"))
obs = [i for i in PEInput.columns if not i in ['Unnamed: 0','Time']]
myTime = PEInput['Time'].iloc[0]

data = pd.read_excel(os.path.join(work_dir_U,cmdDict['data']),
                     engine='openpyxl')     

'normalized to ERK'
data = extractBlock('average per experiment',
                    os.path.join(work_dir_U,cmdDict['data']))
data = parsIndex(data)
t = pd.read_csv(os.path.join(work_dir_U,cmdDict['dataKey']),header=None)
t = t.set_index(0).squeeze().to_dict()
data = data.drop(columns=[i for i in data.columns if not (i in t.keys())])
data = data.rename(columns=t)
t = list(set(data.index.get_level_values(1)))
t.sort()
t = {j:i for i,j in enumerate(t)}
data = data/data.mean()
data = {k+"_DVT"+str(t[i[1]])+"_"+i[0]:v for i,j in data.iterrows() 
        for k,v in j.items()}
#print(data)

data = extractBlock('normalized to ERK',
                    os.path.join(work_dir_U,cmdDict['data']))
data = parsIndex(data)
dataKey = pd.read_csv(os.path.join(work_dir_U,cmdDict['dataKey']),header=None)
dataKey = dataKey.set_index(0).squeeze().to_dict()

elem = loadPick(["data",myDir,cmdDict['elemOut']], relative=True)
estimate = [i for i in elem['pre elements']['kineticParams'] 
            if not i in cmdDict["doNotEstimate"].split(":")]

data = data.drop(columns=[i for i in data.columns if (not i[0] in dataKey) 
            or dataKey[i[0]][:-3] not in elem['pre elements']['metabolites']])

exp_cond = {"stab":{"AA":0, "Insulin":0}, "dyn":{"AA":1, "Insulin":1}}
exp_cond = pd.DataFrame(exp_cond).transpose()
t = []
for i in pd.unique(data.index.get_level_values(0)):
    t.append(exp_cond.copy())
    t[-1]["conditionName"] = [j+" for "+i for j in exp_cond.index]
    t[-1].index = [j+"_"+i for j in t[-1].index]
    for j in elem['pre elements']['metabolites']:
        t[-1][j] = j+"_"+i
        estimate.append(j+"_"+i)
exp_cond = pd.concat(t)
exp_cond.index = exp_cond.index.rename("conditionId")

t = {k:v for k,v in dataKey.items() if k in data.columns.get_level_values(0)}
estimate.extend([i[:-2]+"scale" for _,i in t.items()])
t = [{"observableId":"obs_"+k.replace("-","_").replace("/","_").replace(" ","_"),
      "observableFormula":v[:-2]+"scale*"+v,
      "noiseFormula":"noiseParameter1_"+v} for k,v in t.items()]
obs = pd.DataFrame(t)

meas = []
times = pd.unique(data.index.get_level_values(1))
for i,j in data.iterrows():
    for k,v in j.items():
        if pd.isna(v):
            continue
        t = "obs_"+k[0].replace("-","_").replace("/","_").replace(" ","_")
        if t in list(obs["observableId"]):
            s = ("sd_"+k[0].replace("-","_").replace("/","_").replace(" ","_"))
            row = {"simulationConditionId":"dyn_"+i[0],
                   "time":i[1]*cmdDict['timeMul'],
                   "observableId":t,
                   "measurement":v,
                   "noiseParameters":s}
            if not s in estimate:
                estimate.append(s)
            meas.append(row)
            if (i[1]*cmdDict['timeMul']==0):
                for p in times:
                    row = {"simulationConditionId":"stab_"+i[0],
                           "time":p*cmdDict['timeMul'],
                           "observableId":t,
                           "measurement":v,
                           "noiseParameters":s}
                    meas.append(row)
meas = pd.DataFrame(meas)

estimate = pd.DataFrame({'parameterId': estimate})
estimate['parameterScale'] = "log10"
estimate['lowerBound'] = cmdDict['parameterLB']
estimate['upperBound'] = cmdDict['parameterUB']
estimate['nominalValue'] = None
estimate['estimate'] = 1

myYaml = """
format_version: 1
parameter_file: parameters.tsv
problems:
  - condition_files:
    - experimental_conditions.tsv
    measurement_files:
    - measurement_data.tsv
    observable_files:
    - observables.tsv
    sbml_files:
    - sbml_model.xml
"""   

forceDir(["PETab",myDir],relative=True) 

estimate.to_csv(resolvePath(["PETab",myDir,"parameters.tsv"],relative=True), 
                sep="\t")
meas.to_csv(resolvePath(["PETab",myDir,"measurement_data.tsv"],relative=True), 
                sep="\t")
obs.to_csv(resolvePath(["PETab",myDir,"observables.tsv"],relative=True), 
                sep="\t")
exp_cond.to_csv(resolvePath(["PETab",myDir,"experimental_conditions.tsv"],
                            relative=True), sep="\t")
f = open(resolvePath(["PETab",myDir,"sbml_model.xml"], relative=True), "w")
f.write(sbmlStr)
f.close()
f = open(resolvePath(["PETab",myDir,"petab.yaml"], relative=True), "w")
f.write(myYaml)
f.close()

