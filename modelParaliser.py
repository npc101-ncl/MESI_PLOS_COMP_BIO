#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 13:37:11 2022

@author: peter
"""

from python.utilityTools import *
from python.pycotoolsHelpers import *
import json

cmdDefaults = {
        #"json":"torParaliser.json"
        #"json":"redParaliser.json"
        #"json":"sashaParaliser.json"
        #"json":"sashaParaliserIns.json"
        #"json":"sashaBrokenParaliser.json"
        #"json":"sashaParaliserBackground.json"
        #"json":"torMCF7Paraliser.json"
        "json":"torZR75Paraliser.json"
        #"json":"toyParaliser.json"
        #"json":"redParaliserZR75.json"
        }

cmdDict, cmdFlags = getCmdLineArgs()
cmdDict = parseCmdDict(cmdDict,cmdDefaults)

myJ = """
{
  "preAntStr":"antString2.txt",
  "postAntStr":"antString2M.txt",
  "outAntStr":"antString2MP3.txt",
  "CDID":"multiTest",
  "modelName":"AktModelWithMMRateLaws",
  "modelVariations":["MCF7","ZR751","T47D"],
  "bigNum":10000,
  "path":"/Applications/copasi"
}
"""

f = open(cmdDict["json"],'r')
myJ = f.read()
f.close()

myJ = json.loads(myJ)

if "path" in cmdDict:
    myJ["path"]=cmdDict["path"]

"""
preAntStr = "antString2.txt"
postAntStr = "antString2M.txt"
outAntStr = "antString2MP3.txt"
CDID = "multiTest"
modelName = "AktModelWithMMRateLaws"
modelVariations = ["MCF7","ZR751","T47D"]
bigNum = 10000
"""

preAntStr = myJ["preAntStr"]
postAntStr = myJ["postAntStr"]
outAntStr = myJ["outAntStr"]
CDID = myJ["CDID"]
modelName = myJ["modelName"]
modelVariations = myJ["modelVariations"]
if isinstance(modelVariations,dict):
    overrideParams = modelVariations
    modelVariations = list(modelVariations.keys())
else:
    overrideParams = {}
bigNum = myJ["bigNum"]

alternate = True

if "path" in myJ:
    addCopasiPath(myJ["path"])
if "setScale" in myJ:
    setScale = myJ["setScale"]
else:
    setScale = {}

preText = loadTxt([preAntStr], relative=True)
preModel = modelRunner(antString=preText, 
                       run_dir=resolvePath(["copasiRuns", CDID],relative=True))
preModel.clearRunDirectory()

preEliments = preModel.getModelEliments()
preParam = preModel.extractModelParam()

postText = loadTxt([postAntStr], relative=True)
postModel = modelRunner(antString=postText, 
                       run_dir=resolvePath(["copasiRuns", CDID],relative=True))
postModel.clearRunDirectory()

postEliments = postModel.getModelEliments()

linkedParams = preEliments["kineticParams"]
if "hidenParams" in myJ:
    hidenParams = myJ["hidenParams"]
else:
    hidenParams = []
linkedParams = list(set(linkedParams)- set(hidenParams))

outputs = list(set(postEliments["kineticParams"])-set(linkedParams))
if len(overrideParams)>0:
    linkedParams = list(set(linkedParams)-
                        set.union(*[set(i.keys()) for i 
                                    in overrideParams.values()]))
outputs = {i:[j for j in outputs if j.startswith(i)] 
           for i in preEliments["assignments"]}
outputs = {k:v for k,v in outputs.items() if len(v)>0}
outputsL = [j for _,i in outputs.items() for j in i]

pusudoTime = list(set(postEliments["metabolites"])-
                  set(preEliments["metabolites"]))[0]

events = [i.strip() for i in postText.splitlines() 
          if i.strip().startswith("at ")]
events = [i.split(":") for i in events if any([j in i for j in outputsL])]
events = [i[0].split("(") for i in events if len(i)>0]
events = [i[1].split(")") for i in events if len(i)==2]
events = [i[0] for i in events if len(i)==2]
events = [[j for j in i.split(pusudoTime) if len(j)>1] for i in events]
events = [i[0] for i in events if len(i)==1]
events = ["".join([j for j in i if j in "0123456789."]) for i in events]
events = [float(i) for i in events]
events = max(events)

inputs = list(set(postEliments["kineticParams"])-
              set(preEliments["kineticParams"])-set(outputsL))
inputs = [i for i in inputs 
          if any([i.startswith(j) for j in preEliments["metabolites"]])]

if not "joinInputs" in myJ:
    pass
elif isinstance(myJ["joinInputs"],bool):
    if myJ["joinInputs"]:
        linkedParams = linkedParams+inputs
        preParam.update({i:"1" for i in inputs if not i in preParam})
        inputs = []
elif isinstance(myJ["joinInputs"],list):
    for i in myJ["joinInputs"]:
        if i in inputs:
            linkedParams.append(i)
            inputs.remove(i)
            if not i in preParam:
                preParam[i]="1"


print(inputs)
print(outputs)
print(pusudoTime)
print(events)
print(linkedParams)

modelImps = [i+": "+modelName+"();" for i in modelVariations]
paramiterLinks = [j+"."+i+" is "+i+";" for i in linkedParams for j 
                  in modelVariations]
paramiterSets = [i+"="+preParam[i]+";" for i in linkedParams]
paramiterSets += [i+"=1.0" for _,i in setScale.items()]
overideSets = [k+"."+i+"="+str(j)+";" for k,v in overrideParams.items() 
               for i,j in v.items()]
inputDefs = [j+"."+i+" is "+i+"_"+j+";" for i in inputs for j 
             in modelVariations]
outputLinks = [j+"."+i+" is "+i+"__"+j+";" for k,v in outputs.items() for i
               in v for j in modelVariations]



if alternate:
    meanDefs = [k+"_mean"+"=("+"+".join([j+"__"+i for j in v 
                                         for i in modelVariations])+")/"+
                str(len(v)*len(modelVariations)) for k,v in outputs.items()
                if not k in setScale]
else:
    meanDefs = [k+"_"+i+"=("+"+".join([j+"__"+i for j in v])+")/"+str(len(v)) 
               for k,v in outputs.items() for i in modelVariations 
               if k not in setScale]
meanEvt = " && ".join(["("+pusudoTime+"_"+i+">"+str(events+1)+")" for i 
                       in modelVariations])
if len(meanDefs)==0:
    meanEvt = ""
else:
    meanEvt = "at ("+meanEvt+"): "+", ".join(meanDefs)
if alternate:
    meanDefs = [k+"_mean = "+str(bigNum)+" ;" for k,v in outputs.items()
                if k not in setScale]
else:
    meanDefs = [k+"_"+i+" = "+str(bigNum)+" ;" for k,v 
                in outputs.items() for i in modelVariations 
                if k not in setScale]

psudoTimeLinks = [j+"."+pusudoTime+" is "+pusudoTime+"_"+j+";" 
                  for j in modelVariations]
if alternate:
    outDefs = [i+"_"+j+"="+i+"__"+j+"/"+k+"_mean" for k,v in outputs.items()
               for i in v for j in modelVariations if k not in setScale]
else:
    outDefs = [i+"_"+j+"="+i+"__"+j+"/"+k+"_"+j for k,v in outputs.items()
               for i in v for j in modelVariations if k not in setScale]

if not alternate:
    rWeightLine = "rWeight = 1;"
    ratioDefs = [k+"_"+i+"_ratio = rWeight*"+k+"_"+i+"/("+
                 "+".join([k+"_"+j for j in modelVariations])+")" 
                 for k,_ in outputs.items() for i in modelVariations]
else:
    rWeightLine = ""
    ratioDefs = []
normEvt = " && ".join(["("+pusudoTime+"_"+i+">"+str(events+2)+")" for i 
                       in modelVariations])
if len(outDefs+ratioDefs)==0:
    normEvt = ""
else:
    normEvt = "at ("+normEvt+"): "+", ".join(outDefs+ratioDefs)
outDefs = [i+"_"+j+" = "+str(bigNum)+" ;" for k,v in outputs.items() for i 
           in v for j in modelVariations if k not in setScale]
outDefs += [i+"_"+j+":="+i+"__"+j+"*"+setScale[k] for k,v in outputs.items() 
            for i in v for j in modelVariations if k in setScale]
if not alternate:
    ratioDefs = [k+"_"+i+"_ratio = " + str(bigNum)+" ;" 
                 for k,_ in outputs.items() for i in modelVariations]

pModelLines = (modelImps+paramiterLinks+paramiterSets+overideSets+
               [rWeightLine]+outputLinks+psudoTimeLinks+inputDefs+meanDefs+
               outDefs+ratioDefs+[meanEvt,normEvt])

pModelLines = ["model paralized_model()"]+["\t"+i for i in pModelLines]+["end"]
pModelLines = "\n".join(pModelLines)

outText = postText+"\n\n"+pModelLines

f = open(resolvePath([outAntStr],relative=True), "w")
f.write(outText)
f.close()