#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 15:55:55 2022

@author: peter
"""

from python.antimonyTools import *
from python.utilityTools import *
import json

cmdDefaults = {"antStr":"antString2.txt",
               #"antStr":"antStringRedF.txt",
               #"antStr":"modelV-arsenite.txt",
               #"antStr":"modelVI-arsenite.txt",
               #"antStr":"antStringToy.txt",
               "stabAdj":"stabAdj2.json",
               #"stabAdj":"stabAdjRed.json",
               #"stabAdj":"stabAdjSash.json",
               #"stabAdj":"stabAdjSashBackground.json",
               #"stabAdj":"stabToy.json",
               "outStr":"antString2M.txt"
               #"outStr":"antStringRedM.txt"
               #"outStr":"modelV-arseniteM.txt"
               #"outStr":"modelVI-arseniteM.txt"
               #"outStr":"antStringToyM.txt"
               }

cmdDict, cmdFlags = getCmdLineArgs()
cmdDict = parseCmdDict(cmdDict,cmdDefaults)
    
f = open(cmdDict["antStr"],'r')
antStr = f.read()
f.close()

myJ = """
{
  "simp":{"AA":0,"Insulin":0,"Cell":1},
  "ODEcull":["S6K","IRS1","PI3K", "Akt","TSC2pT1462","FourEBP1","PRAS40"],
  "totalSubs":{"S6K":"S6K_T-S6KpT389",
               "IRS1":"IRS1_T-IRS1a-IRS1pS636_639",
               "PI3K":"PI3K_T-pPI3K",
               "Akt":"Akt_T-AktpS473-AktpT308-AktpT308S473",
               "TSC2pT1462":"TSC2_T-TSC2",
               "FourEBP1":"FourEBP1_T-FourEBP1pT37_46",
               "PRAS40":"PRAS40_T-PRAS40pT246"},
  "solSteps":[{"eqs":["S6KpT389"],
               "sols":["S6KpT389"],
               "retain":["pmTORC1"]},
              {"eqs":["IRS1a","IRS1pS636_639"],
               "sols":["IRS1a","IRS1pS636_639"],
               "retain":["pmTORC1"]},
              {"eqs":["pPI3K"],
               "sols":["pPI3K"],
               "retain":["pmTORC1"]},
              {"eqs":["AktpS473","AktpT308","AktpT308S473"],
               "sols":["AktpS473","AktpT308","AktpT308S473"],
               "retain":["pmTORC1","pmTORC2"]},
              {"eqs":["PRAS40pT246"],
               "sols":["PRAS40pT246"],
               "retain":["pmTORC1","pmTORC2"]},
              {"eqs":["TSC2"],
               "sols":["TSC2"],
               "retain":["pmTORC1","pmTORC2"]},
              {"eqs":["FourEBP1pT37_46"],
               "sols":["FourEBP1pT37_46"],
               "retain":["pmTORC1"]}],
  "guessConstraints":{"mTOR_T":["mTOR", "mTORC1cyt", "mTORC1lys", "pmTORC1", 
                                "mTORC2", "pmTORC2"],
                      "RICTOR_T":["RICTOR", "mTORC2", "pmTORC2"],
                      "RPTOR_T":["RPTOR", "mTORC1cyt", "mTORC1lys", 
                                 "pmTORC1"]},
  "varToTrack":["Akt_wb", "AktpT308_wb", "AktpS473_wb", "PRAS40_wb", 
                "PRAS40pT246_wb", "S6K_wb", "S6KpT389_wb", "TSC2_wb", 
                "TSC2pT1462_wb", "IRS1_wb", "IRS1pS636_639_wb", "FourEBP1_wb", 
                "FourEBP1pT37_46_wb"],
  "stimulation":{"AA":1, "Insulin":1},
  "preStim":{"AA":0, "Insulin":0},
  "times":[0, 900, 1800, 3600, 5400, 7200]
}
"""

f = open(cmdDict["stabAdj"],'r')
myJ = f.read()
f.close()

myJ = json.loads(myJ)

test = antToSympyReader(antStr)
myODEs = {k:v.subs([(k,v) for k,v in myJ["simp"].items()]) 
          for k,v in test.ODEs.items()}
myODEs = {k:v for k,v in myODEs.items() if not k in myJ["ODEcull"]}

"""
totals = {"S6K":"S6K_T-S6KpT389",
          "IRS1":"IRS1_T-IRS1a-IRS1pS636_639",
          "PI3K":"PI3K_T-pPI3K",
          "Akt":"Akt_T-AktpS473-AktpT308-AktpT308S473",
          "TSC2pT1462":"TSC2_T-TSC2",
          "FourEBP1":"FourEBP1_T-FourEBP1pT37_46",
          "PRAS40":"PRAS40_T-PRAS40pT246"}
"""

totals = myJ["totalSubs"]
test2 = stabSolver(myODEs, additionalCons = totals)
"""
test2.subSolveSim(["S6KpT389"],["S6KpT389"],["pmTORC1"])
test2.subSolveSim(["IRS1a","IRS1pS636_639"],["IRS1a","IRS1pS636_639"],
                        ["pmTORC1"])
test2.subSolveSim(["pPI3K"],["pPI3K"],["pmTORC1"])
test2.subSolveSim(["AktpS473","AktpT308","AktpT308S473"],
                        ["AktpS473","AktpT308","AktpT308S473"],
                        ["pmTORC1","pmTORC2"])
test2.subSolveSim(["PRAS40pT246"],["PRAS40pT246"],["pmTORC1","pmTORC2"])
test2.subSolveSim(["TSC2"],["TSC2"],["pmTORC1","pmTORC2"])
test2.subSolveSim(["FourEBP1pT37_46"],["FourEBP1pT37_46"],["pmTORC1"])
"""
for i in myJ["solSteps"]:
    test2.subSolveSim(i["eqs"],i["sols"],i["retain"],debug=True)

"""
initCondDicts = {"mTOR_T":["mTOR", "mTORC1cyt", "mTORC1lys", "pmTORC1", 
                           "mTORC2", "pmTORC2"],
                 "RICTOR_T":["RICTOR", "mTORC2", "pmTORC2"],
                 "RPTOR_T":["RPTOR", "mTORC1cyt", "mTORC1lys", "pmTORC1"]}
"""

initCondDicts = myJ["guessConstraints"]

test2.asignICGuesses(initCondDicts)
newAntStr = test2.amendAntStr(antStr)

print(newAntStr)

"""
varToTrack = ["Akt_wb", "AktpT308_wb", "AktpS473_wb", "PRAS40_wb", 
              "PRAS40pT246_wb", "S6K_wb", "S6KpT389_wb", "TSC2_wb", 
              "TSC2pT1462_wb", "IRS1_wb", "IRS1pS636_639_wb", "FourEBP1_wb", 
              "FourEBP1pT37_46_wb"]
"""
varToTrack = myJ["varToTrack"]
#stimulation = {"AA":1, "Insulin":1}
stimulation = myJ["stimulation"]
#preStim = {"AA":0, "Insulin":0}
preStim = myJ["preStim"]
#times = [i*60 for i in [0,15,30,60,90,120]]
times = myJ["times"]
if "stabMod" in myJ:
    stabMod = myJ["stabMod"]
else:
    stabMod = {}

antStr2 = genRunInAntStr(myODEs,varToTrack,times,newAntStr,stimulation,preStim,
                         stabilityModifyer=stabMod)

f = open(cmdDict["outStr"],'w')
f.write(antStr2)
f.close()