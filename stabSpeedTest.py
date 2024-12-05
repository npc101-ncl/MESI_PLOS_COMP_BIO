#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 23:28:03 2024

@author: peter
"""

from python.pycotoolsHelpers import *
from python.utilityTools import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn

params = ["kIRS1Act","kIRS1Inact","kIRS1Phos","kPI3KPhos","kPI3KDephos",
          "kAktPhos_kcat","kAktDephos","kAktByTor2Phos","kAktByTor2Dephos",
          "kTSC2Phos","kTSC2PhosBoost","kTSC2Dephos","kmTORC1cytToLys",
          "kmTORC1LysToCyt","kmTORC1Phos","kmTORC1Dephos","kPras40Phos",
          "kPras40Dephos","kFourEBP1Phos","kFourEBP1Dephos","kS6KPhos",
          "kS6KDephos","kRAS","kmTORC2Dephos","kmTORC2DephosByS6K",
          "kmTOR_RIC_Comb","kmTORC2_Dis","kmTOR_RPT_Comb","kmTORC1_Dis",
          "InsulinB","AAB","RICTOR_T","FourEBP1_T","S6K_T","PRAS40_T",
          "RPTOR_T","IRS1_T","mTOR_T","PI3K_T","TSC2_T","Akt_T"]
maxTime = 10
timeRes = 0.01
generate = False

def getTimes(myTC):
    sympyTimes = []
    for j in myTC:
        k=j[j["PTSpeed"]==0]
        if len(k)==len(j):
            sympyTimes.append(None)
            continue
        sympyTimes.append(k.iloc[-1]['Time'])
    return sympyTimes

if generate:
    myCond = pd.DataFrame(np.random.uniform(0.01,100,(200,len(params))),
                          columns=params)
    
    addCopasiPath("/Applications/copasi")
              
    
    myModel = modelRunner(antString=loadTxt(["antStringMSympy.txt"],
                                            relative=True),
                          run_dir=resolvePath(["copasiRuns","sympyTest"],
                                              relative=True))
    myModel.clearRunDirectory()
    timeCourse = myModel.runTimeCourse(maxTime,
                                       adjustParams=myCond,
                                       stepSize=timeRes)
    
    sympyTimes = getTimes(timeCourse)
    myCond["With Sympy Adjusment"] = sympyTimes
    
    myModel = modelRunner(antString=loadTxt(["antStringMNoSympy.txt"],
                                            relative=True),
                          run_dir=resolvePath(["copasiRuns","sympyTest"],
                                              relative=True))
    myModel.clearRunDirectory()
    timeCourse = myModel.runTimeCourse(maxTime,
                                       adjustParams=myCond,
                                       stepSize=timeRes)
    
    noSympyTimes = getTimes(timeCourse)
    myCond["Without Sympy Adjusment"] = noSympyTimes

    myCond.to_csv("sympyTest.csv",index=False)
else:
    myCond = pd.read_csv("sympyTest.csv",index_col=False)
 
seaborn.set(style='whitegrid')
 
myPlot = seaborn.scatterplot(x="Without Sympy Adjusment",
                    y="With Sympy Adjusment",
                    data=myCond)
myPlot.set(xscale="log", yscale="log", aspect=1)
myPlot.set_xlabel('Without Sympy Adjusment (s)')
myPlot.set_ylabel('With Sympy Adjusment (s)')
fig = myPlot.get_figure()
fig.tight_layout()
fig.savefig("sympyTest.png")