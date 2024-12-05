#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 21:43:21 2024

@author: peter
"""

import pandas as pd
import scipy.stats as stats
import math

myObs = pd.read_csv("zr75_obs.csv",index_col=["cellLine","time"])

t = []
myKey = {'4E-BP1':'FourEBP1_wb',
         '4E-BP1pT37/46':'FourEBP1pT37_46_wb',
         'Akt':'Akt_wb',
         'AktpS473':'AktpS473_wb',
         'AktpT308':'AktpT308_wb',
         'IRS1':'IRS1_wb',
         'IRS1pS636/639':'IRS1pS636_639_wb',
         'PRAS40':'PRAS40_wb',
         'PRAS40pS183':'PRAS40pS183_wb',
         'PRAS40pT246':'PRAS40pT246_wb',
         'S6K':'S6K_wb',
         'S6KpT229':'S6KpT229_wb',
         'S6KpT389':'S6KpT389_wb',
         'TSC2':'TSC2_wb',
         'TSC2pT1462':'TSC2pT1462_wb'}

for i,r in myObs.iterrows():
    for k,v in r.items():
        j = k.split(".")[0]
        t.append({"cellLine":i[0],"time":i[1],"species":myKey[j],"value":v})
myObs = t

t = pd.DataFrame(myObs)
meanObs = t.groupby(["cellLine","time","species"]).mean()
for i in t.index:
    j = t.loc[i]
    j = meanObs.loc[tuple(j[["cellLine","time","species"]])]["value"]
    t.loc[i,"value"]=(t.loc[i,"value"]-j)**2
t = t[["species","value"]].groupby(["species"]).mean()
stdObs = t**(0.5)

t = [('figures/GPMCF7/exporttest7.csv',35+13),
     ('figures/paper/exporttest6.csv',39+13),
     ('figures/redPaper/exporttest4.csv',17+12),
     ('figures/GPMCF7/exporttestZR751.csv',4)]

def getChiSq (fPath, myObs, stdObs):
    myExp = pd.read_csv(fPath,
                        index_col=["cellLine","time","species"])
    chiSq = 0
    RSS = 0
    comCount = 0
    for i in myObs:
        if (i["cellLine"],i["time"],i["species"]) in myExp.index:
            j = myExp.loc[(i["cellLine"],i["time"],i["species"])]['value']
            k = stdObs.loc[i["species"]]['value']
            chiSq += ((i['value']-j)/k)**2
            RSS += (i['value']-j)**2
            comCount += 1
    return chiSq, comCount, RSS

def getChiSq2 (fPath, myObs, stdObs):
    myExp = pd.read_csv(fPath,
                        index_col=["cellLine","time","species"])
    RSS = 0
    comCount = 0
    for i in myObs:
        if (i["cellLine"],i["time"],i["species"]) in myExp.index:
            j = myExp.loc[(i["cellLine"],i["time"],i["species"])]['value']
            k = stdObs.loc[i["species"]]['value']
            #chiSq += ((i['value']-j)/k)**2
            RSS += (i['value']-j)**2
            comCount += 1
    df = pd.merge(myExp,pd.DataFrame(myObs),
                  on=["cellLine","time","species"],
                  suffixes=("_exp","_obs"))
    df = df.groupby(["cellLine","time","species"]).agg(
            {'value_exp': 'first', 'value_obs': ['mean',"std"]})
    df[('chi_sq','term')] = ((df[('value_obs','mean')]-
                              df[('value_exp','first')])/
                             df[('value_obs','std')])**2
    df = df[[('chi_sq','term')]]
    df.columns = ['chi_sq']
    df = df.groupby(["cellLine","species"]).agg(['sum','count'])
    df = df.reset_index()
    df = df.pivot(index='species', columns='cellLine', 
                  values=[('chi_sq','sum'), ('chi_sq','count')])
    return df, comCount, RSS
    
        

v = []
for myPath, params in t:
    #chiSq, comCount, RSS = getChiSq (myPath, myObs, stdObs)
    chiSq, comCount, RSS = getChiSq2 (myPath, myObs, stdObs)
    v.append(comCount-params)
    AIC = (2*params+comCount*math.log(RSS/comCount)+
           comCount*math.log(2*math.pi)+comCount)
    #print(myPath,comCount-params,chiSq,stats.chi2.sf(chiSq,df=comCount-params),
    #      AIC)
    print(myPath, comCount, "{:.2e}".format(RSS), "{:.2e}".format(AIC), "\n\n")
    print(chiSq,"\n\n")
    
import matplotlib.pyplot as plt
for j in v:
    x = [(i/10,stats.chi2.sf(i/10,df=j)) for i in range(3000,4000)]
    x,y = zip(*x)
    plt.plot(x, y)

