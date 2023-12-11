#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 12:33:29 2022

@author: peter
"""

from python.utilityTools import *
import pandas as pd

def getBlock(fPathList,cornerName,relative=True):
    df = pd.read_excel (resolvePath(fPathList, relative=True),
                        engine='openpyxl', header=None)
    potentials = df[df==cornerName].stack().index.tolist()
    retList = []
    for row, col in potentials:
        colEnd = df.loc[row, col:].isnull()
        colEnd = colEnd[colEnd]
        if len(colEnd)>0:
            colEnd = colEnd.index[0]-1
        else:
            colEnd = len(df.columns)-1            
        rowEnd = df.loc[row:, col].isnull()
        rowEnd = rowEnd[rowEnd]
        if len(rowEnd)>0:
            rowEnd = rowEnd.index[0]-1
        else:
            rowEnd = len(df.index)-1
        if row>0:
            if not all(list(df.loc[row-1, col:colEnd].isnull())):
                continue
        if rowEnd<len(df.index)-1:
            if not all(list(df.loc[rowEnd+1, col:colEnd].isnull())):
                continue
            
        if col>0:
            if not all(list(df.loc[row:rowEnd, col-1].isnull())):
                continue
        if colEnd<len(df.columns)-1:
            if not all(list(df.loc[row:rowEnd, colEnd+1].isnull())):
                continue
        df = df.loc[row:rowEnd,col:colEnd]
        df.columns = [i.strip() for i in df.iloc[0]]
        df = df.drop(df.index[0])
        df = df.set_index(cornerName)
        retList.append(df)
    if len(retList)==1:
        retList = retList[0]
    return retList

def listSepToTupals(iterable, sep = " ", discardInvarients = True):
    myList = [i.split(sep) for i in iterable]
    t = min([len(i) for i in myList])
    myList = [i[:2] for i in myList] # we now force 2 arguments
    """
    myList = [i[:t] for i in myList]
    if discardInvarients and len(myList)>0:
        t = [len(set([j[i] for j in myList]))>1 for i in range(t)]
        myList = [[j for j,k in zip(i,t) if k] for i in myList]
    """
    myList = [tuple(i) for i in myList]
    print(myList)
    return myList

def importData(fPathList, blockName = "average per experiment",
               timeMultiplyer = 60):
    df = getBlock(fPathList,blockName)
    df.index = pd.MultiIndex.from_tuples(listSepToTupals(df.index))
    df = {i:df.loc[i] for i in set(df.index.get_level_values(0))}
    for i in df:
        df[i].index = [float(j)*timeMultiplyer for j in df[i].index]
        df[i].index.name = "Time"
    return df

def scale(dictOfDf,scales):
    t = {k:v.copy() for k,v in dictOfDf.items()}
    for dk in t.keys():
        for k,v in scales.items():
            if k in t[dk].columns:
                t[dk][k] = t[dk][k]*v
    return t

def tableToDict(df,indexKey="index",indexVal=0,keyVar="dataKey"):
    if not indexKey in df.columns:
        return None
    df2 = df[df[indexKey]==indexVal].drop(columns=[indexKey])
    df2 = {i:df2[df2[keyVar]==i].drop(columns=[keyVar]) 
           for i in set(df2[keyVar])}
    df2 = {k:v.reset_index(drop=True) for k,v in df2.items()}
    return df2

def normalize(dictOfDf, acrossLines = False, noScale = []):
    t = {k:v.copy() for k,v in dictOfDf.items()}
    u = {k:{} for k,_ in dictOfDf.items()}
    for dk in t.keys():
        for k in [i for i in t[dk].columns if not i in noScale]:
            u[dk][k] = t[dk][k].mean()
            if not acrossLines:
                t[dk][k] = t[dk][k]/u[dk][k]
    if acrossLines:
        v = set.union(*[set(i.columns) for _,i in t.items()])
        v = {i:[k for _,j in t.items() for _,k in j[i].items()
                if i in j.columns] for i in v}
        v = {i:sum(j)/len(j) for i,j in v.items()}
        for dk in t.keys():
            for k in [i for i in t[dk].columns if not i in noScale]:
                t[dk][k] = t[dk][k]/v[k]
    return t, u