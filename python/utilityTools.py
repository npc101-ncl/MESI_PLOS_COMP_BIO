#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 12:05:19 2020

@author: peter
"""

import sys, os, pickle, time, re

cmdLineArgU = sys.argv[1:]
work_dir_U = os.path.dirname(os.path.abspath(__file__))
work_dir_U = os.path.dirname(work_dir_U)

def getCmdLineArgs():
    outSet = set({})
    outDict = {}
    for arg in cmdLineArgU:
        temp = arg.split(":", 1)
        if len(temp)==2:
            outDict[temp[0]] = temp[1]
        else:
            outSet.add(arg)
    return outDict, outSet

def convertType(value,tagetType):
    if tagetType==int:
        return int(value)
    elif tagetType==float:
        return float(value)
    else:
        return str(value)
    
def listify(arg):
    if arg is None:
        out = None
    elif isinstance(arg,list):
        out = arg
    else:
        out = [arg]
    return out

def parseCmdDict(cmdDict,defaults):
    outList = {}
    for key, default in defaults.items():
        if not key in cmdDict:
            outList[key] = default
            continue
        myType = type(default)
        outList[key] = convertType(cmdDict[key],myType)
    for key, value in cmdDict.items():
        if key not in defaults:
            outList[key] = value
    return outList

def extractCompositArgument(cmdArg,sep=":",dictSep=None,forceType=str):
    myOut = cmdArg.split(sep)
    if sep == dictSep:
        myOut = {myOut[2*i]:myOut[2*i+1] for i in range(len(myOut)//2)}
    elif dictSep is not None:
        myOut = [i.split(dictSep,1) for i in myOut]
        myOut = {i[0]:i[1] for i in myOut}
    if isinstance(myOut,list):
        myOut = [convertType(i,forceType) for i in myOut]
    elif isinstance(myOut,dict):
        myOut = {k:convertType(v,forceType) for k,v in myOut.items()}
    return myOut

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

def resolvePath(fPath,requireExt=None,relative=False):
    if isinstance(fPath,str):
        pathList=[fPath]
    elif isinstance(fPath,list):
        pathList = fPath
    else:
        return None
    if relative:
        pathList.insert(0,work_dir_U)
    pathList = os.path.join(*pathList)
    if isinstance(requireExt,str):
        if os.path.splitext(pathList)[1] != "."+requireExt:
            return None
    return pathList

def savePick(fPath, data, relative=False, mkDir=True):
    pathList = resolvePath(fPath,requireExt="p",relative=relative)
    if pathList is None:
        return None
    dirName = os.path.dirname(pathList)
    if mkDir and not os.path.isdir(dirName):
        os.makedirs(dirName)
    file = open(pathList, 'wb')
    pickle.dump(data, file)
    file.close()
    
def loadPick(fPath, relative=False):
    pathList = resolvePath(fPath,requireExt="p",relative=relative)
    if pathList is None:
        return None
    if not os.path.isfile(pathList):
        return None
    file = open(pathList, 'rb')
    data = pickle.load(file)
    file.close()
    return data

def loadTxt(fPath, relative=False):
    pathList = resolvePath(fPath,requireExt="txt",relative=relative)
    if pathList is None:
        return None
    if not os.path.isfile(pathList):
        return None
    f = open(pathList, "r")
    text = f.read()
    f.close()
    return text

def getFilesIn(fPath, relative = False):
    if isinstance(fPath,str):
        pathList = [fPath]
    elif isinstance(fPath,list):
        pathList = fPath
    else:
        return None
    if relative:
        pathList.insert(0,work_dir_U)
    pathList = os.path.join(*pathList)
    return [f for f in os.listdir(pathList)
            if os.path.isfile(os.path.join(pathList, f))]
    
def getFileTargets(source_name, endings, names):
    if isinstance(endings, str):
        myEnds = [endings]
    elif isinstance(endings, list):
        myEnds = endings
    else:
        return None
    myList = []
    for myEnd in myEnds:
        pattern = re.compile("^"+source_name+"-"+"(.*)"+"-"+myEnd+".p$")
        myList2 = []
        for name in names:
            myMatch = pattern.fullmatch(name)
            if myMatch is not None:
                myList2.append(myMatch.group(1))
        myList.append(set(myList2))
    return list(set.intersection(*myList))

def getEndTime(seconds = 60*60*47):
    return time.time()+seconds

def repairPath(path):
    refDir = os.path.split(work_dir_U)[1]
    brokePath = os.path.split(path)
    fixedPath = []
    while brokePath[0]!="/":
        if (brokePath[1]==refDir):
            break
        else:
            fixedPath=[brokePath[1]]+fixedPath
            brokePath = os.path.split(brokePath[0])
    else:
        return None
    fixedPath=[work_dir_U]+fixedPath
    fixedPath = os.path.join(*fixedPath)
    return fixedPath

def forceDir(fPath,relative=False):
    myDir = resolvePath(fPath,relative=relative)
    if not os.path.isdir(myDir):
        os.makedirs(myDir)
        
        