#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 16:43:11 2020

@author: peter
"""

import pickle, os, re, math
from python.visualisationTools import *
from python.analysisTools import *
from python.utilityTools import *
from python.pycotoolsHelpers import extractAntReactions
from scaleFinderXPrimer import importData, makeInterpolatedTable, applyScale
import numpy as np

working_directory = os.path.dirname(os.path.abspath(__file__))
source_name = "red5"
maxPerGrid = 20

RS = loadPick(["data",source_name,"runSwitches.p"], relative=True)

"tor1-S6K-T47D-MCF7-params.p"

files = getFilesIn(["data",source_name], relative = True)

targets = getFileTargets(source_name,
                         ["params", "timeCourse", "timeCourseSS"], files)

dispTables = {}

maxCurves = 3

isCombo = len([1 for k in RS.keys() if k.startswith("Combo-")])>0

class myHTMLtable:
    def __init__(self,name=None):
        self.df = pd.DataFrame()
        self.name = name
    
    def add(self,col,row,entry):
        if not col in self.df.columns:
            self.df[col]=None
        if not row in self.df.index:
            self.df = self.df.append(pd.Series(name=row))
        self.df.at[row, col] = entry
        
    def print(self):
        htmlSTR = "<table style='width:100%'>\n"
        htmlSTR += "<th>"+self.name+"</th>\n"
        for col in self.df.columns:
            htmlSTR += "<th>"+col+"</th>\n"
        for index, row in self.df.iterrows():
            htmlSTR += "<tr>\n"
            htmlSTR += "<td>"+index+"</td>\n"
            for _, item in row.items():
                htmlSTR += "<td>"+item+"</td>\n"
            htmlSTR += "</tr>\n"
        htmlSTR += "</table>"
        return htmlSTR
    
    def fillNA(self,fillValue=""):
        self.df = self.df.fillna(value=fillValue)
        
class imgTag:
    def __init__(self,src):
        self.src = src
        
    def print(self):
        return "<img src='"+self.src+"'>"
                

def getRawClustorNames(names,clustor,reactionPrefix=None):
    clustorNames = [name for name in names if name.startswith(clustor) and
                    not (name.endswith("_wb") or name.endswith("_sum"))]
    if isinstance(reactionPrefix,str):
        reactionNames = [reactionPrefix]
    else:
        reactionNames = reactionPrefix
    if isinstance(reactionNames,list):
        reactionNames = [name for name in names if
                         any([name.startswith(prefix) for prefix
                              in reactionNames])]
        reactionNames = {name:(name.split("__")[1]) for name
                         in reactionNames}
        reactions = extractAntReactions(RS["antimony_string"])
        killKeys = []
        for k in reactionNames.keys():
            candName = [r for r in reactions if r["name"]==reactionNames[k]]
            if len(candName)!=1:
                killKeys.append(k)
                continue
            else:
                candName = candName[0]
            if any([i["var"].startswith(clustor) for i in candName["LHS"]]):
                pass
            elif any([i["var"].startswith(clustor) for i in candName["RHS"]]):
                pass
            else:
                killKeys.append(k)
                continue
            if candName['irreversible']:
                candName = ("+".join([i["var"] for i in candName["LHS"]])
                            + "=>" +
                            "+".join([i["var"] for i in candName["RHS"]]))
            else:
                candName = ("+".join([i["var"] for i in candName["LHS"]])
                            + "->" +
                            "+".join([i["var"] for i in candName["RHS"]]))
            reactionNames[k] = candName
        for k in killKeys:
            del(reactionNames[k])
    else:
        reactionNames = None
    return clustorNames, reactionNames
            

def pathToLocal(path,working_directory):
    stopDir = os.path.split(working_directory)[1]
    pathList = [working_directory]
    remPath = path
    while remPath != "":
        remPath, testDir = os.path.split(remPath)
        if testDir!=stopDir:
            pathList.append(testDir)
        else:
            return os.path.join(*pathList)
    return None

class myHTMLpage:
    def __init__(self,title):
        self.title = title
        self.body = []
        self.head = []
        
    def addToBody(self,item):
        self.body.append(item)
        
    def addToHead(self,item):
        self.head.append(item)
    
    def print(self):
        return """<!DOCTYPE html>
        <html>
        <head>
        <title>"""+self.title+"""</title>
        """+"\n".join(self.head)+"""
        </head>
        <body>\n"""+"\n".join(self.body)+"""\n</body>
        </html>"""

clustors = ["IRS", "Akt", "TSC2", "S6K"]

wTable = myHTMLtable("Waterfall")
cluTable = {c:myHTMLtable(c) for c in clustors}
wbTable = myHTMLtable("Western Blot")
torTable = myHTMLtable("mTOR")
depTable = myHTMLtable("mTOR depleation")
scaleTable = myHTMLtable("Scales")
violinTable = myHTMLtable("Parameter Distrobution")

figPath = os.path.join(working_directory, "figures", source_name)
if not os.path.isdir(figPath):
    os.makedirs(figPath)

print(targets)
if isCombo:
    subDatas = [{"MCF7T47D":["MCF7","T47D"],
                 "MCF7ZR75":["MCF7","ZR75"]}[target.split("-")[1]]
                for target in targets]
    targets = [[target for _ in subData] for target, subData 
               in zip(targets,subDatas)]
    targets = [j for i in targets for j in i]
    subDatas = [j for i in subDatas for j in i]
else:
    subDatas = [None for _ in target]

for target, subData in zip(targets,subDatas):
    if isCombo:
        pathInsert = target.split("-")[0]+"-"+subData
    else:
        pathInsert = target
    params = loadPick(["data",source_name,
                      source_name+"-"+target+"-params.p"], relative=True)
    PEVis = parameterEstimationVisualiser(params)
    figPath = (source_name+"-"+pathInsert+"-waterfall.png")
    myImg = imgTag(figPath)
    wTable.add(pathInsert,"",myImg.print())
    figPath = os.path.join(working_directory, "figures", source_name,
                           figPath)
    PEVis.waterFall(save = figPath)
    timeCourse = loadPick(["data",source_name,
                           source_name+"-"+target+"-timeCourse.p"],
                          relative=True)
    timeCourseSS = loadPick(["data",source_name,
                             source_name+"-"+target+"-timeCourseSS.p"],
                            relative=True)
    if subData is None:
        pass
    else:
        if subData == "MCF7":
            toDrop = ["ModB_","B_IrRevRe__"]
            prefixToSwap = {"ModA_":"",
                            "A_IrRevRe__":"IrRevRe__"}
        else:
            toDrop = ["ModA_","A_IrRevRe__"]
            prefixToSwap = {"ModB_":"",
                            "B_IrRevRe__":"IrRevRe__"}
        for TC2 in [timeCourse,timeCourseSS]:
            for TC in TC2:
                for i in toDrop:
                    TC.drop(columns=[col for col in TC.columns
                                     if col.startswith(i)],
                            inplace=True)
                for k, v in prefixToSwap.items():
                    TC.rename(columns={col:v+col[len(k):] for col
                                       in TC.columns
                                       if col.startswith(k)},
                              inplace=True)
        print(subData,timeCourse[0].columns)
    data = {key:value for key, value in RS["data_filenames"].items()
            if str.lower(key) in str.lower(target)}
    if len(data)==1:
        if isCombo:
            theScale = RS["Combo-scale-"+source_name+"-"+
                          target.split("-")[0]]
            dataName = {"MCF7":"MCF-7","T47D":"T47D",
                        "ZR75":"ZR-75-1"}[subData]
        else:
            theScale = RS["scale-"+source_name+"-"+target.split("-")[0]]
            dataName = {"MCF7":"MCF-7","T47D":"T47D",
                        "ZR75":"ZR-75-1"}[target.split("-")[1]]
        data = pathToLocal(GFID(data),working_directory)
        data = importData(dataName,data,RS["columns_from_data"])
        data = applyScale(data, theScale)
    else:
        data = None
    try:
        firstCluster = RSSClusterEstimation(GFID(params))[0]
        BestClusterSize = firstCluster["size"]
        firstRSSdivide = firstCluster["maxRSS"]
        maxCurves = min(maxCurves,BestClusterSize)
    except:
        BestClusterSize = maxCurves
    maxRSS = GFID(params).iloc[maxCurves-1]["RSS"]
    showList = list(range(maxCurves))
    if data is not None:
        timeCourseVis = timeCourseVisualiser(timeCourse)
        varSel = list(timeCourse[0].columns)
        varSel = [i for i in varSel if i in data.columns]
        for i in range(math.ceil(len(timeCourse[0].columns)/maxPerGrid)):
            varSel = varSel[i*maxPerGrid:(i+1)*maxPerGrid]
            if len(varSel) == 0:
                break
            figPath = (source_name+"-"+pathInsert+"-timeCourseWB"+str(i)+
                       ".png")
            myImg = imgTag(figPath)
            wbTable.add(pathInsert,"dynamic "+str(i),myImg.print())
            figPath = os.path.join(working_directory, "figures", source_name,
                                   figPath)
            timeCourseVis.multiPlot(indexSelect=showList, varSelect=varSel,
                                    compLines=data, save=figPath)
        dataClip = data.iloc[[0]]
        timeCourseVis = timeCourseVisualiser(timeCourseSS)
        varSel = list(timeCourseSS[0].columns)
        varSel = [i for i in varSel if i in data.columns]
        for i in range(math.ceil(len(timeCourseSS[0].columns)/maxPerGrid)):
            varSel = varSel[i*maxPerGrid:(i+1)*maxPerGrid]
            if len(varSel) == 0:
                break
            figPath = source_name+"-"+pathInsert+"-timeCourseWBSS"+str(i)+".png"
            myImg = imgTag(figPath)
            wbTable.add(pathInsert,"steady "+str(i),myImg.print())
            figPath = os.path.join(working_directory, "figures", source_name,
                                   figPath)
            timeCourseVis.multiPlot(indexSelect=showList, varSelect=varSel,
                                    compLines=dataClip, save=figPath)
            figPath = source_name+"-"+pathInsert+"-barT0WBSS"+str(i)+".png"
            figPath = os.path.join(working_directory, "figures", source_name,
                                   figPath)
            varSelB = ["Akt_wb","PRAS40_wb","S6K_wb","TSC2_wb","IRS1_wb",
                       "FourEBP1_wb"]
            timeCourseVis.barChart(0, indexSelect = showList,
                                   varSelect = varSel,
                                   compLines = dataClip,
                                   save=figPath)
            figPath = source_name+"-"+pathInsert+"-barTEWBSS"+str(i)+".png"
            figPath = os.path.join(working_directory, "figures", source_name,
                                   figPath)
            timeCourseVis.barChart(list(dataClip.index)[-1],
                                   indexSelect = showList,
                                   varSelect = varSel, compLines = dataClip,
                                   save=figPath)
    timeCourseVis = timeCourseVisualiser(timeCourse)
    varSel = list(timeCourse[0].columns)
    varSel = [i for i in varSel if "TOR" in i and (not i.startswith("k")) 
              and (not i.endswith("_F")) and (not i.endswith("_z"))]
    for i in range(math.ceil(len(timeCourse[0].columns)/maxPerGrid)):
        varSel = varSel[i*maxPerGrid:(i+1)*maxPerGrid]
        if len(varSel) == 0:
            break
        figPath = source_name+"-"+pathInsert+"-timeCourseMTOR"+str(i)+".png"
        myImg = imgTag(figPath)
        torTable.add(pathInsert,"dynamic "+str(i),myImg.print())
        figPath = os.path.join(working_directory, "figures", source_name,
                               figPath)
        timeCourseVis.multiPlot(indexSelect=showList, varSelect=varSel,
                                save=figPath)
    timeCourseVis = timeCourseVisualiser(timeCourseSS)
    varSel = list(timeCourseSS[0].columns)
    varSel = [i for i in varSel if "TOR" in i and not i.startswith("k")
              and (not i.endswith("_F")) and (not i.endswith("_z"))]
    for i in range(math.ceil(len(timeCourseSS[0].columns)/maxPerGrid)):
        varSel = varSel[i*maxPerGrid:(i+1)*maxPerGrid]
        if len(varSel) == 0:
            break
        figPath = source_name+"-"+pathInsert+"-timeCourseMTORSS"+str(i)+".png"
        myImg = imgTag(figPath)
        torTable.add(pathInsert,"steady "+str(i),myImg.print())
        figPath = os.path.join(working_directory, "figures", source_name,
                               figPath)
        timeCourseVis.multiPlot(indexSelect=showList, varSelect=varSel,
                                save=figPath)
    for clustor in clustors:
        varSel, rVarSel = getRawClustorNames(list(timeCourse[0].columns),
                                             clustor,
                                             reactionPrefix = "IrRevRe__")
        repeats={}
        varSel = [i for i in varSel if (not i.endswith("_F"))
                  and (not i.endswith("_z"))]
        for k in rVarSel.keys():
            if rVarSel[k] in repeats.keys():
                repeats[rVarSel[k]] = repeats[rVarSel[k]] + 1
                rVarSel[k] += "#"+str(repeats[rVarSel[k]])
            else:
                repeats[rVarSel[k]] = 0
        print(repeats)
        print(rVarSel)
        varSel.extend([val for _,val in rVarSel.items()])
        tempTC = [TC.copy() for TC in timeCourse]
        for i in range(len(tempTC)):
            tempTC[i] = tempTC[i].rename(columns=rVarSel)
        timeCourseVis = timeCourseVisualiser(tempTC)
        for i in range(math.ceil(len(varSel)/maxPerGrid)):
            varSel2 = varSel[i*maxPerGrid:(i+1)*maxPerGrid]
            if len(varSel2) == 0:
                break
            figPath = (source_name+"-"+pathInsert+"-timeCourse-"+clustor+
                       str(i)+".png")
            myImg = imgTag(figPath)
            cluTable[clustor].add(pathInsert,"dynamic "+str(i),myImg.print())
            figPath = os.path.join(working_directory, "figures",
                                   source_name, figPath)
            timeCourseVis.multiPlot(indexSelect=showList, varSelect=varSel2,
                                    save=figPath)
        tempTC = [TC.copy() for TC in timeCourseSS]
        for i in range(len(tempTC)):
            tempTC[i] = tempTC[i].rename(columns=rVarSel)
        timeCourseVis = timeCourseVisualiser(tempTC)
        for i in range(math.ceil(len(varSel)/maxPerGrid)):
            varSel2 = varSel[i*maxPerGrid:(i+1)*maxPerGrid]
            if len(varSel2) == 0:
                break
            figPath = (source_name+"-"+pathInsert+"-timeCourseSS-"+clustor+
                       str(i)+".png")
            myImg = imgTag(figPath)
            cluTable[clustor].add(pathInsert,"steady "+str(i),myImg.print())
            figPath = os.path.join(working_directory, "figures",
                                   source_name, figPath)
            timeCourseVis.multiPlot(indexSelect=showList, varSelect=varSel2,
                                    save=figPath)
    timeCourse = loadPick(["data",source_name, "mTorDep-"+source_name+"-"+
                           target+"-timeCourse.p"], relative=True)
    timeCourseSS = loadPick(["data",source_name, "mTorDep-"+source_name
                             +"-"+target+"-timeCourseSS.p"], relative=True)
    if (timeCourse is not None) and (timeCourseSS is not None):
        for PECase in showList:
            TC = [TC[PECase] for TC in timeCourse]
            timeCourseVis = timeCourseVisualiser(TC)
            varSel = list(TC[0].columns)
            varSel = [i for i in varSel if "TOR" in i and (not i.startswith("k")) 
                      and (not i.endswith("_F")) and (not i.endswith("_z"))]
            for i in range(math.ceil(len(TC[0].columns)/maxPerGrid)):
                varSel = varSel[i*maxPerGrid:(i+1)*maxPerGrid]
                if len(varSel) == 0:
                    break
                figPath = ("mTorDep-"+source_name+"-"+pathInsert
                           +"-timeCourseMTORdep"+str(i)+"-"+str(PECase)+
                           ".png")
                myImg = imgTag(figPath)
                depTable.add(pathInsert,"dynamic "+str(i)+"-"+str(PECase),
                             myImg.print())
                figPath = os.path.join(working_directory, "figures",
                                       source_name,
                                       figPath)
                timeCourseVis.multiPlot(indexSelect=showList,
                                        varSelect=varSel,
                                        save=figPath)
            TC = [TC[PECase] for TC in timeCourseSS]
            timeCourseVis = timeCourseVisualiser(TC)
            varSel = list(TC[0].columns)
            varSel = [i for i in varSel if "TOR" in i and (not i.startswith("k")) 
                      and (not i.endswith("_F")) and (not i.endswith("_z"))]
            for i in range(math.ceil(len(TC[0].columns)/maxPerGrid)):
                varSel = varSel[i*maxPerGrid:(i+1)*maxPerGrid]
                if len(varSel) == 0:
                    break
                figPath = ("mTorDep-"+source_name+"-"+pathInsert+
                           "-timeCourseSSMTORdep"+str(i)+"-"+str(PECase)+
                           ".png")
                myImg = imgTag(figPath)
                depTable.add(pathInsert,"steady "+str(i)+"-"+str(PECase),
                             myImg.print())
                figPath = os.path.join(working_directory, "figures",
                                       source_name,
                                       figPath)
                timeCourseVis.multiPlot(indexSelect=showList,
                                        varSelect=varSel,
                                        save=figPath)
    if isCombo:
        for key, value in RS["Combo-scale-"+source_name+"-"+
                             target.split("-")[0]].items():
            scaleTable.add(target,key,str(value))
    else:
        for key, value in RS["scale-"+source_name+"-"+
                             target.split("-")[0]].items():
            scaleTable.add(target,key,str(value))
    
    paramSelectList = breakSeriesByScale(GFID(params).iloc[
            list(range(BestClusterSize))].mean(), maxRunLength=10)
    for paramSelection, i in zip(paramSelectList,range(len(paramSelectList))):
        figPath = (source_name+"-"+pathInsert+ "-Violin-"+str(i)+".png")
        myImg = imgTag(figPath)
        figPath = os.path.join(working_directory, "figures", source_name,
                               figPath)
        PEVis.violinPlot(paramSelect=paramSelection,
                         RSSSelect=maxRSS,
                         save = figPath)
        violinTable.add(pathInsert, "paramiters "+str(i), myImg.print())
 
myStyle = """<style>
table {
  border-collapse: collapse;
}

table, th, td {
  border: 1px solid black;
}
</style>"""
           
myPage = myHTMLpage("report")
myPage.addToHead(myStyle)
myPage.addToBody(wTable.print())
myPage.addToBody(wbTable.print())
myPage.addToBody(torTable.print())
if len(depTable.df)!=0:
    myPage.addToBody(depTable.print())
for _, item in cluTable.items():
    myPage.addToBody(item.print())
myPage.addToBody(scaleTable.print())
violinTable.fillNA()
myPage.addToBody(violinTable.print())
figPath = os.path.join(working_directory, "figures", source_name,
                       "report.html")
text_file = open(figPath, "wt")
text_file.write(myPage.print())
text_file.close()