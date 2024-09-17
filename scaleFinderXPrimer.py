#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 11:12:39 2019

@author: npc101
"""

import os
import pickle
import pandas as pd
import time
import sys
from python.pycotoolsHelpers import *
from python.utilityTools import *

cmdLineArg = sys.argv[1:]

CLDict, CLSet = getCmdLineArgs()

working_directory = os.path.dirname(os.path.abspath(__file__))
RS={}
RS["is_interpolated"] = False
RS["is_data_constant_shifted"] = False
RS["mcf7_normed"] = "mcf7_normed" in CLSet

if len(cmdLineArg)>0:
    data_dir = os.path.join(working_directory, 'data', cmdLineArg[0])
    RS["run_name_prefix"]=cmdLineArg[0]
elif RS["is_interpolated"]:
    data_dir = os.path.join(working_directory, 'data', 'interpolatedX')
    RS["run_name_prefix"]="naiveFitX"
else:
    data_dir = os.path.join(working_directory, 'data', 'noninterpolatedX')
    RS["run_name_prefix"]="naiveFitX"
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)  

antimony_string =  """
    function MM(km, Vmax, S)
        Vmax * S / (km + S)
    end

    function MMWithKcat(km, kcat, S, E)
        kcat * E * S / (km + S)
    end

    function NonCompetitiveInhibition(km, ki, Vmax, n, I, S)
        Vmax * S / ( (km + S) * (1 + (I / ki)^n ) )
    end

    function NonCompetitiveInhibitionWithKcat(km, ki, kcat, E, n, I, S)
        kcat * E * S / ( (km + S) * (1 + (I / ki)^n ) )
    end

    function NonCompetitiveInhibitionWithKcatAndExtraActivator(km, ki, kcat, E1, E2, n, I, S)
        kcat * E1 * E2 * S / ( (km + S) * (1 + (I / ki)^n ) )
    end


    function MA1(k, S)
        k * S
    end

    function MA2(k, S1, S2)
        k * S1 * S2
    end

    function MA1Mod(k, S, M)
        k * S * M
    end

    function MA2Mod(k, S1, S2, M)
        k * S1 * S2 * M
    end

    function CompetitiveInhibitionWithKcat(km, ki, kcat, E, I, S)
        kcat * E * S / (km + S + ((km * I )/ ki)  )
    end

    function CompetitiveInhibition(Vmax, km, ki, I, S)
        Vmax * S / (km + S + ((km * I )/ ki)  )
    end

    function Hill(km, kcat, L, S, h)
        kcat * L * (S / km)^h  /   1 + (S / km)^h
    end

    model AktModelWithMMRateLaws()
        compartment             Cell = 1;
        var IRS1                in Cell;
        var IRS1a               in Cell;
        var IRS1pS636_639       in Cell;
        var Akt                 in Cell;
        var AktpT308            in Cell;
        var AktpS473            in Cell;
        var AktpT308S473        in Cell;
        var TSC2                in Cell;
        var TSC2pT1462          in Cell;
        var PRAS40              in Cell;
        var PRAS40pT246         in Cell;
        var S6K                 in Cell;
        var S6KpT389            in Cell;
        var FourEBP1            in Cell;
        var FourEBP1pT37_46     in Cell;
        var PI3K                in Cell;
        var pPI3K               in Cell;
        var pmTORC1            in Cell;
        var mTORC1cyt          in Cell;
        var mTORC1lys          in Cell;
        var pmTORC2            in Cell;
        var mTORC2             in Cell;       
        var mTORC1rap          in Cell;
        var rap                in Cell;
        const Insulin          in Cell;
        const AA               in Cell;
        const InsulinB         in Cell;
        
        IRS1_sum       := IRS1+IRS1a+IRS1pS636_639
        PI3K_sum       := PI3K+pPI3K
        Akt_sum        := Akt+AktpT308+AktpS473+AktpT308S473
        TSC2_sum       := TSC2+TSC2pT1462
        PRAS40_sum     := PRAS40+PRAS40pT246
        FourEBP1_sum   := FourEBP1+FourEBP1pT37_46
        S6K_sum        := S6K+S6KpT389
        
        Akt_wb             := Akt_sum
        AktpT308_wb        := AktpT308+AktpT308S473
        AktpS473_wb        := AktpS473+AktpT308S473
        PRAS40_wb          := PRAS40_sum
        PRAS40pT246_wb     := PRAS40pT246
        S6K_wb             := S6K_sum
        S6KpT389_wb        := S6KpT389
        TSC2_wb            := TSC2_sum
        TSC2pT1462_wb      := TSC2pT1462
        IRS1_wb            := IRS1_sum
        IRS1pS636_639_wb   := IRS1pS636_639
        FourEBP1_wb        := FourEBP1_sum
        FourEBP1pT37_46_wb := FourEBP1pT37_46

        // kinetic parameters
        kIRS1Act                = 0.1;
        kIRS1Inact              = 0.1;
        kIRS1Phos               = 0.1;
        kPI3KPhos               = 0.1;
        kPI3KDephos             = 0.1;
        kAktPhos_kcat           = 0.1;
        kAktDephos              = 0.1;
        
        // to be determined
        kAktByTor2Phos          = 0.1;
        kAktByTor2Dephos        = 0.1;
        
        kTSC2Phos               = 0.1;
        
        // to be determined
        kTSC2PhosBoost          = 2;
        
        kTSC2Dephos             = 0.1;
        kmTORC1cytToLys         = 0.1;
        kmTORC1LysToCyt         = 0.1;
        kmTORC1Phos             = 0.1;
        kmTORC1Dephos           = 0.1;
        kPras40Phos             = 0.1;
        kPras40Dephos           = 0.1;
        kFourEBP1Phos           = 0.1;
        kFourEBP1Dephos         = 0.1;
        kS6KPhos                = 0.1;
        kS6KDephos              = 0.1;
        
        // paramiters for candidate reactions
        kRAS                    = 0;
        kmTORC2Dephos           = 0.1;
        kmTORC2DephosByS6K      = 0;
        
        // raptor
        kmTORC1bind             = 1;
        kmTORC1unbind           = 1;

        // reactions // MMWithKcat(km, kcat, S, E)
        // function CompetitiveInhibitionWithKcat(km, ki, kcat, E, I, S)
        // function MM(km, Vmax, S)
        // function NonCompetitiveInhibitionWithKcat(km, ki, kcat, E, n, I, S)
        R1f     : IRS1 => IRS1a                         ; Cell * kIRS1Act*IRS1*(Insulin+InsulinB)
        R2f     : IRS1a => IRS1pS636_639                ; Cell * kIRS1Phos*IRS1a*S6KpT389
        R2b     : IRS1pS636_639 => IRS1                 ; Cell * kIRS1Inact*IRS1pS636_639
        R3f     : PI3K => pPI3K                         ; Cell * kPI3KPhos*PI3K*IRS1a
        R3b     : pPI3K => PI3K                         ; Cell * kPI3KDephos*pPI3K
        R4f     : Akt => AktpT308                       ; Cell * kAktPhos_kcat*Akt*pPI3K
        R4b     : AktpT308 => Akt                       ; Cell * kAktDephos*AktpT308
        R4p1f   : AktpS473 => AktpT308S473              ; Cell * kAktPhos_kcat*AktpS473*pPI3K
        R4p1b   : AktpT308S473 => AktpS473              ; Cell * kAktDephos*AktpT308S473
        R4p2f   : Akt => AktpS473                       ; Cell * kAktByTor2Phos*Akt*pmTORC2
        R4p2b   : AktpS473 => Akt                       ; Cell * kAktByTor2Dephos*AktpS473
        R4p3f   : AktpT308 => AktpT308S473              ; Cell * kAktByTor2Phos*AktpT308*pmTORC2
        R4p3b   : AktpT308S473 => AktpT308              ; Cell * kAktByTor2Dephos*AktpT308S473
        R5f     : TSC2 => TSC2pT1462                    ; Cell * kTSC2Phos*TSC2*AktpT308
        R5p1f   : TSC2 => TSC2pT1462                    ; Cell * kTSC2Phos*kTSC2PhosBoost*TSC2*AktpT308S473
        R5b     : TSC2pT1462 => TSC2                    ; Cell * kTSC2Dephos*TSC2pT1462
        R6f     : mTORC1cyt => mTORC1lys                ; Cell * kmTORC1cytToLys*mTORC1cyt*AA
        R6b     : mTORC1lys => mTORC1cyt                ; Cell * kmTORC1LysToCyt*mTORC1lys
        R7f     : mTORC1lys => pmTORC1                  ; Cell * kmTORC1Phos*mTORC1lys
        R7b     : pmTORC1 => mTORC1lys                  ; Cell * kmTORC1Dephos*TSC2*pmTORC1
        R8f     : PRAS40 => PRAS40pT246                 ; Cell * kPras40Phos*PRAS40*(AktpT308+AktpT308S473)
        R8b     : PRAS40pT246 => PRAS40                 ; Cell * kPras40Dephos*PRAS40pT246
        R9f     : FourEBP1 => FourEBP1pT37_46           ; Cell * kFourEBP1Phos*FourEBP1*pmTORC1
        R9b     : FourEBP1pT37_46 => FourEBP1           ; Cell * kFourEBP1Dephos*FourEBP1pT37_46
        R10f    : S6K => S6KpT389                       ; Cell * kS6KPhos*S6K*pmTORC1
        R10b    : S6KpT389 => S6K                       ; Cell * kS6KDephos*S6KpT389
        // candidte reactions
        R11f3   : mTORC2 => pmTORC2                     ; Cell * kRAS*mTORC2*(Insulin+InsulinB)
        R11b1   : pmTORC2 => mTORC2                     ; Cell * kmTORC2Dephos*pmTORC2
        R11b2   : pmTORC2 => mTORC2                     ; Cell * kmTORC2DephosByS6K*pmTORC2*S6KpT389
        
        R12f1   : mTORC1cyt + rap => mTORC1rap          ; Cell * kmTORC1bind*mTORC1cyt*rap
        R12f2   : mTORC1lys + rap => mTORC1rap          ; Cell * kmTORC1bind*mTORC1lys*rap
        R12f3   : pmTORC1 + rap => mTORC1rap            ; Cell * kmTORC1bind*pmTORC1*rap
        R12b1   : mTORC1rap => mTORC1cyt + rap          ; Cell * kmTORC1unbind*mTORC1rap
        
        // global variables
        Insulin = 1;
        InsulinB = 0;
        AA = 1;

        IRS1            = 10;
        IRS1a           = 0;
        IRS1pS636_639   = 0;
        PI3K            = 10;
        pPI3K           = 0;
        AktpT308S473    = 0;
        AktpS473        = 0;
        AktpT308        = 0;
        Akt             = 10;
        TSC2pT1462      = 10;
        TSC2            = 0;
        pmTORC1         = 0;
        mTORC1cyt       = 10;
        mTORC1lys       = 0;
        PRAS40          = 10;
        PRAS40pT246     = 0;
        FourEBP1pT37_46 = 0;
        FourEBP1        = 10;
        S6KpT389        = 0;
        S6K             = 10;
        pmTORC2         = 0;
        mTORC2          = 10;
        mTORC1rap       = 0;
        rap             = 0;

    end

    """
    
columns_from_data = {'IRS1pS636/639':'IRS1pS636_639_wb',
                     '4E-BP1':'FourEBP1_wb',
                     '4E-BP1pT37/46':'FourEBP1pT37_46_wb',
                     'Akt':'Akt_wb',
                     'PRAS40':'PRAS40_wb',
                     'S6K':'S6K_wb',
                     'TSC2':'TSC2_wb',
                     'IRS1':'IRS1_wb',
                     'AktpT308':'AktpT308_wb',
                     'AktpS473':'AktpS473_wb',
                     'PRAS40pT246':'PRAS40pT246_wb',
                     'S6KpT389':'S6KpT389_wb',
                     'TSC2pT1462':'TSC2pT1462_wb',
                     'FourEBP1pT37_46':'FourEBP1pT37_46_wb'}

def makeInterpolatedTable(df,expansion=3,kind="cubic"):
    newTime = [df.index.values[0]]
    for i in df.index.values[1:]:
        newTime.extend([newTime[-1]+j*(i-newTime[-1])/expansion
                        for j in range(1,expansion+1)])
    returnDF = pd.DataFrame(index=newTime)
    for (columnName, columnData) in df.iteritems():
        f = interp1d(df.index.values, columnData.values, kind=kind)
        returnDF[columnName]=f(returnDF.index.values)
    returnDF.index.name = df.index.name
    return returnDF

def applyScale(df,scale):
    if isinstance(df,list):
        df2 = []
        for dfs in df:
            df2.append(dfs.copy())
            for col, sVal in scale.items():
                if col in df2[-1].columns:
                    df2[-1][col] = df2[-1][col]*sVal
    else:
        df2 = df.copy()
        for col, sVal in scale.items():
            if col in df2.columns:
                df2[col] = df2[col]*sVal
    return df2
    
def importData(startsWith,fileName,columns,sheet_name='Sheet2',skiprows=51,
               index_col=0,usecols=21,nrows=12,addConstColls=None):
    if isinstance(usecols,int):
        myUsecols = list(range(usecols+1))
    else:
        myUsecols = usecols
    if fileName.endswith(".xlsx"):
        myEng = "openpyxl"
    else:
        myEng = None
    df = pd.read_excel (fileName,sheet_name=sheet_name,skiprows=skiprows,
                        index_col=index_col, engine=myEng,
                        usecols=myUsecols,nrows=nrows)
    temp = {col:col[:-2] for col in df.columns if col.endswith(".1")}
    df.rename(columns=temp, inplace=True)
    print(df.columns)
    df = df[df.index.str.startswith(startsWith)]
    
    df.index = [float(row[1])*60 for row in df.index.str.split(' ')]
    df.index.name = "Time"
    df.rename(columns=columns, inplace=True)
    discardCols = [col for col in df.columns if not col in
                   [val for _, val in columns.items()]]
    df=df.drop(discardCols, axis=1, errors='ignore')
    if isinstance(addConstColls,dict):
        for col, val in addConstColls.items():
            if not col in df.columns:
                df[col] = val
    return df

def genConstShift (dfForScale, dfForClash, relationsDict,
                   targetOveride = None):
    dfForScaleOut = dfForScale.copy()
    dfForClashOut = dfForClash.copy()
    for sumCol, partsList in relationsDict.items():
        valToAdj = dfForScaleOut.iloc[0][sumCol]
        if targetOveride is not None:
            if sumCol in targetOveride.keys():
                target = targetOveride[sumCol]
            else:
                target = 1
        else:
            target = 1
        dfForScaleOut[sumCol] = dfForScaleOut[sumCol]+target-valToAdj
        dfForScale[sumCol]=dfForScale[sumCol].apply(lambda x : 0
                  if x<0 else x)
        dfForClash[sumCol] = dfForClash[sumCol]+target-valToAdj
        dfForClash[sumCol]=dfForClash[sumCol].apply(lambda x : 0
                  if x<0 else x)
        for part in partsList:
            valToAdj = dfForScaleOut.iloc[0][part]
            if targetOveride is not None:
                if part in targetOveride.keys():
                    target = targetOveride[part]
                else:
                    target = 0.1
            else:
                target = 0.1
            dfForScaleOut[part] = dfForScaleOut[part]+target-valToAdj
            dfForScale[part]=dfForScale[part].apply(lambda x : 0
                      if x<0 else x)
            dfForClash[part] = dfForClash[part]+target-valToAdj
            dfForClash[part]=dfForClash[part].apply(lambda x : 0
                      if x<0 else x)
    return dfForScaleOut, dfForClash

def genScale2(dfForClash, dfForScale, relationsDict, scale=2,
              scaleTarget={},initialTarget={}):
    #each relation scaling problem broken down into its own dict entry
    problemChunks = {parientCol:{
            "clashTable":dfForClash[[parientCol]+childrenCol].copy(),
            "outTable":dfForScale[[parientCol]+childrenCol].copy()} 
            for parientCol, childrenCol in relationsDict.items()}
    # remaining non relational scailings in their own df
    remaindingCols = [colL for _, colL in relationsDict.items()]
    remaindingCols = [j for i in remaindingCols for j in i]
    remaindingCols.extend(relationsDict.keys())
    remaindingCols = [col for col in dfForScale.columns
                      if not col in remaindingCols]
    remaindingCols = dfForScale[remaindingCols].copy()
    # calculating max allowed scale
    for pCol, cCols in relationsDict.items():
        problemChunks[pCol]["maxScales"] = problemChunks[pCol][
                clashTable].copy()
        for cCol in cCols:
            problemChunks[pCol]["maxScales"][cCol]/=problemChunks[
                    pCol]["maxScales"][pCol]
        problemChunks[pCol]["maxScales"] = problemChunks[pCol][
                "maxScales"].drop([pCol], axis=1)
        problemChunks[pCol]["maxScales"] = problemChunks[pCol][
                "maxScales"].min()

def genScaleSub(colsToScale, targetCols=None, targetsSS=None,
                parentCols=None, limitColPairs=None, scale=2):  
    # promote to lists
    if isinstance(colsToScale,list):
        CTS = colsToScale
    else:
        CTS = [colsToScale]
    if isinstance(targetCols,list):
        TC = targetCols
    elif targetCols is not None:
        TC = [targetCols]
    if isinstance(targetsSS,list):
        TSS = targetsSS
    elif targetsSS is not None:
        TSS = [targetsSS]
    if isinstance(parentCols,list):
        PC = parentCols
    elif parentCols is not None:
        PC = [parentCols]
    if isinstance(limitColPairs,list):
        LCP = limitColPairs
    elif limitColPairs is not None:
        LCP = [limitColPairs]
    # calculate maximum allowed scale
    maxScale = []
    if parentCols is not None:
        if len(PC)!=len(CTS):
            print("parentCols not match colsToScale")
            return None
        for child, parent in zip(CTS, PC):
            maxScale.append((parent/child).min())
    if limitColPairs is not None:
        for pair in LCP:
            maxScale.append((pair["parent"]/pair["child"]).min())
    if len(maxScale)>0:
        maxScale = min(maxScale)
    else:
        maxScale = None
    # calculate scale based on targets
    top = 0
    bottom = 0
    if targetCols is not None:
        if len(TC)!=len(CTS):
            print("targetCols not match colsToScale")
            return None
        for col, targetCol in zip(CTS, TC):
            top += (col*targetCol).sum()
            bottom += (col*col).sum()
    if targetsSS is not None:
        for col, targetSS in zip(CTS, TSS):
            top += len(col)*(col.iloc[0])*(targetSS.iloc[-1])
            bottom += len(col)*col.iloc[0]*col.iloc[0]
    if bottom!=0:
        outScale = top/bottom
    elif parentCols is not None:
        # calculate scale based on child standard deviation
        mySD = pd.concat(CTS, ignore_index=True).std()
        myDist = float("inf")
        for child, parent in zip(CTS, PC):
            tempDist = (parent-child).min()
            if tempDist>=myDist:
                continue
            myDist = tempDist
            myIdx = (parent-child).idxmin()
            myMin = child[myIdx]
            myMax = parent[myIdx]
        outScale = (myMin-scale*mySD)/myMax
    else: 
        outScale = 1
    # apply scale ceiling to ensure child never larger than parent
    # (warning may give rise to negative scales)
    if maxScale is not None:
        outScale = min(maxScale,outScale)
    return outScale

def genScale3(dfForClash, dfForScale, relationsDict, scale=2,
              scaleTarget=None,initialTarget=None):
    outScale = {}
    dfFCList = listify(dfForClash)
    dfFSList = listify(dfForScale)
    STList = listify(scaleTarget)
    ITList = listify(initialTarget)
    for sumCol, partCols in relationsDict.items():
        for part in partCols:
            if STList is not None:
                tempTarg = [(lambda i: i[part] if part in i else None)(x)
                            for x in STList]
            else:
                tempTarg = None
            if ITList is not None:
                tempSSTarg = [(lambda i: i[part] if part in i else None)(x)
                              for x in ITList]
            else:
                tempSSTarg = None
            if dfFCList is not None:
                LCPs = [{"parent":i[sumCol], "child":i[part]}
                         for i in dfFCList]
            else:
                LCPs = None
            outScale[part] = genScaleSub([i[part] for i in dfFSList], 
                    targetCols=tempTarg, targetsSS=tempSSTarg,
                    parentCols=[i[sumCol] for i in dfFSList],
                    limitColPairs=LCPs, scale=scale)
            if outScale[part]<0:
                outScale[part] = 0
    return outScale

def genScale(dfForClash, dfForScale, relationsDict, scale=2,
             scaleTarget={},initialTarget={}):
    if scaleTarget == {}:
        scaleTarget1 = None
    else:
        scaleTarget1 = scaleTarget
    if initialTarget == {}:
        initialTarget1 = None
    else:
        initialTarget1 = initialTarget
    return genScale3(dfForClash, dfForScale, relationsDict, scale=scale,
                     scaleTarget=scaleTarget1,
                     initialTarget=initialTarget1)
    # overriding function with new code
    myScale = {}
    closestPoint = {}
    df2 = pd.DataFrame()
    for sumCol, partCols in relationsDict.items():
        if isinstance(partCols,list):
            df2 = pd.concat([df2, dfForClash[[sumCol]+partCols].copy()],
                             axis=1)
            for part in partCols:
                df2[part+"_distance"] = dfForClash[sumCol] - dfForClash[part]
                closestPoint[part]=df2[part+"_distance"].idxmin()
        else:
            df2 = pd.concat([df2, dfForClash[[sumCol]+[partCols]].copy()],
                             axis=1)
            df2[partCols+"_distance"] = (dfForClash[sumCol] -
               dfForClash[partCols])
            closestPoint[partCols] = df2[partCols+"_distance"].idxmin()
    sumMin={}
    partMax={}
    partSd={}
    for part, index in closestPoint.items():
        sumCol = [k for k, v in relationsDict.items() if part in v]
        sumCol = sumCol[0]
        sumMin[part]=df2[sumCol][index]
        partMax[part]=df2[part][index]
        partSd[part] = dfForScale[part].std()
        if (part in scaleTarget.keys()) and (part in initialTarget.keys()):
            myScale[part] = (((dfForScale[part]*scaleTarget[part]).sum()+
                   len(dfForScale[part])*(dfForScale[part].iloc[[0]]*
                      initialTarget[part]).sum())
            /((dfForScale[part]*dfForScale[part]).sum()+
              len(dfForScale[part])*(dfForScale[part].iloc[0]**2)))
            myScale[part] = min(sumMin[part]/partMax[part],myScale[part])
        elif part in scaleTarget.keys():
            myScale[part] = (((dfForScale[part]*scaleTarget[part]).sum())
            /((dfForScale[part]*dfForScale[part]).sum()))
            myScale[part] = min(sumMin[part]/partMax[part],myScale[part])
        elif scaleTarget!={}:
            raise RuntimeError('missing scale value')
        else:
            myScale[part] = (sumMin[part]-scale*partSd[part])/partMax[part]
        if myScale[part]<0:
            myScale[part] = 0
    return myScale

def genRunIn(df,working_directory,runName,i):
    dfPrequil = df.copy()
    dfPrequilB = df.iloc[:-1].copy()
    for col in dfPrequilB.columns:
        dfPrequilB[col]=df.iloc[0][col]
    timeOffset = list(df.index)[-1]
    dfPrequil.index = dfPrequil.index+timeOffset
    dfPrequil = pd.concat([dfPrequilB,dfPrequil])
    print(dfPrequil)
    fname4 = os.path.join(working_directory,
                          runName+'_MCF7_exp_PQ'+str(i)+'_data.csv')
    dfPrequil.to_csv(fname4)
    return timeOffset, fname4

def duplicator(antStr,cpfilep,paramToJoin = None):
    decompStr = getModelsAndFunctions(antStr)
    myModel = modelRunner(antString = antStr, run_dir = cpfilep)
    myE = myModel.getModelEliments()
    myP = myModel.extractModelParam()
    myPTJ = paramToJoin
    if paramToJoin is None:
        myPTJ = myE["kineticParams"] 
    elif not isinstance(myPTJ,list):
        myPTJ = []
    myModelName = decompStr["models"][0]["name"]
    notToJoin = [i for i in myE["kineticParams"]
                 if not i in myPTJ]
    myE["kineticParams"] = [i for i in myE["kineticParams"] 
                            if i in myPTJ]
    print(notToJoin)
    print(myE["kineticParams"])
    newLines = []
    for myMet in myE["metabolites"]:
        newLines.append("\tvar ModA_"+myMet+";")
        newLines.append("\tvar ModB_"+myMet+";")
    newLines.append("\tA: "+myModelName+"();")
    newLines.append("\tB: "+myModelName+"();")
    for myMet in myE["metabolites"]:
        newLines.append("\tA."+myMet+" is ModA_"+myMet+";")
        newLines.append("\tB."+myMet+" is ModB_"+myMet+";")
    for myAss in myE["assignments"]:
        newLines.append("\tA."+myAss+" is ModA_"+myAss+";")
        newLines.append("\tB."+myAss+" is ModB_"+myAss+";")
    for myPar in myE["kineticParams"]:
        newLines.append("\tA."+myPar+" is "+myPar+";")
        newLines.append("\tB."+myPar+" is "+myPar+";")
        newLines.append("\t"+myPar+" = "+myP[myPar]+";")
    for myPar in notToJoin:
        newLines.append("\tA."+myPar+" is ModA_"+myPar+";")
        newLines.append("\tB."+myPar+" is ModB_"+myPar+";")
    
    antStr = antStr +  "model combiner\n"
    antStr = antStr + "\n".join(newLines) + "\nend"
    return antStr

myUpperBound=10
myCopyNum=200
mySuperComputer=True
RS["suppressTSC2Data"] = False
RS["setUnmodReactToConst"] = "nonWB"
RS["overrideTotalsOnly"] = "ORTotalsOnly" in CLSet

RS["autoEstVar"] = "autoEst" in CLSet

if "estExPat" in CLDict:
    RS["estimationBlockPattern"] = CLDict["estExPat"]
if "estInPat" in CLDict:
    RS["estimationForcePattern"] = CLDict["estInPat"]
if "estRegExPat" in CLDict:
    RS["estimationRegExPattern"] = extractCompositArgument(
            CLDict["estRegExPat"])
if "estEx" in CLDict:
    RS["estimationBlockList"] = extractCompositArgument(CLDict["estEx"])
if "estIn" in CLDict:
    RS["estimationForceList"] = extractCompositArgument(CLDict["estIn"])

RS["antimony_string"] = antimony_string
RS["columns_from_data"] = columns_from_data

kineticParams = ["kIRS1Act", "kIRS1Inact", "kIRS1Phos", "kPI3KPhos",
                 "kPI3KDephos", "kAktPhos_kcat",
                 "kAktDephos", "kTSC2Phos", "kTSC2Dephos",
                 "kmTORC1cytToLys", "kmTORC1LysToCyt", "kmTORC1Phos",
                 "kmTORC1Dephos", "kPras40Phos", "kPras40Dephos",
                 "kFourEBP1Phos", "kFourEBP1Dephos", "kS6KPhos",
                 "kS6KDephos", "kAktByTor2Phos", "kAktByTor2Dephos",
                 "kTSC2PhosBoost", "kRAS", "kmTORC2Dephos"]

IRS1SynthParam = ["kIRS1Synth"]

kineticParamsOpt = ["kmTORC2DephosByS6K"]

unmodReactNonWB = ["kmTORC1LysToCyt","kmTORC1Phos","kPI3KDephos"]
unmodReactWB = ["kIRS1Inact","kPras40Dephos","kFourEBP1Dephos","kS6KDephos",
                "kAktDephos","kAktByTor2Dephos"]
if RS["suppressTSC2Data"]:
    unmodReactNonWB.append("kTSC2Dephos")
else:
    unmodReactWB.append("kTSC2Dephos")
if RS["setUnmodReactToConst"]=="nonWB":
    for name in unmodReactNonWB:
        kineticParams.remove(name)
elif RS["setUnmodReactToConst"]=="all":
    for name in unmodReactNonWB:
        kineticParams.remove(name)
    for name in unmodReactWB:
        kineticParams.remove(name)    
        
fullTORmod = ["kmTOR_RIC_Comb", "kmTORC2_Dis", "kmTOR_RPT_Comb",
              "kmTORC1_Dis", "kmTORC1LysToCyt", "kmTORC1cytToLys", 
              "kmTORC1Phos", "kmTORC1Dephos", "kRAS", "kmTORC2Dephos"]

if "fullTOR" in CLSet:
    kineticParams = list(set(kineticParams+fullTORmod))
    kineticParams.remove("kTSC2PhosBoost")
    
kineticBGRParams = ["kBGR_IRS1Phos", "kBGR_PI3KPhos", "kBGR_AktPhos_kcat_1",
                    "kBGR_AktPhos_kcat_2", "kBGR_AktByTor2Phos_1", 
                    "kBGR_AktByTor2Phos_2", "kBGR_TSC2Phos",
                    "kBGR_mTORC1Dephos", "kBGR_Pras40Phos",
                    "kBGR_FourEBP1Phos", "kBGR_S6KPhos"]
    
metabolites= ["IRS1", "IRS1a", "IRS1pS636_639", "PI3K", "pPI3K",
              "AktpT308", "Akt", "AktpS473", "AktpT308S473", "TSC2pT1462",
              "TSC2", "pmTORC1", "mTORC1cyt", "mTORC1lys", "PRAS40",
              "PRAS40pT246", "FourEBP1pT37_46", "FourEBP1", "S6KpT389",
              "S6K", "mTORC2", "pmTORC2", "InsulinB", "AAB", "RPTOR", 
              "RICTOR"]


# determins how inital conditions are force set from western blot data
# in paramiter estimations
overrideDict = {'IRS1_wb':'IRS1_sum',
                'IRS1pS636_639_wb':'IRS1pS636_639',
                'FourEBP1_wb':'FourEBP1_sum',
                'FourEBP1pT37_46_wb':'FourEBP1pT37_46',
                'Akt_wb':'Akt_sum',
                'AktpT308_wb':'AktpT308_wb',
                'AktpS473_wb':'AktpS473_wb',
                'PRAS40_wb':'PRAS40_sum',
                'PRAS40pT246_wb':'PRAS40pT246',
                'S6K_wb':'S6K_sum',
                'S6KpT389_wb':'S6KpT389',
                'TSC2_wb':'TSC2_sum',
                'TSC2pT1462_wb':'TSC2pT1462'}
if RS["suppressTSC2Data"]:
    del overrideDict["TSC2pT1462_wb"]
    
if RS["overrideTotalsOnly"]:
    overrideDict = {k:v[:-3]+"T" for k,v in overrideDict.items() if v.endswith("_sum")}
    
RS["kineticParams"] = kineticParams
RS["kineticParamsOpt"] = kineticParamsOpt
RS["kineticBGRParams"] = kineticBGRParams
RS["IRS1SynthParam"] = IRS1SynthParam
RS["metabolites"] = metabolites
RS["overrideDict"] = overrideDict
RS["intiTime"] = time.time()
RS["data_dir"] = data_dir

if RS["mcf7_normed"]:
    FNS = "norm"
else:
    FNS = "data"
    
data_filenames = {"t47d":os.path.join(working_directory,
                                      't47d_'+FNS+'.xlsx'),
                  "zr75":os.path.join(working_directory,
                                      'zr75_'+FNS+'.xlsx')}

RS["data_filenames"] = data_filenames
   
tempDict = {'Akt_wb': ['AktpT308_wb', 'AktpS473_wb'],
            'PRAS40_wb': ['PRAS40pT246_wb'],
            'S6K_wb': ['S6KpT389_wb'],
            'TSC2_wb': ['TSC2pT1462_wb'],
            'IRS1_wb': ['IRS1pS636_639_wb'],
            'FourEBP1_wb': ['FourEBP1pT37_46_wb']}

if RS["suppressTSC2Data"]:
    tempDict.pop('TSC2_wb',None)
    
RS["data_relations"] = tempDict

if __name__ == "__main__":
    if "ant" in CLDict.keys():
        RS["antimony_string"] = loadTxt(CLDict["ant"], relative=True)
    elif len(cmdLineArg)>1:
        RS["antimony_string"] = loadTxt(cmdLineArg[1], relative=True)
    
    """
    f = open(resolvePath([data_dir,'refAntStr.txt']), "w")
    f.write(RS["antimony_string"])
    f.close()
    """
    savePick([data_dir,'runSwitches.p'], RS)