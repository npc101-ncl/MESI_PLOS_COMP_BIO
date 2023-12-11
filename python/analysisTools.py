#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 16:39:30 2019

@author: peter
"""

from sklearn.cluster import DBSCAN, MeanShift, estimate_bandwidth
from sklearn import preprocessing as pp
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd

def clusterParameterEstimation(PEDataFrame,testDistance=5,minInCluster=2):
    """Clusters estimated paramiters
    
    based on their 'nearness' in paramiter space using DBSCAN

    Args:
       PEDataFrame (DataFrame):  paramiter estimations to cluster

    Kwargs:
       testDistance (float): eps value used in DBSCAN. Basicly the
           minimum 'distance' alowed between clusters.
       minInCluster (int): minimum size alowed for cluster.

    Returns:
       Dict of Lists: Each cluster is returned as an intiger indexed
       List except for the unclustered points which are returned as a
       list in the -1 dict entry.
    """
    myScaledData = pp.StandardScaler().fit_transform(PEDataFrame.drop(
            columns="RSS", errors='ignore'))
    clustering = DBSCAN(eps=testDistance,
                        min_samples=minInCluster).fit(myScaledData)
    returnDict = {}
    for i in set(clustering.labels_):
        returnDict[str(i)]=[]
    for i in range(len(clustering.labels_)):
        returnDict[str(clustering.labels_[i])].append(i)
    return returnDict

def pcaParameterEstimation(df,components='mle'):
    """
    
    """
    myScaledData = pp.StandardScaler().fit_transform(df.drop(
            columns="RSS", errors='ignore'))
    clustering = PCA(n_components = components).fit(myScaledData)
    return clustering

def RSSClusterEstimation(PEDataFrame):
    """Clustered paramiter estimationes
    
    bassed purely on their RSS values using MeanShift.
    
    Args:
       PEDataFrame (DataFrame):  paramiter estimations to cluster
       
    Returns:
       Lists of Dicts: Clusteres are returned as dictionaries in order of
       the size of their largest RSS value (assending). dict structure is::
           
           {"id":arbitrary_id_assigned_to_cluster,
           "maxRSS":largest_RSS_value_in_cluster,
           "size":number_of_members_of_cluster}
           
    """
    myRSS = PEDataFrame.filter(items=["RSS"])
    myRSS = myRSS.replace([np.inf, -np.inf], np.nan)
    myRSS = myRSS.dropna()
    myScaledData = pp.StandardScaler().fit_transform(myRSS)
    bandwidth = estimate_bandwidth(myScaledData)
    if bandwidth==0:
        return [{"id":0, "maxRSS":PEDataFrame.iloc[0]['RSS'],
                 "size":len(PEDataFrame)}]
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(myScaledData)
    returnList = []
    for i in range(len(set(ms.labels_))):
        returnList.append({"id":i, "maxRSS":0, "size":0})
    for i in range(len(ms.labels_)):
        returnList[ms.labels_[i]]["maxRSS"] = max(
                returnList[ms.labels_[i]]["maxRSS"],
                PEDataFrame["RSS"].iloc[i])
        returnList[ms.labels_[i]]["size"]=returnList[ms.labels_[i]]["size"]+1
    returnList.sort(key=lambda cl:cl["maxRSS"])
    return returnList

def profileLikelyhood(TimeCourses,experamentMeans,StandardDeviations,
                      calculateConstant = True):
    """ gives -2*ln of profile likelyhood
    
    assuming normaly distributed noise given by the StandardDeviations data.
    TimeCourses, experamentMeans and StandardDeviations are 3 aligned lists
    who's enteries corisponde to each other.
    
    Args:
       TimeCourses (list of DataFrame or dict): A list each eliment of which
           can be a data frame or dict describing either a simulated time
           course or steady state. If time course must have
       experamentMeans (list of DataFrame or dict):
       StandardDeviations (list of DataFrame or dict):
           
    Kwargs:
       calculateConstant (bool): if false will not add constant based on SD to
       the weighted sum of square of residuals. If the same experamental data
       is used for two cases this constant will cancel out when you take the
       difrence anyway
       
    Returns:
       float: -2*ln of profile likelyhood
    """
    total = 0
    for TC, EM, SD in zip(TimeCourses, experamentMeans, StandardDeviations):
        if (isinstance(TC, pd.DataFrame) and isinstance(EM, pd.DataFrame)
            and isinstance(SD, pd.DataFrame)):
            if "Time" in EM.columns:
                times = list(EM["Time"])
            else: 
                times = list(EM.index)
            if "Time" in SD.columns:
                if times != list(SD["Time"]):
                    return None
            elif times != list(SD.index):
                return None
            cols = list(set(EM.columns).intersection(SD.columns))
            cols = [col for con in cols if col != "Time"]
            if len(cols)==0:
                return None
            elif not all([col in TC.columns for col in cols]):
                return None
            aligned = []
            for time in times:
                row = (TC[TC["Time"]==time])[cols].copy()
                if len(row)!=1:
                    print("Time course index match confusion")
                    return None
                else:
                    alignedTC.append(row.squeeze())
            aligned = pd.DataFrame(alignedTC, columns=cols)
            aligned = ((alignedTC-EM[cols])/SD[cols])**2
            total = total + aligned.sum().sum()
            if calculateConstant:
                total = total + np.log((2*np.pi*
                                        (SD[cols]**2)).prod().prod())
        elif (isinstance(TC, dict) and isinstance(EM, dict)
            and isinstance(SD, dict)):
            keys = list(set(EM.keys()).intersection(SD.keys()))
            if len(keys)==0:
                return None
            if not all([key in TC.keys() for key in keys]):
                return None
            total = total + sum([(TC[key]-EM[key])/SD[key] for key in keys])
            if calculateConstant:
                total = total + np.log(prod([2*np.pi*(SD[key]**2)
                                             for key in keys]))
        else:
            return None
    return total

def paramPCA(params,components='mle'):
    returnData = {}
    myArray = np.array(params.drop(columns="RSS"))
    print(myArray.shape)
    myPCA = PCA(n_components=components)
    myPCA = myPCA.fit(myArray)
    names = [col for col in params.columns if col!="RSS"]
    returnData["logVarDrop"] = np.diff(np.log(myPCA.explained_variance_))
    returnData["newPoints"] = myPCA.transform(myArray)
    returnData["newPoints"] = pd.DataFrame(returnData["newPoints"])
    returnData["vectors"] = [(pd.Series(data=comp, index=names)) for comp
              in myPCA.components_]
    returnData["cumVarExp"] = myPCA.explained_variance_ratio_.cumsum()
    returnData["PCA"] = myPCA
    return returnData