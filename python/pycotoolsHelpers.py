#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 14:00:40 2019

@author: peter
"""

import site, os, time, re
from pycotools3 import model, tasks, viz
import tellurium as te
import pandas as pd
import logging, pickle, sys
from random import randint
from subprocess import getoutput
from math import gcd
import gc


def addCopasiPath(copasiPath):
    """adds the path to copasi to pythons working path
    
    Copasi is often not in the path python uses. This will check if its
    not and add it.

    Args:
       copasiPath (str):  a string with the path to copasi on the
       relivent machine.
    """
    if not re.search(copasiPath, os.environ["PATH"]):
        os.environ["PATH"] += os.pathsep + copasiPath
        
def antStrFromCopasi(path):
    myModel = model.Model(path)
    return myModel.to_antimony()
        
def extractAntReactions(antString):
    antLines=antString.splitlines()
    reactions = []
    for antLine in antLines:
        reaction={}
        temp=antLine.split("->")
        if len(temp) != 2:
            temp=antLine.split("=>")
            if len(temp) != 2:
                continue
            else:
                reaction['irreversible']=True
        else:
            reaction['irreversible']=False
        temp2=temp[0].split(":")
        if len(temp2)!=2:
            if len(temp2)==1:
                LHS=temp2[0].strip()
            else:
                continue
        else:
            reaction['name']=temp2[0].strip()
            LHS=temp2[1].strip()
        temp2=temp[1].split(";")
        if len(temp2)==1:
            RHS=temp2[0].strip()
        else:
            reaction['formula']=temp2[1].strip()
            RHS=temp2[0].strip()
        LHS=[i.strip() for i in LHS.split("+")]
        if LHS==['']:
            LHS=[]
        for i in range(len(LHS)):
            if LHS[i][0]=='$':
                LHS[i]={'fixed':True, 'var':LHS[i][1:]}
            else:
                LHS[i]={'fixed':False, 'var':LHS[i]}
        reaction['LHS']=LHS
        RHS=[i.strip() for i in RHS.split("+")]
        if RHS==['']:
            RHS=[]
        for i in range(len(RHS)):
            if RHS[i][0]=='$':
                RHS[i]={'fixed':True, 'var':RHS[i][1:]}
            else:
                RHS[i]={'fixed':False, 'var':RHS[i]}
        reaction['RHS']=RHS
        reactions.append(reaction)
    return reactions

def getModelsAndFunctions(antStr):
    myFuncts = re.findall(r'(?<=\bfunction\b)(.*?)(?=\bend\b)', antStr,
                          flags=re.DOTALL)
    myFuncts = [re.split(r'\bfunction\b', i)[-1] for i in myFuncts]
    myModels = re.findall(r'(?<=\bmodel\b)(.*?)(?=\bend\b)', antStr,
                          flags=re.DOTALL)
    myFunctsLines = [myFun.splitlines() for myFun in myFuncts]
    myModelsLines = [myMod.splitlines() for myMod in myModels]
    myFunctNames = [(lambda x: x.group(1) if x is not None
                     else None)(re.match(r'[\s*]*(\w+)',myFun[0]))
                    for myFun in myFunctsLines]
    myModelNames = [(lambda x: x.group(1) if x is not None
                     else None)(re.match(r'[\s*]*(\w+)',myMod[0]))
                    for myMod in myModelsLines]
    for i, name in enumerate(myFunctNames):
        arguments = re.search(r'\((.*?)\)',myFunctsLines[i][0]).group(1)
        arguments = arguments.split(",")
        arguments = [A.strip() for A in arguments]
        myFuncts[i] = {"lines":myFunctsLines[i], "name":name,
                "text":"function"+myFuncts[i]+"end",
                "arguments":arguments}
    for i, name in enumerate(myModelNames):
        myModels[i] = {"lines":myModelsLines[i], "name":name,
                "text":"model"+myModels[i]+"end"}
    return {"functions":myFuncts, "models":myModels}

def subFunIntoReact (reaction,functionList):
    rNew = reaction["formula"]
    for F in functionList:
        subs = []
        formula = F["lines"][1].split(";")[0].strip()
        for M in re.finditer(F["name"]+r'\s*\((.*?)\)', rNew):
            args = M.group(1)
            args = [A.strip() for A in args.split(",")]
            if len(args)!=len(F["arguments"]):
                continue
            subs.append({"s":M.start(),"e":M.end(),"a":args})
            "text[:M.start()] + newText + text[M.end():]"
        for S in reversed(subs):
            ins = formula
            rep = {A:B for A, B in zip(F["arguments"], S["a"])}
            pattern = re.compile("|".join(rep.keys()))
            ins = pattern.sub(lambda m: rep[m.group(0)], ins)
            rNew = rNew[:S["s"]] + "(" +ins +")"+ rNew[S["e"]:]
    return rNew

def replaceVariable(theString,oldName,newName):
    """changes the name of a variable in an antimony string
    
    It's helpfull to be able to rename variables in antimony strings on
    the fly to utalise pycotools 'prefix' feature. This quick a dirty
    methiod uses regular expresions to do this and works suprisingly
    reliably.

    Args:
       theString (str):  the antimony string to parse
       oldName (str):  name of variable to replace
       newName (str):  name of variable to be the replacement

    Returns:
       str: the new modified antimony string.
    """
    return newName.join(re.split('(?<![\w_])'+oldName+'(?![\w_])',
                                 theString))
    
def renameCSVColumns(CSVPathOrList,targetDir,renameDict,indepToAdd=None):
    """preproceses CSV data files for paramiter estimations
    
    Allows us to take one of more CSV files and rename their columns (for
    example to take advndage of prefix feture) and optionaly add new
    columns bearing the "_indep" sufix used to indicate a column should
    be treated as an indipendent variable in the paramiter estimation.
    
    Args:
       CSVPathOrList (str or list of str):  path or list of paths to csv
           files to modify.
       targetDir (str):  path to directory in which modified csv files
           should be saved.
       renameDict (dict):  dictioary of str key value pairs that indicate
           the replacments to perform ineach csv files column names.
       
    Kwargs:
       indepToAdd (dict or list of dicts): key value pairs representing
           the independed variables that should be added to each data set
           and the corisponding values they should be set to. Ordering of 
           list should match CSVPathOrList.
           
    Returns:
       List or str: path or list of paths to saved modified files.
       (ordering will match CSVPathOrList)
    """
    if not isinstance(CSVPathOrList, list):
        CSVPathOrList = [CSVPathOrList]
    if indepToAdd is not None:
        if not isinstance(indepToAdd, list):
            indepToAdd = [indepToAdd]
        if len(CSVPathOrList)!=len(indepToAdd):
            print("each data set needs independed variables if you use them")
            return None
    outList=[]
    for CSVPath, i in zip(CSVPathOrList,range(len(CSVPathOrList))):
        CSVData = pd.read_csv(CSVPath)
        CSVData=CSVData.rename(columns=renameDict)
        if indepToAdd is not None:
            for myKey, myVal in indepToAdd[i].items():
                CSVData[myKey+"_indep"]=myVal
        newPath = os.path.join(targetDir,os.path.split(CSVPath)[1])
        CSVData.to_csv(newPath,index=False)
        outList.append(newPath)
    if len(outList)==1:
        outList=outList[0]
    return outList
  
class modelRunner:
    def __init__(self, antString=None, run_dir=None, CopasiFile = None):
        if run_dir is None:
            # run_dir = os.path.dirname(os.path.abspath(''))
            run_dir = os.path.abspath('')
            run_dir = os.path.join(run_dir, 'copasiRuns')
            i=1
            nameFree=False
            while not nameFree:
                temp_dir = os.path.join(run_dir, 'run'+str(i))
                if not os.path.isdir(temp_dir):
                    os.makedirs(temp_dir)
                    self.run_dir=temp_dir
                    nameFree=True
                else:
                    i=i+1
        else:
            if os.path.isdir(run_dir):
                self.run_dir=run_dir
            else:
                os.makedirs(run_dir)
                self.run_dir=run_dir
        if antString is None:
            if CopasiFile is None:
                self.antString = '''
                model empty_model()
                end
                '''
            else:
                myModel = model.Model(CopasiFile)
                self.antString = myModel.to_antimony()
        else:
            self.antString=antString
        self.methodDict = {
                "particle_swarm_test":{
                        "method":"particle_swarm",
                        "swarm_size":5,
                        "iteration_limit":100},
                "particle_swarm_default":{
                        "method":"particle_swarm",
                        "swarm_size":50,
                        "iteration_limit":2000},
                "particle_swarm_rigorous":{
                        "method":"particle_swarm",
                        "swarm_size":100,
                        "iteration_limit":4000},
                "particle_swarm_aggressive":{
                        "method":"particle_swarm",
                        "swarm_size":200,
                        "iteration_limit":6000},
                "particle_swarm_heroic":{
                        "method":"particle_swarm",
                        "swarm_size":300,
                        "iteration_limit":6000},
                "particle_swarm_ridiculous":{
                        "method":"particle_swarm",
                        "swarm_size":350,
                        "iteration_limit":7000}}
                
    def checkRuningProcesses(self, job=None, user=None):
        path = self.run_dir
        allparts = []
        while True:
            parts = os.path.split(path)
            if parts[1] == "":  # sentinel for absolute paths
                break
            elif parts[1] == path: # sentinel for relative paths
                allparts.insert(0, parts[1])
                break
            else:
                path = parts[0]
                allparts.insert(0, parts[1])
        path = "_"+"_".join(allparts)+"_sge_job_file.sh"
        if isinstance(job,str):
            path = " -n " + job
        else:
            if not isinstance(user,str):
                path = " -n " + path
            else:
                path = ""
        if isinstance(user,str):
            path = path + " -u " + user
        comand = "squeue" + path + " | wc -l"
        try:
            return int(getoutput(comand))-1
        except:
            if job is None:
                return float('inf')
            else:
                return 0
        
    def extractModelParam(self):
        """a way of getting the model paramiters / initial conditions out of a
        model you've loaded for refrence
        
        Returns:
            dict: keys are paramiter / variable names values are model values
        """
        copasi_filename = self.genPathCopasi("extractor")
        self.recentModel = model.loada(self.antString, copasi_filename)
        return self.recentModel.parameters.copy().squeeze().to_dict()
    
    def getReactions(self):
        copasi_filename = self.genPathCopasi("extractor")
        self.recentModel = model.loada(self.antString, copasi_filename)
        myReact = self.recentModel.reactions
        myReact = [{"name":R.name, "expression":R.expression,
                    "rate_law":R.rate_law, "reversible":R.reversible,
                    "simulation_type":R.simulation_type}
                   for R in myReact]
        myReact = pd.DataFrame(myReact)
        myFun = self.recentModel.functions
        myFun = [{"name":F.name, "key":F.key, "expression":F.expression}
                 for F in myFun]
        lines = self.antString.splitlines()
        for i, F in enumerate(myFun):
            myRegEx = re.compile(r"^\s*function\s+"+F["name"]+
                                 "\s*\(([\s,\w])\)\s*$")
            dec = [myRegEx.search(L) for L in lines]
            if i ==0 :
                print(F["name"])
                print(lines)
            dec = [M.group() for M in dec if M is not None]
            if len(dec)==1:
                myFun[i]["arguments"]=dec[0]
        myFun = pd.DataFrame(myFun)
        rr = te.loada(self.antString)
        myTR = [rr[rid] for rid in rr.getReactionIds()]
        return {"reactions":myReact, "functions":myFun}
                    
    def getModelEliments(self):
        copasi_filename = self.genPathCopasi("extractor")
        self.recentModel = model.loada(self.antString, copasi_filename)
        tEl = self.recentModel.get_variable_names(which='a',
                                                  include_assignments=True)
        mEl = self.recentModel.get_variable_names(which='m',
                                                  include_assignments=False)
        gEl = self.recentModel.get_variable_names(which='g', 
                                                  include_assignments=False)
        lEl = self.recentModel.get_variable_names(which='l', 
                                                  include_assignments=False)
        cEl = self.recentModel.get_variable_names(which='c', 
                                                  include_assignments=False)
        tEl = list(set(tEl)-set(mEl)-set(gEl)-set(lEl)-set(cEl))
        kEl = list(set(gEl) | set(lEl))
        return({"metabolites":mEl, "kineticParamsGlobal":gEl,
                "kineticParamsLocal":lEl, "kineticParams":kEl,
                "compartments":cEl, "assignments":tEl})
            
    def genPathCopasi(self,nameBase,suffix=".cps"):
        """generates a copasi path that isn't being used yet
        
        we want to avoid acidently useing the same copasi file for difrent
        calculations as it seems they don't always get overwriten. This
        function generates a path to an unused file name in the classes
        run directory by trying difrent sufixes.
        
        Args:
           nameBase (str):  base to use in trying to creat file name.
           
        Kwargs:
           suffix (str): optional string for file extention incase we ever
               want to generate a path to non copasi file someday. defaults
               to ".cps"
               
        Returns:
           str: path to unused filename in objects run directory.
        """
        i=0
        nameFree=False
        while not nameFree:
            copasi_filename = os.path.join(self.run_dir,nameBase+
                                           str(i)+suffix)
            nameFree = not os.path.exists(copasi_filename)
            i=i+1
        return copasi_filename
    
    def clearRunDirectory(self):
        """empties the objects run directory
        
        empties the objects run directory of all '.cps', '.txt', '.sbml',
        '.csv' extention files and any empty directories. It seems
        pycotools creats a lot of files in paramiter estimations it has a
        habit of reusing. Clearing the run directory is one way to prevent
        this.
        """
        for root, dirs, files in os.walk(self.run_dir, topdown=False):
            for name in files:
                if name.lower().endswith(('.cps', '.txt', '.sbml', '.csv')):
                    os.remove(os.path.join(root, name))
            for name in dirs:
                if len(os.listdir(os.path.join(root, name)))==0:
                    os.rmdir(os.path.join(root, name))
            
    def genPrefixAntString(self,estimatedVar,prefix="_"):
        """Creates a new antimony string with added prefixs
        
        creats a new antimony string (stored seperatly from the refrence
        string) where certain variables in the string have had a prefix
        added. Pycotools uses prefixes to determine which variables should
        be estimated in a parameter estimation.
        
        Args:
           estimatedVar (list of str):  the variables in the antimony
               string (and hence the model) that should be prefixed.
               
        Kwargs:
           prefix (str): optionaly the prefix to add
        """
        self.prefixAntString = self.antString
        for name in estimatedVar:
            self.prefixAntString = replaceVariable(self.prefixAntString,
                                                   name,prefix+name)
            
    def genReactionAntString(self, revTag = "RevRe__",
                             iRevTag = "IrRevRe__"):
        """Creates a new antimony string with added prefixs
        
        creates a new antimony string (stored seperatly from the refrence
        string) with aditional asigned variables to represent each reaction
        in the model. These variables can be used to create diagrams showing
        the 'flow' of reactions.
        
        Kwargs:
           revTag (str):  the prefix used to represent reactions with
               revesable flow.
           iRevTag (str):  the prefix used to represent reactions with
               irevesable flow.
        """
        
        lines = self.antString.splitlines()
        lines = [line.split("#")[0] for line in lines]
        rLines = [line.split(":") for line in lines if
                  len(line.split(":"))==2]
        rLines = [[line[0]]+line[1].split(";") for line in rLines
                  if len(line[1].split(";"))>=2]
        rLines = [[part.strip() for part in line] for line in rLines]
        rLines = [line for line in rLines if ("->" in line[1]) or
                  ("=>" in line[1])]
        rLines = [[line[0], "->" in line[1], line[2]] for line in rLines]
        rLines = [[revTag+line[0], line[1], line[2]] if line[1] else
                  [iRevTag+line[0], line[1], line[2]] for line in rLines]
        rLines = [line[0]+" := "+line[2]+";" for line in rLines]
        primed = False
        for i, line in zip(range(len(lines)),lines):
            if line.strip().startswith("model"):
                primed = True
            if (line.strip() == "end") and primed:
                break
        print("line "+str(i))
        indent = ""
        while indent == "" and i>0:
            i = i-1
            indent = re.search(r'^\s*', lines[i]).group()
        rLines = [indent+line for line in rLines]
        self.reactionAntString = "\n".join(lines[:i+1]+rLines+lines[i+1:])
        
    def genColumnWeights(self,df,method = "Mean Square"):
        """Atempts to replicate the weights copasi generates for experamental
        data in calibration.
        
        Args:
           df (DataFrame or Series):  Results of a spicific experament (time
             course or steady state)
               
        Kwargs:
           method (str): methiod used to calculate weights, curently only
             suports 'Mean Square'
             
        Returns:
           Serise: the weights for the experamental file used for calculating
             RSS values
        """
        if method == "Mean Square":
            if isinstance(df,pd.Series):
                df2 = pd.DataFrame([df.to_dict()])
            else:
                df2 = df.copy()
            if "Time" in df2.columns:
                df2 = df2.drop(columns="Time")
            df2 = df2*df2
            df2 = df2.mean(axis = 0)
            df2 = 1/df2
            df2 = df2/(df2.sum())
            return df2
        return None
    
    def getRSSforTimeCourse(self, TC, df, colWeights=None,
                            method = "Mean Square"):
        """Atempts to calculate the RSS associated with a time course in the
        same way copasi does in its paramiter estimations.
        
        Args:
           TC (DataFrame):  The time course to calculated RSS for
           df (DataFrame):  Results of a asociated experament (time course)
               
        Kwargs:
           colWeights (Serise):  used to overide the weights as produced by
             genColumnWeights
           method (str): methiod used to calculate weights, curently only
             suports 'Mean Square'
        
        Returns:
           float: the RSS value.
        """
        if colWeights is None:
            colWeights = self.genColumnWeights(df, method)
        mySum = 0
        for col, wei in colWeights.items():
            subSum=0
            for _, row in df.iterrows():
                TCR = TC[TC["Time"] == row["Time"]].squeeze()
                subSum = subSum + (row[col]-TCR[col])**2
            mySum = mySum + wei*subSum
        return mySum
    
    def genRefCopasiFile(self, filePath = None, adjustParams = None,
                         setIndex = None, sbmlPath = None):
        """generates a copasi file from the model optionaly adjusted by a
        paramiter estimation. Useful in testing things in copasi before
        coding.
        
        Kwargs:
           filePath (str): path to where you want the copasi file. if you omit
             is created in model run directory
           adjustParams (DataFrame): paramiters (probably from a paramiter
             estimation) to be used to overide the models value. (must also
             have setIndex)
           setIndex (int): index of the paramiters in adjustParams you wantto
             use to overide the models paramiters.
        """
        if filePath is None:
            copasi_filename = self.genPathCopasi("refFile")
        else:
            copasi_filename = filePath
        
        self.recentModel = model.loada(self.antString, copasi_filename)
        
        if isinstance(adjustParams,pd.DataFrame): 
            if (setIndex>=0 and setIndex < len(adjustParams)):
                model.InsertParameters(self.recentModel,
                                       df=adjustParams,
                                       index=setIndex,inplace=True)
        if sbmlPath is not None:
            self.recentModel.to_sbml(sbml_file=sbmlPath)
    
    def runSteadyStateFinder_TC(self, params=None, rocket=False, duration=1,
                                genReactions=False, cycleMax=1,
                                threshold=0):
        workingParams = params
        for i in range(cycleMax): 
            timeCourses = self.runTimeCourse(duration,
                                             adjustParams=workingParams,
                                             rocket=rocket,
                                             genReactions=genReactions)
            oldWorkingParams = workingParams
            if isinstance(timeCourses,list):
                workingParams = pd.DataFrame([timeCourse.iloc[-1]
                                              for timeCourse
                                              in timeCourses])
                if ((workingParams-oldWorkingParams)**2
                    ).sum(axis=0).min()>threshold:
                    break
            else:
                acum = 0
                workingParams = timeCourses.iloc[-1].to_dict()
                for key, value in oldWorkingParams.items():
                    acum = acum + (workingParams[key]-value)**2
                if acum>threshold:
                    break
        if isinstance(timeCourses,list):            
            return [row.to_dict() for _, row in workingParams.iterrows()]
        else:
            return workingParams
        
    # paralel version broken
    def runSteadyStateFinder(self, params=None, rocket=False, chunks=1):
        if isinstance(params, pd.DataFrame) and rocket:
            list_df = [params[i:i+chunks] for i
                       in range(0, params.shape[0], chunks)]
            jobName = "SS"+str(randint(0, 99999))
            for index, df in zip(range(len(list_df)),list_df):
                file = open(os.path.join(self.run_dir,
                                         'SteadyState'+str(index)+'.p'), 'wb')
                pickle.dump({"param":df, "run_dir":self.run_dir,
                             "antimony_string":self.antString},file)
                file.close()
                myScriptName = self.genPathCopasi("SSSlurmScript",
                                                      suffix = ".sh")
                shellString = ("#!/bin/bash\n "+
                               '#SBATCH --job-name="'+jobName+'"\n'+
                               "python -u " + __file__ + 
                               " runSteadyStateFinder" +
                               ' SteadyState'+str(index)+'.p')
                f = open(myScriptName, 'w')
                f.write(shellString)
                f.close()
                os.system("sbatch "+myScriptName)
            while self.checkRuningProcesses(job = jobName)>0:
                time.sleep(60)
            df2 = pd.DataFrame()
            for index in range(len(list_df)):
                file = open(os.path.join(self.run_dir,
                                         'SteadyState'+str(index)+'.p'), 'rb')
                df3 = pickle.load(file)
                file.close()
                if isinstance(df3, pd.DataFrame):
                    df2 = pd.concat([df2,df3])
                else:
                    for _ in range(list_df[index].shape[0]):
                        df2.append(pd.Series(), ignore_index=True)
            return df2
        elif isinstance(params, pd.DataFrame):
            dictList = [self.runSteadyStateFinder(params=row.to_dict())
                        for _, row in params.iterrows()]
            return pd.DataFrame(dictList)
        r = te.loada(self.antString)
        outValues = r.getFloatingSpeciesIds()
        boundValues = r.getBoundarySpeciesIds()
        outValues.extend(boundValues)
        r.conservedMoietyAnalysis = True
        r.setSteadyStateSolver('nleq1')
        r.getSteadyStateSolver().maximum_iterations = 1000
        if isinstance(params,dict):
            for key, value in params.items():
                temp = r[key]
                r[key] = value
                # print(key,": ",temp," -> ",r[key])
        try:
            r.steadyState()
        except:
            pass
        rState = {key:r[key] for key in outValues}
        return rState
    
    # if you reuse the same PEName twise the results get over writen by
    # the first
    def runParamiterEstimation(self,expDataFP,PEName=None,
                               copyNum=1,estimateIC=True,estimatedVar=None,
                               prefix='_',rocket=False,upperParamBound=None,
                               upperParamBoundOverride=None,
                               lowerParamBound=None,
                               lowerParamBoundOverride=None,
                               method="particle_swarm_default",
                               overrideParam=None,
                               indepToAdd=None, endTime=None, chain=None,
                               randStartVal = True, throttle = None):
        """runs a paramiter estimations on the loaded model
        
        Performs a paramiter estimation based on the model loaded
        
        Args:
           expDataFP (str or list of str): either a path to a file or a list
           of paths where each file represents a distinct experament to be 
           calibrated against. Should be CSV files and if a "Time" column is
           not defined only the 1st row will be used and treated as a steady
           state.
           
        Kwargs:
           indepToAdd (dict or list of dicts): key value pairs representing
           the independed variables that should be assumed for each
           experament. Used for overriding kinetic paramiters or initial
           conditions. Ordering of list should match expDataFP.
           PEName (str): provide a distinct name for the calculation if you
           don't want to use the auto assigned one.
           copyNum (int): the number of copys of the paramiter estimation to
           do.
           estimateIC (bool): whether the initial condtions (eg metabolites)
           should be estimated.
           estimatedVar (list of str): if provided defines the variables that
           should be estimated.
           prefix (str): can be used to avoid clashes with the naming
           conventions you use for variables / paramiters. Modify the prefix
           to ensure adding this charicter as a symbol wont create ambiguity in
           your names.
           rocket (bool): determines weather the the function should paralise
           the task and submit it as seperate slurm jobs.
           upperParamBound (float): the upper limit alowed for paramiters in
           the search.
           lowerParamBound (float): the lower bound for paramiters in the
           search.
           method (str or dict): a string that represents the methiod to use
           in paramiter estimation. Optionaly pass a dictionary instead where
           the key : value pairs are passed to the pycotools function
           context.set to set the methiod.
           overrideParam (dict or DataFrame): if a dict overrides the value of
           paramiters / metabolites etc in acordence with the dicts entries
           prior to computation. If a DataFrame will split the estimation into
           multipul estimations for each row and overide the values for each
           estimation based on that row. (so number of estimations is
           len(overrideParam)*copyNum) This option can be used for
           profilelikelyhoods or to do the second part in a global chaser
           estimation if randStartVal is set to False.
           randStartVal (bool): should each estimation have the value of the
           variables to be estimated randomised prior to estimation? this will
           erase any value set by overrideParam but not indepToAdd.
           endTime (float): a time (seconds since epoch) that the calculation
           should end prematirly at if reached.
           chain (int): don't use this, its used internaly to help the
           function call itself when spliting up paramiter estimation tasks
           when overrideParam is a DataFrame.
               
        Returns:
           dict: a dict with one eliment which is a DataFrame of results. this
           format is for consistancy with pycotools. If overrideParam was a
           DataFrame the row order should be 1st overide copy 1, ..., 1st
           overide copy n, ..., last overide copy 1, ..., last overide copy n. 
        """
        
        if isinstance(overrideParam, pd.DataFrame):
            return self.chainRPEhandeler(expDataFP,PEName,copyNum,
                                         estimateIC,estimatedVar,
                                         prefix,rocket,upperParamBound,
                                         lowerParamBound,method,
                                         overrideParam,indepToAdd,
                                         endTime,randStartVal,throttle)
        if chain is None:
            if PEName is None:
                copasi_filename = self.genPathCopasi("paramiterEstimation")
            else:
                copasi_filename = os.path.join(self.run_dir, PEName)
        else:
            if PEName is None:
                copasi_filename = self.genPathCopasi("paramiterEstimation")
            else:
                copasi_filename = os.path.join(self.run_dir,
                                               str(chain).join(os.path.splitext(PEName)))
        
        overlap = list(set(estimatedVar).intersection(set([j for i 
                       in [list(k.keys()) for k 
                           in indepToAdd] for j in i])))
        if len(overlap)>0:
            print("You my not estimate independent variables:",*overlap)
            return None
        if estimatedVar is not None:
            if chain is None or chain == 0:
                self.genPrefixAntString(estimatedVar)
            self.recentModel = model.loada(self.prefixAntString,
                                           copasi_filename)
            colRenameDict = {i:(prefix+i) for i in estimatedVar}
            if chain is None or chain == 0:
                expDataFP = renameCSVColumns(expDataFP,self.run_dir,
                                             colRenameDict,
                                             indepToAdd=indepToAdd)
        else:
            self.recentModel = model.loada(self.antString, copasi_filename)
        if overrideParam is not None:
            rOverrideParam = {(lambda x: prefix+x if x 
                               in estimatedVar else x)(k):v for k,v
                              in overrideParam.items()}
            model.InsertParameters(self.recentModel,
                                   parameter_dict=rOverrideParam,
                                   inplace=True)
        
        if estimateIC:
            PEParams='gm'
        else:
            PEParams='g'
            
        if rocket==True:
            runMode='slurm'
        else:
            runMode=True
   
        with tasks.ParameterEstimation.Context(self.recentModel, expDataFP, 
                                               context='s', 
                                               parameters=PEParams) as context:
            context.set('randomize_start_values', randStartVal)
            context.set('separator', ',')
            if estimatedVar is not None:
                context.set('prefix', prefix)
            if upperParamBound is not None:
                context.set('upper_bound', upperParamBound)
            if lowerParamBound is not None:
                context.set('lower_bound', lowerParamBound)
            if isinstance(method,dict):
                for key, value in method.items():
                    context.set(key, value)
            else:
                for key, value in self.methodDict[method].items():
                    context.set(key, value)
            context.set('run_mode', runMode) 
            context.set('pe_number', 1) 
            context.set('copy_number', copyNum)
            config = context.get_config()
            if not randStartVal:
                for key, value in rOverrideParam.items():
                    if key in config.items.fit_items.keys():
                        config.items.fit_items[key].start_value = value
            if isinstance(upperParamBoundOverride,dict):
                for key, value in upperParamBoundOverride.items():
                    if key in config.items.fit_items.keys():
                        config.items.fit_items[key].upper_bound = value
            if isinstance(lowerParamBoundOverride,dict):
                for key, value in lowerParamBoundOverride.items():
                    if key in config.items.fit_items.keys():
                        config.items.fit_items[key].lower_bound = value
        if isinstance(throttle,dict):
            if ("user" in throttle.keys() and
                "jobLimit" in throttle.keys()):
                while (self.checkRuningProcesses(user = throttle["user"])
                >= throttle["jobLimit"]):
                    time.sleep(60)
        self.recentPE = tasks.ParameterEstimation(config)
        if chain is not None:
            return {"chain":chain, "expDataFP":expDataFP,
                    "PE":self.recentPE}
        #logging.disable(logging.WARNING)
        while rocket and self.checkRuningProcesses()>0:
            time.sleep(60)
            if endTime is not None and rocket:
                if endTime<=time.time():
                    print("out of time")
                    break
        try:
            parse_object = viz.Parse(self.recentPE)
            if not isinstance(parse_object[list(parse_object)[0]],
                                           pd.DataFrame):
                print("failiour to parse (unexpected return structure)")
                return None
        except:
            print("failiour to parse (throws unknown error)")
            return None      
        #logging.disable(logging.NOTSET)
        if estimatedVar is not None:
            colRenameDict = {(prefix+i):i for i in estimatedVar}
        return_data={}
        for key in list(parse_object):
            if estimatedVar is None:
                return_data[key]=parse_object[key].copy()
            else:
                return_data[key]=parse_object[key].copy()
                return_data[key].rename(columns=colRenameDict, inplace=True)
            if overrideParam is not None:
                for newParam in overrideParam.keys():
                    if not newParam in return_data[key].columns:
                        return_data[key][newParam]=overrideParam[newParam]
        return return_data
    
    def chainRPEhandeler(self, expDataFP, PEName, copyNum, estimateIC,
                         estimatedVar, prefix, rocket, upperParamBound,
                         lowerParamBound, method, overrideParam, indepToAdd,
                         endTime, randStartVal, throttle):
        i = 0
        workingExpDataFP = expDataFP
        k=[]
        for index, rowOP in overrideParam.iterrows():
            j = self.runParamiterEstimation(workingExpDataFP,
                                            PEName=PEName,
                                            copyNum=copyNum,
                                            estimateIC=estimateIC,
                                            estimatedVar=estimatedVar,
                                            prefix=prefix,
                                            rocket=rocket,
                                            upperParamBound=upperParamBound,
                                            lowerParamBound=lowerParamBound,
                                            method=method,
                                            overrideParam=rowOP.to_dict(),
                                            indepToAdd=indepToAdd,
                                            chain=i,
                                            randStartVal=randStartVal,
                                            throttle=throttle)
            workingExpDataFP = j["expDataFP"]
            j["overrideParam"] = rowOP
            i = i+1
            k.append(j)
        while rocket and self.checkRuningProcesses()>0:
            time.sleep(60)
            if endTime is not None and rocket:
                if endTime<=time.time():
                    print("out of time")
                    break
        parse_object_list = []
        for myIt in k:
            try:
                df = viz.Parse(myIt["PE"])
                df = df[next(iter(df))]
                if isinstance(df,pd.DataFrame):
                    if len(df)<copyNum:
                        print("parse incompleat")
                        df = pd.concat([df,pd.DataFrame(
                                [{} for _ in range(copyNum-len(df))])],
                                       ignore_index=True)
                else:
                    print("parse bad return")
                    df = pd.DataFrame([{} for _ in range(copyNum)])
            except:
                print("parse throw")
                df = pd.DataFrame([{} for _ in range(copyNum)])
            parse_object_list.append(df)
        print(parse_object_list)
        if estimatedVar is not None:
            colRenameDict = {(prefix+i):i for i in estimatedVar}
        return_data=[]
        for parse_object, myIt in zip(parse_object_list, k):
            df = parse_object.copy()
            if estimatedVar is not None:
                df.rename(columns=colRenameDict, inplace=True)
            for newParam, value in myIt["overrideParam"].items():
                    if not newParam in df.columns:
                        df[newParam]=value
            return_data.append(df)
        return_data = pd.concat(return_data, ignore_index=True)
        try:
            return {next(iter(viz.Parse(k[0]["PE"]))):return_data}
        except:
            return {" ":return_data}
    
    def preProcessParamEnsam(self,Params):
        """fills in the NAs in a data frame based on model
        
        Fills in NAs in a data frame by matching the column names to the named
        variables / paramiters in the antimony string and using those values
        to do the filling in.
        
        Args:
           Params (DataFrame): A DataFrame sutable for the overrideParam
           variable in runParamiterEstimation except that it has NA values in
           it.
               
        Returns:
           DataFrame: data frame with NAs filled in.
        """
        copasi_filename = self.genPathCopasi("adjuster")
        self.recentModel = model.loada(self.antString, copasi_filename)
        adjust = {}
        for colName in list(Params.columns):
            baseVal = self.recentModel.get('metabolite', colName, by='name')
            if baseVal == []:
                baseVal = self.recentModel.get('global_quantity', colName,
                                               by='name')
                if baseVal!=[]:
                    baseVal = baseVal.to_df()
                    adjust[colName] = baseVal["Value"]["initial_value"]
            else:
                baseVal = baseVal.to_df()
                adjust[colName] = baseVal["Value"]["concentration"]
        return Params.fillna(value=adjust)
    
    def TCendState(self, timeCourse, variables = None):
        if isinstance(timeCourse,list):
            myTC=timeCourse
        else:
            myTC=[timeCourse]
        myTC = [i.tail(1).squeeze().to_dict() for i in myTC]
        myTC = pd.DataFrame(myTC)
        myTC = myTC.drop(columns=["Time"])
        if variables == "metabolites":
            copasi_filename = self.genPathCopasi("getter")
            self.recentModel = model.loada(self.antString, copasi_filename)
            myVar = [i.name for i in self.recentModel.metabolites]
        if variables is not None:
            myTC = myTC[myVar]
        return myTC
    
    def runTimeCourse(self,duration,stepSize=None,intervals=None,
                      TCName=None,adjustParams=None,subSet=None,
                      rocket=False,genReactions=False,throttle=None):
        """Run one or more time courses (simulations)
        
        Takes the model and runs a simulation. Useful for getting time course
        data for paramiter ensabels so you can see how difrent members of the
        ensambel behave and match difrently to the calibration data.
        
        Args:
           duration (float or list of float): duration to simulate over. (list
           can be used if multipul simulations being done with difrent
           durations)
           
        Kwargs:
           stepSize (float): time step size between returned data points. does
           not effect stap size used in ODE solution.
           intervals (int): over all number of time steps that should be in
           return data.
           TCName (str): a Name for the calculation if you don't want to use
           an automaticly generated default.
           adjustParams (DataFrame): A dataframe such as those outputted by
           runParamiterEstimation that overides the paramiter / metabolite
           values for the simulations
           subSet (list of ints): a list of which members of adjustParams
           should be used for simulation.
           rocket (bool): determines weather the the function should paralise
           the task and submit it as seperate slurm jobs.
           genReactions (bool): if true will return data for reaction rates as
           well.
               
        Returns:
           DataFrame or list of DataFrames: the time course as a data frame of
           if multipul time courses a list of dataframes each its own
           timecourse.
        """
        if (stepSize is None) and (intervals is None):
            stepSize = 0.01
        if adjustParams is not None:
            if subSet is None:
                subSet = range(len(adjustParams.index))
            results = []
            if not isinstance(duration,list):
                duration = [duration for i in subSet]
            if len(duration)!=len(subSet):
                return False
            if rocket:
                jobName = "TC"+str(randint(0, 99999))
                timeCourses = []
                for setIndex, myDur in zip(subSet,duration):
                    myScriptName = os.path.join(self.run_dir,
                                                jobName+str(setIndex))
                    if not os.path.isdir(myScriptName):
                        os.mkdir(myScriptName)
                    copasi_filename = os.path.join(myScriptName,
                                                   "timeCourse.cps")
                    if genReactions!=False:
                        if isinstance(genReactions,dict):
                            self.genReactionAntString(*genReactions)
                        else:
                            self.genReactionAntString()
                        self.recentModel = model.loada(self.reactionAntString,
                                                       copasi_filename)
                    else:
                        self.recentModel = model.loada(self.antString,
                                                       copasi_filename)
                    model.InsertParameters(self.recentModel,
                                           df=adjustParams,
                                           index=setIndex,inplace=True)
                    if stepSize is not None:
                        self.recentTimeCourse = tasks.TimeCourse(
                                self.recentModel,end=myDur,
                                step_size=stepSize,
                                run=False)
                    else:
                        self.recentTimeCourse = tasks.TimeCourse(
                                self.recentModel,end=myDur,
                                step_size=myDur/intervals,
                                run=False)
                    timeCourses.append(self.recentTimeCourse)
                    myScriptName = os.path.join(myScriptName,
                                                jobName + ".sh")
                    shellString = ("#!/bin/bash\nCopasiSE "+
                                   copasi_filename)
                    f = open(myScriptName, 'w')
                    f.write(shellString)
                    f.close()
                    if isinstance(throttle,dict):
                        if ("user" in throttle.keys() and
                            "jobLimit" in throttle.keys()):
                            while(self.checkRuningProcesses(
                                    user=throttle["user"]) >=
                            throttle["jobLimit"]):
                                time.sleep(60)
                    os.system("sbatch "+myScriptName)
                while(self.checkRuningProcesses(job=jobName+".sh") > 0):
                    time.sleep(60)
                for theTimeCourse in timeCourses:
                    sucsessful=False
                    logging.disable(logging.WARNING)
                    while not sucsessful:
                        try:
                            parse_object = viz.Parse(theTimeCourse)
                            results.append(parse_object.data.copy())
                            sucsessful = True
                        except:
                            time.sleep(60)
                    logging.disable(logging.NOTSET)
                return results
            else:
                if TCName is None:
                    copasi_filename = self.genPathCopasi("timeCourse")
                else:
                    copasi_filename = os.path.join(self.run_dir, TCName)
                if genReactions!=False:
                    if isinstance(genReactions,dict):
                        self.genReactionAntString(*genReactions)
                    else:
                        self.genReactionAntString()
                    self.recentModel = model.loada(self.reactionAntString,
                                                   copasi_filename)
                else:
                    self.recentModel = model.loada(self.antString,
                                                   copasi_filename)
                for setIndex, myDur in zip(subSet,duration):
                    model.InsertParameters(self.recentModel,
                                           df=adjustParams,
                                           index=setIndex,inplace=True)
                    if stepSize is not None:
                        self.recentTimeCourse = tasks.TimeCourse(
                                self.recentModel,end=myDur,
                                step_size=stepSize)
                    else:
                        self.recentTimeCourse = tasks.TimeCourse(
                                self.recentModel,end=myDur,
                                step_size=myDur/intervals)
                    results.append(viz.Parse(
                            self.recentTimeCourse).data.copy())
            return results
        else:
            if TCName is None:
                copasi_filename = self.genPathCopasi("timeCourse")
            else:
                copasi_filename = os.path.join(self.run_dir, TCName)
                
            if genReactions!=False:
                if isinstance(genReactions,dict):
                    self.genReactionAntString(*genReactions)
                else:
                    self.genReactionAntString()
                self.recentModel = model.loada(self.reactionAntString,
                                               copasi_filename)
            else:
                self.recentModel = model.loada(self.antString,
                                               copasi_filename)
            if stepSize is not None:
                self.recentTimeCourse = tasks.TimeCourse(
                        self.recentModel,end=duration,
                        step_size=stepSize)
            else:
                self.recentTimeCourse = tasks.TimeCourse(
                        self.recentModel,end=duration,
                        step_size=duration/intervals)
            return viz.Parse(self.recentTimeCourse).data.copy()
        
    def getRSSforParamiters(self, param, IC, exp, method = "Mean Square",
                            step = None, GRRSName = None, rocket=False):
        if not isinstance(IC,list):
            IC = [IC]
        if not isinstance(exp,list):
            exp = [exp]
        if len(IC)!=len(exp):
            return None
        superParam = []
        for case in IC:
            df = param.copy()
            for key, val in case.items():
                df[key]=val
            superParam.append(df)
        superParam = pd.concat(superParam, sort=True, ignore_index=True)
        if "RSS" in superParam.columns:
            superParam = superParam.drop(columns="RSS")
        superParam = self.preProcessParamEnsam(superParam)
        duration = 0
        myStep = []
        pointStats = []
        for path in exp:
            df = pd.read_csv(path)
            pointStats.append({"len":len(df),
                               "cols":[i for i in df.columns
                                       if i != "Time"]})
            duration = max(duration,df.iloc[-1]["Time"])
            myStep.append(df["Time"].diff().min(skipna=True))
        if step is not None:
            myStep = step
        else:
            myStep = min(myStep)
        timeCourses = self.runTimeCourse(duration, stepSize=myStep,
                      TCName = GRRSName, adjustParams=superParam,
                      rocket=rocket)
        RSS = []
        for i, TC in zip(range(len(timeCourses)),timeCourses):
            if (i%len(param) == 0):
                df = pd.read_csv(exp[i//len(param)])
                myColWeights = self.genColumnWeights(df, method = method)
                pointStats[i//len(param)]["cols"] = len(
                        set(pointStats[i//len(param)]["cols"]).intersection(
                                set(myColWeights.keys())))
                pointStats[i//len(param)] = (pointStats[i//len(param)][
                        "cols"]*pointStats[i//len(param)]["len"])
                # need to add weights to experaments based on data points etc
            if TC.iloc[-1]["Time"] == duration:
                expRSS = self.getRSSforTimeCourse(TC, df,
                                                  colWeights=myColWeights,
                                                  method = method)
                RSS.append(expRSS)
            else:
                RSS.append(None)
        RSS2 = []
        pointStats = [i/sum(pointStats) for i in pointStats]
        for i in range(len(param)):
            subSum = 0
            for j in range(len(IC)):
                if RSS[i+j*len(param)] is None:
                    subSum = None
                    break
                subSum = subSum + RSS[i+j*len(param)]*pointStats[j]
            RSS2.append(subSum)
        return RSS2
    
    def calculateLikelyhood(param, IC, exp, method = "Mean Square"):
        pass
    
    def makeSensAnalDF(self, sensa, params=None):
        if isinstance(params,pd.DataFrame):
            df = params.copy()
            if "RSS" in df.columns:
                df.drop(columns="RSS")
        else:
            df = pd.DataFrame()
        for col in sensa.keys():
            if not (col in df.columns):
                df[col] = pd.NA
        df = self.preProcessParamEnsam(df)
        myList = [df]
        for col, adjs in sensa.items():
            if isinstance(adjs,int):
                df2 = df[df[col]!=0].copy()
                df2[col] = df2[col]*(1+adjs)
                myList.append(df2)
                df2 = df[df[col]!=0].copy()
                df2[col] = df2[col]*(1-adjs)
                myList.append(df2)
            elif isinstance(adjs,list):
                for adj in adjs:
                    df2 = df[df[col]!=0].copy()
                    df2[col] = df2[col]*(1+adj)
                    myList.append(df2)
                    df2 = df[df[col]!=0].copy()
                    df2[col] = df2[col]*(1-adj)
                    myList.append(df2)
        df = pd.concat(myList, ignore_index=True)
        return df 

    def getSensativitys(self, adjustParams=None, setIndex=None):
        if (setIndex is None) and (isinstance(adjustParams,pd.DataFrame)):
            indexs = range(len(adjustParams))
        elif not isinstance(setIndex,list):
            indexs = [setIndex]
        else:
            indexs = setIndex
        if not isinstance(adjustParams,pd.DataFrame):
            indexs = None
        if indexs is not None:
            outList = []
            for index in indexs:
                copasi_filename = self.genPathCopasi("sensativity")
                self.recentModel = model.loada(self.antString,
                                               copasi_filename)
                model.InsertParameters(self.recentModel, df=adjustParams,
                                       index=index, inplace=True)
                try:
                    self.recentSensativity = tasks.Sensitivities(
                            self.recentModel)
                    outList.append(
                            self.recentSensativity.sensitivities.copy())
                except:
                    outList.append(None)
            return outList
        else:
            copasi_filename = self.genPathCopasi("sensativity")
            self.recentModel = model.loada(self.antString,
                                           copasi_filename)
            try:
                self.recentSensativity = tasks.Sensitivities(
                        self.recentModel)
                return self.recentSensativity.sensitivities.copy()
            except:
                return None
    
    def runProfileLikelyhood2(self, expDataFP, adjRange, estimatedVar,
                         estimateIC=True, prefix='_', rocket=False,
                         upperParamBound=None, lowerParamBound=None,
                         overrideParam=None, indepToAdd=None, depth=1,
                         method=None):   
        """runs a profile likelyhood on the loaded model
        
        Performs a profile likelyhood based on the loaded model by running
        multipul paramiter estimations where one paramiter is not estimated
        but instead adjusted so see how changing it leeds to other changes in
        the estimated variables. The methiod used is default hooke_jeeves.
        
        Args:
           expDataFP (str or list of str): either a path to a file or a list
           of paths where each file represents a distinct experament to be 
           calibrated against. Should be CSV files and if a "Time" column is
           not defined only the 1st row will be used and treated as a steady
           state.
           adjRange (list of float): the range of adjustments you want applied
           (by scale) to the default value being adjusted. should include 1.0
           for comtrole.
           estimatedVar (list of str): if provided defines the variables that
           should be re-estimated.
           
        Kwargs:
           indepToAdd (dict or list of dicts): key value pairs representing
           the independed variables that should be assumed for each
           experament. Used for overriding kinetic paramiters or initial
           conditions. Ordering of list should match expDataFP.
           estimateIC (bool): whether the initial condtions (eg metabolites)
           should be estimated.
           prefix (str): can be used to avoid clashes with the naming
           conventions you use for variables / paramiters. Modify the prefix
           to ensure adding this charicter as a symbol wont create ambiguity in
           your names.
           rocket (bool): determines weather the the function should paralise
           the task and submit it as seperate slurm jobs.
           upperParamBound (float): the upper limit alowed for paramiters in
           the search.
           lowerParamBound (float): the lower bound for paramiters in the
           search.
           overrideParam (dict): overrides the value of paramiters /
           metabolites etc in acordence with the dicts entries prior to
           computation.
           method (dict): can be used to overide usual paramiter estimation
           method settings.
               
        Returns:
           tupul: first of which is a dict of data frames where the key is the
           paramiter fixed and adjusted and the dataframe contains the
           paramiter estimation results. The 2nd in tupel is the refrence for
           the origional values the 3rd is legasy and can be ignored.
        """
        if estimatedVar is None:
            return None
        extParam = self.extractModelParam()
        returnData = {}
        returnRef = {}
        self.arrayPE = []
        probeArray = []
        estVar3 = {}
        for adjVar in estimatedVar:
            adjList = []
            estVar2 = [var for var in estimatedVar if var!=adjVar]
            colRenameDict = {i:(prefix+i) for i in estVar2}
            self.genPrefixAntString(estVar2)
            estVar3[adjVar]=estVar2
            for adjustment in adjRange:
                copasi_filename = self.genPathCopasi("profileLikelyhood")
                self.recentModel = model.loada(self.prefixAntString,
                                               copasi_filename)
                expDataFP = renameCSVColumns(expDataFP,self.run_dir,
                                             colRenameDict, 
                                             indepToAdd=indepToAdd)
                if isinstance(overrideParam,dict):
                    overrideParamB = {key:value for key,value
                                      in overrideParam.items()
                                      if key!="RSS"}
                else:
                    overrideParamB = {}
                if adjVar in overrideParamB.keys():
                    returnRef[adjVar] = overrideParamB[adjVar]
                else:
                    returnRef[adjVar] = float(extParam[adjVar])
                overrideParamB[adjVar] = returnRef[adjVar]*adjustment
                overrideParamB = {(lambda k: prefix+k if k
                                   in estVar2 else k)(key):val for key, val 
                                  in overrideParamB.items()}
                model.InsertParameters(self.recentModel,
                                       parameter_dict = overrideParamB,
                                       inplace = True)
                
                if estimateIC:
                    PEParams='gm'
                else:
                    PEParams='g'
                    
                if rocket == True:
                    runMode='slurm'
                else:
                    runMode=True
                
                with tasks.ParameterEstimation.Context(self.recentModel,
                                                       expDataFP, 
                                                       context='s', 
                                                       parameters=PEParams) as context:
                    context.set('randomize_start_values', False)
                    context.set('separator', ',')
                    if estimatedVar is not None:
                        context.set('prefix', prefix)
                    if upperParamBound is not None:
                        context.set('upper_bound', upperParamBound)
                    if lowerParamBound is not None:
                        context.set('lower_bound', lowerParamBound)
                    if isinstance(method,dict):
                        for mKey, mVal in method.items():
                            context.set(mKey, mVal)
                    else:
                        context.set('method', 'hooke_jeeves')
                    context.set('run_mode', runMode) 
                    context.set('pe_number', 1) 
                    context.set('copy_number', depth) 
                    config = context.get_config()
                    self.arrayPE.append(tasks.ParameterEstimation(config))
                probeArray.append({"adjVar":adjVar,
                                   "adjustment":adjustment,
                                   "overrideParamB":overrideParamB[adjVar]})
                gc.collect()
        logging.disable(logging.WARNING)
        if rocket:
            while self.checkRuningProcesses()>0:
                time.sleep(60)
        """
        print("enters adj loop")
        i=0
        for PE, myProbe in zip(self.arrayPE,probeArray):
            notDone = True
            print("searching",myProbe["adjVar"],":",myProbe["adjustment"],":",
                  "number", i)
            i=i+1
            while notDone:
                try:
                    parse_object = viz.Parse(PE)
                    parse_object2 = parse_object[list(parse_object)[0]]
                    if isinstance(parse_object2, pd.DataFrame):
                        notDone = len(parse_object2) < 1
                except:
                    time.sleep(60)
            gc.collect()
        """
        logging.disable(logging.NOTSET) 
        print("exits adj loop")
        for PE, probe in zip(self.arrayPE, probeArray):
            try:
                parse_object = viz.Parse(PE)
                parse_object2 = parse_object[list(parse_object)[0]]
                df = parse_object2.copy()
                if not probe["adjVar"] in df.columns:
                    df[probe["adjVar"]] = probe["overrideParamB"]
                if not probe["adjVar"] in returnData.keys():
                    returnData[probe["adjVar"]]=[]
                returnData[probe["adjVar"]].append(df)
            except:
                print("missing data for", probe)
            gc.collect()
        returnData = {adjVar:pd.concat(adjList,ignore_index=True)
                      for adjVar, adjList in returnData.items()}
        for adjVar, adjList in returnData.items():
            colRenameDict = {(prefix+i):i for i in estVar3[adjVar]}
            adjList.rename(columns=colRenameDict, inplace=True)
        return returnData, returnRef, {}
    
    def runProfileLikelyhood(self, expDataFP, adjRange, estimatedVar,
                             estimateIC=True, prefix='_', rocket=False,
                             upperParamBound=None, lowerParamBound=None,
                             overrideParam=None, indepToAdd=None, depth=1,
                             method=None, throttle=None):   
        """runs a profile likelyhood on the loaded model
        
        Performs a profile likelyhood based on the loaded model by running
        multipul paramiter estimations where one paramiter is not estimated
        but instead adjusted so see how changing it leeds to other changes in
        the estimated variables. The methiod used is default hooke_jeeves.
        
        Args:
           expDataFP (str or list of str): either a path to a file or a list
           of paths where each file represents a distinct experament to be 
           calibrated against. Should be CSV files and if a "Time" column is
           not defined only the 1st row will be used and treated as a steady
           state.
           adjRange (list of float): the range of adjustments you want applied
           (by scale) to the default value being adjusted. should include 1.0
           for comtrole.
           estimatedVar (list of str): if provided defines the variables that
           should be re-estimated.
           
        Kwargs:
           indepToAdd (dict or list of dicts): key value pairs representing
           the independed variables that should be assumed for each
           experament. Used for overriding kinetic paramiters or initial
           conditions. Ordering of list should match expDataFP.
           estimateIC (bool): whether the initial condtions (eg metabolites)
           should be estimated.
           prefix (str): can be used to avoid clashes with the naming
           conventions you use for variables / paramiters. Modify the prefix
           to ensure adding this charicter as a symbol wont create ambiguity in
           your names.
           rocket (bool): determines weather the the function should paralise
           the task and submit it as seperate slurm jobs.
           upperParamBound (float): the upper limit alowed for paramiters in
           the search.
           lowerParamBound (float): the lower bound for paramiters in the
           search.
           overrideParam (dict): overrides the value of paramiters /
           metabolites etc in acordence with the dicts entries prior to
           computation.
           method (dict): can be used to overide usual paramiter estimation
           method settings.
               
        Returns:
           tupul: first of which is a dict of data frames where the key is the
           paramiter fixed and adjusted and the dataframe contains the
           paramiter estimation results. The 2nd in tupel is the refrence for
           the origional values the 3rd is legasy and can be ignored.
        """
        if estimatedVar is None:
            return None
        extParam = self.extractModelParam()
        returnData = {}
        returnRef = {}
        sendPaths = {}
        returnPaths = {}
        scriptPaths = {}
        jobName = "PE"+str(randint(0, 99999))
        for adjVar in estimatedVar:
            if isinstance(overrideParam,dict):
                overrideParamB = {key:value for key,value
                                  in overrideParam.items()
                                  if key!="RSS"}
            else:
                overrideParamB = {}
            if adjVar in overrideParamB.keys():
                returnRef[adjVar] = overrideParamB[adjVar]
                print("from overide",adjVar,overrideParamB[adjVar])
            elif adjVar in extParam.keys():
                returnRef[adjVar] = float(extParam[adjVar])
                print("from overide",adjVar,returnRef[adjVar])
            else:
                continue
            overrideParamC = []
            for adjustment in adjRange:
                overrideParamC.append(overrideParamB.copy())
                overrideParamC[-1][adjVar] = returnRef[adjVar]*adjustment
            overrideParamC = pd.DataFrame(overrideParamC)
            if method is None:
                myMethod = {'method':'hooke_jeeves'}
            else:
                myMethod = method
            estVar2 = [var for var in estimatedVar if var!=adjVar]
            returnPaths[adjVar] = os.path.join(self.run_dir,
                     'subRunPEout_'+adjVar+'.p')
            subRunPEspec = {"expDataFP":expDataFP,
                            "RPEDict":{"copyNum":depth,
                                       "estimateIC":estimateIC,
                                       "estimatedVar":estVar2,
                                       "prefix":prefix,
                                       "rocket":rocket,
                                       "upperParamBound":upperParamBound,
                                       "lowerParamBound":lowerParamBound,
                                       "method":myMethod,
                                       "overrideParam":overrideParamC,
                                       "indepToAdd":indepToAdd,
                                       "randStartVal":False,
                                       "throttle":throttle},
                            "initDict":{"antString":self.antString,
                                        "run_dir":os.path.join(
                                                self.run_dir,
                                                "pl_"+adjVar)},
                            "returnPath":returnPaths[adjVar]}
            sendPaths[adjVar] = os.path.join(self.run_dir,
                     'subRunPEspec_'+adjVar+'.p')
            file = open(sendPaths[adjVar], 'wb')
            pickle.dump(subRunPEspec,file)
            file.close()
            scriptPaths[adjVar] = os.path.join(self.run_dir, "pl_"+adjVar)
            if not os.path.isdir(scriptPaths[adjVar]):
                os.makedirs(scriptPaths[adjVar])
            scriptPaths[adjVar] = os.path.join(scriptPaths[adjVar],
                                               jobName+".sh")
            shellString = ("#!/bin/bash\n "+
                           "python -u " + __file__ + 
                           " runParamiterEstimation " +
                           sendPaths[adjVar])
            f = open(scriptPaths[adjVar], 'w')
            f.write(shellString)
            f.close()
        sys.stdout.flush()
        if rocket:
            for adjVar in estimatedVar:
                if adjVar in scriptPaths.keys():
                    os.system("sbatch " + scriptPaths[adjVar])
            while self.checkRuningProcesses(job=jobName+".sh")>0:
                time.sleep(60)
        else:
            for adjVar, myPath in sendPaths.items():
                file = open(myPath, 'rb')
                subRunPEspec = pickle.load(file)
                file.close()
                myModel = modelRunner(**subRunPEspec["initDict"])
                myParam = myModel.runParamiterEstimation(
                        subRunPEspec["expDataFP"],
                        **subRunPEspec["RPEDict"])
                myParam = myParam[next(iter(myParam))]
                file = open(subRunPEspec["returnPath"], 'wb')
                pickle.dump(myParam,file)
                file.close()
        for adjVar in estimatedVar:
            if adjVar in returnPaths.keys():
                try:
                    file = open(returnPaths[adjVar], 'rb')
                    returnData[adjVar] = pickle.load(file)
                    file.close()
                except:
                    pass
        return returnData, returnRef, {}
        
if __name__ == "__main__" and len(sys.argv[1:])>0:
    cmdLineArg = sys.argv[1:]
    if cmdLineArg[0]=="runSteadyStateFinder" and len(cmdLineArg)>=2:
        file = open(cmdLineArg[1],'rb')
        myDict = pickle.load(file)
        file.close()
        myModel = modelRunner(myDict["antimony_string"],
                              myDict["run_dir"])
        outputs = myModel.runSteadyStateFinder(params=myDict["param"])
        file = open(cmdLineArg[1],'wb')
        pickle.dump(outputs,file)
        file.close()
    elif cmdLineArg[0]=="runParamiterEstimation" and len(cmdLineArg)>=2:
        file = open(cmdLineArg[1], 'rb')
        subRunPEspec = pickle.load(file)
        file.close()
        time.sleep(5*60)
        myModel = modelRunner(**subRunPEspec["initDict"])
        myParam = myModel.runParamiterEstimation(
                subRunPEspec["expDataFP"], **subRunPEspec["RPEDict"])
        myParam = myParam[next(iter(myParam))]
        file = open(subRunPEspec["returnPath"], 'wb')
        pickle.dump(myParam,file)
        file.close()
        