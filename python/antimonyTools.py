#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 11:08:41 2022

@author: peter
"""

import sympy
from itertools import combinations
from python.antimonyTools import *
import re

class antToSympyReader:
    def __init__(self,inAntStr):
        if isinstance(inAntStr,str):
            antStr = inAntStr
        else:
            return None
        antStr = antStr.splitlines()
        antStr = [i.split('#',1)[0] for i in antStr]
        antStr = [i.split('//',1)[0] for i in antStr]
        antStr = [i.strip() for i in antStr]
        antStr = [i for i in antStr if i!='']
        self.antStr = antStr
        self.functions = self.getFunctions(antStr)
        reactions = {i:self.getReaction(r) for i,r in enumerate(antStr)}
        self.reactions = {i:r for i,r in reactions.items() if r is not None}
        var = [[j["variable"] for j in i["ins"]+i["outs"] if not j["fixed"]] 
               for _,i in self.reactions.items()]
        self.var = list(set([i for l in var for i in l]))
        self.ODEs = {var:self.odeFromReactions(var) for var in self.var}
        
    def getFunctions(self,myLines):
        funcs = []
        for i,line in enumerate(myLines[:-2]):
            if not line.startswith("function "):
                continue
            t = line.split(" ",1)[1].strip()
            if not t.endswith(")"):
                continue
            t = t[:-1]
            t = t.split("(")
            if len(t)!=2:
                continue
            p = t[1]
            t = t[0].strip()
            p = [j.strip() for j in p.split(",")]
            if " " in t:
                continue
            if myLines[i+2]!="end":
                continue
            r = [i.strip() for i in myLines[i+1].split(";") 
                 if len(i.strip())>0]
            if len(r)==1:
                r = r[0]
            else:
                continue
            funcs.append({"name":t,"params":p,"exp":r})
        return funcs
    
    def insFunction(self,exp,funcs,tBEnd="@!@"):
        eStr = exp
        notSpent = True
        tokens = {}
        tcount = 0
        #tokenise
        while notSpent:
            for i, c1 in enumerate(eStr):
                if c1 == "(":
                    t = eStr[i+1:].split(")",1)
                    if len(t)<2:
                        return None
                    if "(" in t[0]:
                        continue
                    tokens[tBEnd+str(tcount)+tBEnd] = t[0]
                    eStr = (eStr[:i]+tBEnd+str(tcount)+tBEnd+
                            eStr[i+len(t[0])+2:])
                    tcount += 1
                    break
            else:
                notSpent = False
        tokens = {k:v.split(",") for k,v in tokens.items()}
        notSpent = True
        while notSpent:
            for k,v in tokens.items():
                t = eStr.split(k,1)
                if len(t)<2:
                    continue
                else:
                    t[0] = t[0].rstrip()
                    for func in funcs:
                        if t[0].endswith(func["name"]):
                            if len(v)!=len(func["params"]):
                                return None
                            t2 = t[0][:-len(func["name"])]
                            u =  func["exp"]
                            for i,j in zip(v,func["params"]):
                                u = re.sub(r'\b'+j+r'\b', "("+i+")", u)
                                #u = u.replace(j,"("+i+")")
                            eStr = t2+"("+u+")"+t[1]
                            break
                    else:
                        if len(v)!=1:
                            return None
                        else:
                            eStr = eStr.replace(k,"("+v[0]+")")
                    break
            else:
                notSpent = False
        return eStr

    def getReaction(self,myStr):
        reversable = False
        i = myStr.split(';')
        if len(i)>3:
            return None
        elif len(i)>=2:
            expresion = i[1].strip()
            reversable = True
        i = i[0]
        i = i.split(':')
        if len(i)==2:
            name = i[0].strip()
            i = i[1]
        elif len(i)==1:
            i=i[0]
            name = None
        else:
            return None
        i = i.split("=>")
        if len(i)>2:
            return None
        elif len(i)==2:
            reversable = True
        else:
            i = i[0].split("->")
            if len(i)!=2:
                return None
        ins = i[0].split("+")
        outs = i[1].split("+")
        ins = [[j for j in i.split() if j!=''] for i in ins]
        outs = [[j for j in i.split() if j!=''] for i in outs]
        ins = [i for i in ins if len(i)!=0]
        outs = [i for i in outs if len(i)!=0]
        for i in range(len(ins)):
            if len(ins[i])==1:
                ins[i]={"multiplicity":"1",
                           "variable":ins[i][0]}
            elif len(ins[i])==2:
                ins[i]={"multiplicity":ins[i][0],
                           "variable":ins[i][1]}
            else:
                return None
            try:
                ins[i]["multiplicity"]=int(ins[i]["multiplicity"])
            except:
                return None
            if ins[i]["variable"][0]=='$':
                ins[i]["variable"] = ins[i]["variable"][1:]
                ins[i]["fixed"] = True
            else:
                ins[i]["fixed"] = False
        for i in range(len(outs)):
            if len(outs[i])==1:
                outs[i]={"multiplicity":"1",
                           "variable":outs[i][0]}
            elif len(outs[i])==2:
                outs[i]={"multiplicity":outs[i][0],
                           "variable":outs[i][1]}
            else:
                return None
            try:
                outs[i]["multiplicity"]=int(outs[i]["multiplicity"])
            except:
                return None
            if outs[i]["variable"][0]=='$':
                outs[i]["variable"] = outs[i]["variable"][1:]
                outs[i]["fixed"] = True
            else:
                outs[i]["fixed"] = False   
        expresion = expresion.replace("^", "**")
        expresion = self.insFunction(expresion,self.functions)
        try:
            expresion = sympy.parsing.sympy_parser.parse_expr(expresion)
        except:
            return None
        return {"name":name, "reversable":reversable, "ins":ins, "outs":outs,
                "expresion":expresion}
        
    def odeFromReactions(self,var):
        ins = [sum([j["multiplicity"] for j in i["ins"]
                    if j["variable"]==var and not j["fixed"]])*i["expresion"]
               for _,i in self.reactions.items()]
        outs = [sum([j["multiplicity"] for j in i["outs"]
                    if j["variable"]==var and not j["fixed"]])*i["expresion"]
               for _,i in self.reactions.items()]
        return sum(outs)-sum(ins)
    
class stabSolver:
    def __init__(self,ODEs,additionalCons = {}):
        """
            ODEs = ODEs to consider in solutions
            additionalCons = aditional considerations / substitutions to use
        """
        self.ODEs = ODEs
        self.extSymb = set()
        for _,i in self.ODEs.items():
            self.extSymb = self.extSymb.union(i.free_symbols)
        self.runningSubs = {k:sympy.parsing.sympy_parser.parse_expr(v) 
                            for k,v in additionalCons.items()}
        self.runningSubs = {sympy.parsing.sympy_parser.parse_expr(k):v 
                            for k,v in self.runningSubs.items()}
        self.newSymb = set()
        for _,i in self.runningSubs.items():
            self.newSymb = self.newSymb.union(i.free_symbols)
        self.newSymb = self.newSymb - self.extSymb
        self.subsTagPrefix = "mySubs"
        self.subItir = 0
        self.refSubs = {}
        
    def getPrefixTag(self):
        t = str(self.subItir)
        for k, v in {"0":"A", "1":"B", "2":"C", "3":"D", "4":"E", "5":"F",
                     "6":"G", "7":"H", "8":"I", "9":"J"}.items():
            t = t.replace(k,v)
        self.subItir += 1
        return self.subsTagPrefix+t
        
    def exprCleaner(self,expr,retainers,symbols):
        """
            expr = sympy expresion
            retainers = tupal of sympy symbols
            symbols = sympy.utilities.iterables.numbered_symbols
        """
        retSubs = {}
        if len(expr.args)==0:
            return expr, retSubs
        elif not any([(j in expr.free_symbols) for j in retainers]):
            newSym = next(symbols)
            retSubs[newSym]=expr
            return newSym, retSubs
        elif expr.func in [sympy.core.add.Add, sympy.core.mul.Mul]:
            argsL = list(expr.args)
            argsS = [i for i in argsL if not any([(j in i.free_symbols) 
                                                  for j in retainers])]
            newSym = next(symbols)
            retSubs[newSym]=expr.func(*argsS)
            argsL = [i for i in argsL if any([(j in i.free_symbols) 
                                              for j in retainers])]
            for i in range(len(argsL)):
                argsL[i], tempSubs = self.exprCleaner(argsL[i], retainers,
                                                      symbols)
                retSubs.update(tempSubs)
            argsL = [newSym]+argsL
            exprR = expr.func(*argsL)
            return exprR, retSubs
        else:
            argsL = list(expr.args)
            for i in range(len(argsL)):
                argsL[i], tempSubs = self.exprCleaner(argsL[i], retainers,
                                                      symbols)
                retSubs.update(tempSubs)
            exprR = expr.func(*argsL)
            return exprR, retSubs
        
    def subSolveSim(self,temp,targets,retainInExp,subsTag = None,debug=False):
        """
            temp = list of strings indicating expresions to solve
            retainInExp = list of strings for expresions to retain in solution
            targets = solutions we want as list of str
            self.runningSubs = dictionary of substitutions to put into system 
                               before solving
            subsTag = label to use for naming dummy variables
        """
        ignores = tuple(sympy.parsing.sympy_parser.parse_expr(i) for i
                        in retainInExp)
        if subsTag is None:
            subsTag = self.getPrefixTag()
        temp2 = [self.ODEs[i] for i in temp]
        for i in range(len(temp2)):
            for k, v in self.runningSubs.items():
                temp2[i] = temp2[i].subs(k,v)
        temp2 = sympy.solve(temp2,targets)
        keyStore = [i for i,_ in temp2.items()]
        temp2 = [sympy.collect(sympy.cancel(i),list(ignores)) for _,i 
                 in temp2.items()]
        symGen = sympy.utilities.iterables.numbered_symbols(subsTag)
        #outSubs, temp = sympy.cse(temp,ignore = ignores, symbols = symGen)
        outSubs = {}
        for i in range(len(temp2)):
            temp2[i], temp3 = self.exprCleaner(temp2[i], ignores, symGen)
            outSubs.update(temp3)
        temp2 = {k:v for k,v in zip(keyStore,temp2)}
        if debug:
            print(temp2)
        self.runningSubs.update(temp2)
        self.refSubs.update(outSubs)
        return temp2
    
    def decomposeSets(self,listOfSets,maxOverlap):
        residues = listOfSets
        pullOuts = []
        for i in range(maxOverlap,1,-1):
            temp = residues
            for j in combinations(residues, i):
                t = set.intersection(*j)
                if len(t)>0:
                    pullOuts.append(t)
                    for k in range(len(temp)):
                        temp[k] = temp[k]-t
            residues = temp
        return pullOuts+residues
    
    def decomposeList(self,myList,dictOfSets):
        mySet = set(myList)
        myExp = [k for k,v in dictOfSets.items() if v.issubset(mySet)]
        return sympy.core.add.Add(*myExp)
    
    def asignICGuesses(self,initCondDicts):
        setPeices = self.decomposeSets([set(v) for _,v 
                                        in initCondDicts.items()],
                                       len(initCondDicts))
        symGen = sympy.utilities.iterables.numbered_symbols("myPeices")
        setPeices = {k:v for k,v in zip(symGen,setPeices)}
        linToSolve = {sympy.parsing.sympy_parser.parse_expr(k):
                      self.decomposeList(v, setPeices) 
                      for k,v in initCondDicts.items()}
        addNew = set(linToSolve.keys())
        testDict = linToSolve.copy()
        reductionOrder = {i:sum([i in v.args for _,v in linToSolve.items()])
                          for i in setPeices.keys()}
        reductionOrder = [[k for k,v in reductionOrder.items() if v==i] 
                          for i in range(len(linToSolve),0,-1)]
        
        assinments = {}   
        for i,insertions in zip(range(len(linToSolve),0,-1),reductionOrder):
            for k in list(linToSolve.keys()):
                if isinstance(linToSolve[k],sympy.core.symbol.Symbol):
                    assinments[linToSolve[k]] = k
                    del linToSolve[k]
            if i>1:
                denominator = (sympy.core.add.Add(*list(linToSolve.keys()))
                               )**(i-1)
                temp = {}
                for peice in insertions:
                    inclusion = {k:(peice in v.args) 
                                 for k,v in linToSolve.items()}
                    numerator = [k for k,v in inclusion.items() if v]
                    numerator = sympy.core.mul.Mul(*numerator)
                    temp[peice] = numerator/denominator
                for peice, mySub in temp.items():
                    for k in list(linToSolve.keys()):
                        if peice in linToSolve[k].args:
                            linToSolve[k-mySub] = sympy.core.add.Add(
                                    *[i for i in linToSolve[k].args 
                                      if i!=peice])
                            del linToSolve[k]
                assinments.update(temp)
        assinments = {k:sympy.simplify(v) for k,v in assinments.items()}
        for peice, value in assinments.items():
            testDict = {k:v.subs(peice,value) for k,v in testDict.items()}
        testDict = {k:sympy.simplify(v) for k,v in testDict.items()}
        for k,v in testDict.items():
            if k!=v:
                return None
        outVals = {}
        for mySym, peices in setPeices.items():
            for i in peices:
                outVals[i] = assinments[mySym]/len(peices)
        outVals = {sympy.parsing.sympy_parser.parse_expr(k):v 
                   for k,v in outVals.items()}
        self.runningSubs.update(outVals)
        self.newSymb = self.newSymb.union(addNew)
        return outVals
    
    def amendAntStr(self,antStr):
        mylines = antStr.splitlines()
        toReplace = [(lambda x: x[0] if len(x)==2 
                      and (not x[1].startswith(">")) 
                      and (not x[0].endswith(":")) else None)(i.split("=",1)) 
                     for i in mylines]
        toReplace = {i:j.strip() for i,j in enumerate(toReplace) 
                     if j is not None}
        toReplace = {i:j for i,j in toReplace.items() if not " " in j}
        myIndentation = {k:mylines[k].split(v)[0] for k,v in toReplace.items()}
        myReplacments = {str(k):str(v).replace("**","^") for k,v 
                         in self.runningSubs.items()}
        toReplace = {k:myIndentation[k]+v+" = "+myReplacments[v] + " ;" 
                     for k,v in toReplace.items() if v in myReplacments.keys()}
        print(toReplace)
        for k,v in toReplace.items():
            mylines[k] = v
        indentCount = {}
        insPoint = min([i for i in toReplace])
        for k,v in myIndentation.items():
            if v in indentCount:
                indentCount[v] += 1
            else:
                indentCount[v] = 1
        i = 0
        for k,v in indentCount.items():
            if v>i:
                myIndentation = k
                i=v
        myInsertions = {str(k):str(v).replace("**","^") for k,v 
                        in self.refSubs.items()}
        myInsertions = [myIndentation+k+" = "+v+" ;" 
                        for k,v in myInsertions.items()]
        mylines = mylines[:insPoint] + myInsertions + mylines[insPoint:]
        myInsertions = [myIndentation+str(i)+" = 1 ;" for i in self.newSymb]
        mylines = mylines[:insPoint] + myInsertions + mylines[insPoint:]
        return "\n".join(mylines)
    
def genRunInAntStr(myODEs,varToTrack,times,antStr,stimulation={},preStim={},
                   dummyVarSuffix = "_DVT",ptSpeedVar = "PTSpeed",
                   ptVar = "PT",reactionLabel = "RTX",
                   stabilityLimit = 0.0001,dummyStart = 10000,
                   stabilityModifyer={}):
    stabCond = ["abs("+str(v).replace("**","^")+")<"+str(stabilityLimit)+
                (lambda x: "*"+str(stabilityModifyer[x]) 
                           if x in stabilityModifyer else "")(k)
                for k,v in myODEs.items()]
    stabCond = " && ".join(stabCond)
    stabCond = ptSpeedVar + "<0.5 && " + stabCond
    initialEvent = "at ("+stabCond+"), t0=false: "+ptSpeedVar+"=1"
    dummySetEvents = []
    for k,v in stimulation.items():
        initialEvent += ", "+k+"="+str(v)
    if times[0]==0:
        for i,j in zip(varToTrack,[i+dummyVarSuffix+"0" for i in varToTrack]):
            initialEvent += ", "+j+"="+i
    else:
        newStr = "at ("+ptVar+">="+str(times[0])+"): "
        newStr += ", ".join([i+dummyVarSuffix+"0="+i for i in varToTrack])
        dummySetEvents.append(newStr)
    for i,t in enumerate(times[1:],1):
        newStr = "at ("+ptVar+">="+str(t)+"): "
        newStr += ", ".join([j+dummyVarSuffix+str(i)+"="+j for j in varToTrack])
        dummySetEvents.append(newStr)
    reactions = [reactionLabel+ ": ->"+ptVar+" ; "+ptSpeedVar,
                 ptVar+" = 0",
                 ptSpeedVar+" = 0"]
    dummyIC = [i+dummyVarSuffix+str(j) for i in varToTrack 
               for j in range(len(times))]
    dummyIC = [i+" = "+str(dummyStart)+" ;" for i in dummyIC]   
    myLines = antStr.splitlines()
    myBegin = [i for i,j in enumerate(myLines) if j.strip().startswith("model ")]
    if len(myBegin)==0:
        return None
    myBegin = myBegin[0]
    myEnd = [i for i,j in enumerate(myLines) if j.strip()=="end"]
    myEnd = [i for i in myEnd if i>myBegin]
    if len(myEnd)==0:
        return None
    myEnd = myEnd[0]
    indents = {}
    for i in myLines[myBegin+1:myEnd]:
        j = len(i)-len(i.lstrip())
        if i[:j] in indents:
            indents[i[:j]] += 1
        else:
            indents[i[:j]] = 1
    indents = max(indents)
    for k,v in preStim.items():
        t = {i:j.split("=") for i,j in enumerate(myLines)}
        t = {i:j[0] for i,j in t.items() if len(j)==2 
             and (not j[0].endswith(":")) and (not j[0].startswith(":")) 
             and j[0].strip()==k}
        t = {i:j[:len(j)-len(j.lstrip())] for i,j in t.items()}
        for i,j in t.items():
            myLines[i] = j+k+" = "+str(v)+" ;"
    newLines = dummyIC+reactions+[initialEvent]+dummySetEvents
    newLines = myLines[:myEnd]+[indents+i for i in newLines]+myLines[myEnd:]
    return "\n".join(newLines)