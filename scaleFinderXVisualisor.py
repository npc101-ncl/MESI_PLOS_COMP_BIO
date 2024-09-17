#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:43:36 2020

@author: peter
"""

import PySimpleGUI as sg
import pickle, os
from python.visualisationTools import *
from python.analysisTools import *
from scaleFinderXPrimer import importData, makeInterpolatedTable, applyScale
import numpy as np

working_directory = os.path.dirname(os.path.abspath(__file__))

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

print("does begin")
sg.theme('DarkAmber')   # Add a touch of color
# All the stuff inside your window.
layout = [  [sg.Text('TimeCourse path:'), sg.Input(key='tcPath'),
             sg.FileBrowse()],
            [sg.Text('Run Switches path:'), sg.Input(key='rsPath'),
             sg.FileBrowse()],
            [sg.Text('Paramiter path:'), sg.Input(key='pPath'),
             sg.FileBrowse()],
            [sg.Button('Load')],
            [sg.Combo([],key='addVar'), sg.Button('Add'),
             sg.Combo([],key='remVar'), sg.Button('Remove')],
            [sg.Text('selected: '), sg.Text(size=(12,1), key='varText')],
            [sg.Text('data: '), sg.Combo([],key='dataVar'),
             sg.Text('scaling: '), sg.Combo([],key='runNameVar'),
             sg.Button('Select')],
            [sg.Text('Figure directory:'), sg.Input(key='fPath'),
             sg.FolderBrowse()],
            [sg.Text('Max Time Courses:'), sg.Spin([i for i in range(1,10)],
                     initial_value=10,key="maxTC")],
            [sg.Button('Save'), sg.Button('Cancel')] ]

# Create the Window
window = sg.Window('Create and save time course plots', layout)

BestClusterSize = None
df = None
runNames=[]
dataNames=[]

# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event in (None, 'Cancel'):   # if user closes window or clicks cancel
        break
    elif event == 'Load':
        # change the "output" element to be the value of "input" element
        if os.path.isfile(values['tcPath']):
            f = open(values['tcPath'], "rb" )
            timeCourse = pickle.load(f)
            f.close()
            varPot = [myStr for myStr in list(timeCourse[0].columns) if
                      myStr!="Time"]
            varSel = []
            window['addVar'].update(values=varPot, set_to_index=0)
            window['addVar'].set_size(size=(max([len(myStr) for myStr in
                  varPot]),None))
            window['remVar'].update(values=[])
            window['remVar'].set_size(size=(max([len(myStr) for myStr in
                  varPot]),None))
            window['varText'].update(", ".join(varSel))
        if os.path.isfile(values['rsPath']):
            f = open(values['rsPath'], "rb" )
            RS = pickle.load(f)
            f.close()
            dataNames = list(RS["data_filenames"].keys())
            runNames = [name[6:] for name in RS.keys() if
                        name.startswith("scale-")]
            window['dataVar'].update(values=dataNames, set_to_index=0)
            window['dataVar'].set_size(size=(max([len(myStr) for myStr in
                  dataNames]),None))
            window['runNameVar'].update(values=runNames, set_to_index=0)
            window['runNameVar'].set_size(size=(max([len(myStr) for myStr in
                  runNames]),None))
        if os.path.isfile(values['pPath']):
            f = open(values['pPath'], "rb" )
            params = pickle.load(f)
            f.close()
            params = {key:val[val["RSS"]!=np.inf] for key, val
                      in params.items()}
            BestClusterSize = RSSClusterEstimation(GFID(params))[0]["size"]
    elif event == "Select":
        if values['dataVar'] in dataNames:
            myPath = pathToLocal(RS["data_filenames"][values['dataVar']],
                                 working_directory)
            df = importData('MCF-7', myPath, RS["columns_from_data"])
            if values['runNameVar'] in runNames:
                df = applyScale(df, RS["scale-"+values['runNameVar']])
        print(df)
    elif event == 'Add':
        varSel.append(values['addVar'])
        window['addVar'].update(values=[var for var in varPot if not var
              in varSel], set_to_index=0)
        window['remVar'].update(values=varSel)
        window['varText'].update(", ".join(varSel))
        window['varText'].set_size(size=(len(", ".join(varSel)),None))
    elif event == 'Remove':
        varSel.remove(values['remVar'])
        window['addVar'].update(values=[var for var in varPot if not var
              in varSel], set_to_index=0)
        window['remVar'].update(values=varSel)
        window['varText'].update(", ".join(varSel))
        window['varText'].set_size(size=(len(", ".join(varSel)),None))
    elif event == 'Save':
        sns.set(context='paper')
        timeCourseVis = timeCourseVisualiser(timeCourse)
        if BestClusterSize is not None:
            showList = list(range(min(BestClusterSize,
                                      int(values["maxTC"]))))
        else:
            showList = list(range(int(values["maxTC"])))
        figPath=None
        if os.path.isdir(values['fPath']):
            figPath = os.path.join(values['fPath'],
                                   os.path.split(values['tcPath'])[1])
        if (df is not None) and (figPath is not None):
            timeCourseVis.multiPlot(indexSelect=showList, varSelect=varSel,
                                    compLines=df, save=figPath)
        elif figPath is not None:
            timeCourseVis.multiPlot(indexSelect=showList, varSelect=varSel,
                                    save=figPath)
        elif df is not None:
            timeCourseVis.multiPlot(indexSelect=showList, varSelect=varSel,
                                    compLines=df)
        else:
            timeCourseVis.multiPlot(indexSelect=showList, varSelect=varSel)
        
        
        

window.close()