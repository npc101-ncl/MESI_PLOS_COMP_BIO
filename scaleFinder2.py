#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 11:12:39 2019

@author: npc101
"""

import site, os, re
import pandas as pd
import time
from pycotools3 import model, tasks, viz

site.addsitedir('/home/npc101/anaconda3/envs/python3p6/lib/python3.6/site-packages/pycotools3')
copasiBinPath = "/home/npc101/COPASI-4.27.217-Linux-64bit/bin"

#site.addsitedir('/Users/peter/opt/anaconda3/envs/ver3p6/lib/python3.6/site-packages/pycotools3')
#copasiBinPath = "/Applications/copasi"

if not re.search(copasiBinPath, os.environ["PATH"]):
    os.environ["PATH"] += os.pathsep + copasiBinPath


working_directory = os.path.dirname(os.path.abspath(__file__))
copasi_filename = os.path.join(working_directory, 'POCModel.cps')
data_filename = os.path.join(working_directory, 't47d_data.xlsx')
output_filename = os.path.join(working_directory, "scaleSearchResults.csv")

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
        var _IRS1                in Cell;
        var _IRS1a               in Cell;
        var _IRS1pS636_639       in Cell;
        var _Akt                 in Cell;
        var _AktpT308            in Cell;
        var _TSC2                in Cell;
        var _TSC2pT1462          in Cell;
        var _PRAS40              in Cell;
        var _PRAS40pT246         in Cell;
        var _S6K                 in Cell;
        var _S6KpT389            in Cell;
        var _FourEBP1            in Cell;
        var _FourEBP1pT37_46     in Cell;
        var _PI3K                in Cell;
        var _pPI3K               in Cell;
        var _pmTORC1            in Cell;
        var _mTORC1cyt          in Cell;
        var _mTORC1lys          in Cell;
        const Insulin          in Cell;
        const AA               in Cell;

        // global variables
        Insulin = 1;
        AA = 1

        _IRS1            = 10;
        _IRS1a           = 0;
        _IRS1pS636_639   = 0;
        _PI3K            = 10;
        _pPI3K           = 0;
        _AktpT308        = 0;
        _Akt             = 10;
        _TSC2pT1462      = 10;
        _TSC2            = 0;
        _pmTORC1         = 0;
        _mTORC1cyt       = 10;
        _mTORC1lys       = 0;
        _PRAS40          = 10;
        _PRAS40pT246     = 0;
        _FourEBP1pT37_46 = 0;
        _FourEBP1        = 10;
        _S6KpT389        = 0;
        _S6K             = 10;
        
        _Akt_scale             = 1;
        _Akt_const             = 0;
        _AktpT308_scale        = 1;
        _AktpT308_const        = 0;
        _PRAS40_scale          = 1;
        _PRAS40_const          = 0;
        _PRAS40pT246_scale     = 1;
        _PRAS40pT246_const     = 0;
        _S6K_scale             = 1;
        _S6K_const             = 0;
        _S6KpT389_scale        = 1;
        _S6KpT389_const        = 0;
        _TSC2_scale            = 1;
        _TSC2_const            = 0;
        _TSC2pT1462_scale      = 1;
        _TSC2pT1462_const      = 0;
        _IRS1_scale            = 1;
        _IRS1_const            = 0;
        _IRS1pS636_639_scale   = 1;
        _IRS1pS636_639_const   = 0;
        _FourEBP1_scale        = 1;
        _FourEBP1_const        = 0;
        _FourEBP1pT37_46_scale = 1;
        _FourEBP1pT37_46_const = 0;
        
        IRS1_sum       := _IRS1+_IRS1a+_IRS1pS636_639
        PI3K_sum       := _PI3K+_pPI3K
        Akt_sum        := _Akt+_AktpT308
        TSC2_sum       := _TSC2+_TSC2pT1462
        PRAS40_sum     := _PRAS40+_PRAS40pT246
        FourEBP1_sum   := _FourEBP1+_FourEBP1pT37_46
        S6K_sum        := _S6K+_S6KpT389
        
        Akt_wb             := _Akt_scale*Akt_sum+_Akt_const
        AktpT308_wb        := _AktpT308_scale*AktpT308+_AktpT308_const
        PRAS40_wb          := _PRAS40_scale*PRAS40_sum+_PRAS40_const
        PRAS40pT246_wb     := _PRAS40pT246_scale*PRAS40pT246+_PRAS40pT246_const
        S6K_wb             := _S6K_scale*S6K_sum+_S6K_const
        S6KpT389_wb        := _S6KpT389_scale*S6KpT389+_S6KpT389_const
        TSC2_wb            := _TSC2_scale*TSC2_sum+_TSC2_const
        TSC2pT1462_wb      := _TSC2pT1462_scale*TSC2pT1462+_TSC2pT1462_const
        IRS1_wb            := _IRS1_scale*IRS1_sum+_IRS1_const
        IRS1pS636_639_wb   := _IRS1pS636_639_scale*IRS1pS636_639+_IRS1pS636_639_const
        FourEBP1_wb        := _FourEBP1_scale*FourEBP1_sum+_FourEBP1_const
        FourEBP1pT37_46_wb := _FourEBP1pT37_46_scale*FourEBP1pT37_46+_FourEBP1pT37_46_const

        // kinetic parameters
        _kIRS1Act                = 0.1;
        _kIRS1Inact              = 0.1;
        _kIRS1Phos               = 0.1;
        _kPI3KPhos               = 0.1;
        _kPI3KDephos             = 0.1;
        _kAktPhos_km            = 0.1;
        _kAktPhos_kcat          = 0.1;
        _kAktDephos              = 0.1;
        _kTSC2Phos               = 0.1;
        _kTSC2Dephos             = 0.1;
        _kmTORC1cytToLys         = 0.1;
        _kmTORC1LysToCyt         = 0.1;
        _kmTORC1Phos             = 0.1;
        _kmTORC1Dephos           = 0.1;
        _kPras40Phos             = 0.1;
        _kPras40Dephos           = 0.1;
        _kFourEBP1Phos           = 0.1;
        _kFourEBP1Dephos         = 0.1;
        _kS6KPhos                = 0.1;
        _kS6KDephos              = 0.1;



        // reactions // MMWithKcat(km, kcat, S, E)
        // function CompetitiveInhibitionWithKcat(km, ki, kcat, E, I, S)
        // function MM(km, Vmax, S)
        // function NonCompetitiveInhibitionWithKcat(km, ki, kcat, E, n, I, S)
        R1f     : IRS1 => IRS1a                         ; Cell * _kIRS1Act*_IRS1*Insulin;
        R2f     : IRS1a => IRS1pS636_639                ; Cell * _kIRS1Phos*_IRS1a*_S6KpT389
        R2b     : IRS1pS636_639 => IRS1                 ; Cell * _kIRS1Inact*_IRS1pS636_639
        R3f     : PI3K => pPI3K                         ; Cell * _kPI3KPhos*_PI3K*_IRS1a
        R3b     : pPI3K => PI3K                         ; Cell * _kPI3KDephos*_pPI3K
        R4f     : Akt => AktpT308                       ; Cell * MMWithKcat(_kAktPhos_km, _kAktPhos_kcat, _Akt, _pPI3K)
        R4b     : AktpT308 => Akt                       ; Cell * _kAktDephos*_pPI3K*_AktpT308
        R5f     : TSC2 => TSC2pT1462                    ; Cell * _kTSC2Phos*_TSC2*_AktpT308
        R5b     : TSC2pT1462 => TSC2                    ; Cell * _kTSC2Dephos*_TSC2pT1462
        R6f     : mTORC1cyt => mTORC1lys                ; Cell * _kmTORC1cytToLys*_mTORC1cyt*AA
        R6b     : mTORC1lys => mTORC1cyt                ; Cell * _kmTORC1LysToCyt*_mTORC1lys
        R7f     : mTORC1lys => pmTORC1                  ; Cell * _kmTORC1Phos*_mTORC1lys*_TSC2
        R7b     : pmTORC1 => mTORC1lys                  ; Cell * _kmTORC1Dephos*_pmTORC1
        R8f     : PRAS40 => PRAS40pT246                 ; Cell * _kPras40Phos*_PRAS40*_AktpT308
        R8b     : PRAS40pT246 => PRAS40                 ; Cell * _kPras40Dephos*_PRAS40pT246
        R9f     : FourEBP1 => FourEBP1pT37_46           ; Cell * _kFourEBP1Phos*_FourEBP1*_pmTORC1
        R9b     : FourEBP1pT37_46 => FourEBP1           ; Cell * _kFourEBP1Dephos*_FourEBP1pT37_46
        R10f    : S6K => S6KpT389                       ; Cell * _kS6KPhos*_S6K*_pmTORC1
        R10b    : S6KpT389 => S6K                       ; Cell * _kS6KDephos*_S6KpT389

    end

    """

if __name__ == "__main__":
    #comm = MPI.COMM_WORLD
    #taskCount = comm.Get_size()
    #taskNumber = comm.Get_rank()
    
    #takes a string turns it into a string buffer, turns that into a fake file
    #then reads it as if it were a csv into a pandas data frame
    df = pd.read_excel (data_filename,sheet_name='Sheet2',skiprows=51,
                        index_col=0,usecols=21,nrows=12)
    df = df[df.index.str.startswith('MCF-7')]
    
    df.index = [float(row[1])*60 for row in df.index.str.split(' ')]
    df.index.name = "Time"
    df.rename(columns={'IRS1pS636/639':'IRS1pS636_639_wb',
                       '4E-BP1':'FourEBP1_wb',
                       '4E-BP1pT37/46':'FourEBP1pT37_46_wb',
                       'Akt':'Akt_wb',
                       'PRAS40':'PRAS40_wb',
                       'S6K':'S6K_wb',
                       'TSC2':'TSC2_wb',
                       'IRS1':'IRS1_wb',
                       'AktpT308':'AktpT308_wb',
                       'PRAS40pT246':'PRAS40pT246_wb',
                       'S6KpT389':'S6KpT389_wb',
                       'TSC2pT1462':'TSC2pT1462_wb',
                       'FourEBP1pT37_46':'FourEBP1pT37_46_wb'}
    , inplace=True)
    df=df.drop(['AktpS473', 'PRAS40pS183', 'S6KpT229', 'GAPDH', 'ERK',
             'ERK-pT202/Y202','p38', 'p38-pT180/Y182', 'ER alpha'], axis=1)
    
    with model.BuildAntimony(copasi_filename) as loader:
        my_model = loader.load(antimony_string)
    fname = os.path.join(working_directory, 'POC_experimental_data.csv')
    df.to_csv(fname)
    my_copy_number=1
    with tasks.ParameterEstimation.Context(my_model, fname,
                                           context='s',
                                           parameters='gm') as context:
        context.set('method', 'particle_swarm') # use just g if you don't 
        context.set('population_size', 10) # want initial consentrations
        context.set('number_of_generations', 350)
        context.set('iteration_limit',6000)
        context.set('swarm_size',100)
        context.set('copy_number', my_copy_number)
        context.set('pe_number', 1)
        context.set('run_mode', True) #'parallel' ... 'slurm'
        #context.set('run_mode', 'slurm')
        context.set('separator', ',')
        context.set('prefix', '_')
        #context.set('working_directory', 'wd', $TMPDIR)
        config = context.get_config()
    pe = tasks.ParameterEstimation(config)
    #pe.models['POCModel'].model.open()
    number_completed=0
    number_last_check=0;
    frameString='POCModel'
    while number_completed<my_copy_number:
        time.sleep(10)
        data=viz.Parse(pe)[frameString]
        if isinstance(data, pd.DataFrame):
            number_completed = data.shape[0]
            print("Done",number_completed,"of",my_copy_number)
            if number_completed>=number_last_check+10:
                data.to_csv(output_filename)
                number_last_check=number_completed
    data.to_csv(output_filename)