this file is auto generated from the comandRecord.p and elements.p files.

A parameter estimation using code to parameterise a single model to two cancer cell cell-lines (MCF7 and ZR75) in such a way that the same kinetic parameters fit both lines but with different initial conditions. The model is based on a small  network which includes a minimal number of elements omitting many unmeasured intermediaries in the network.

filenames and identifyers:

redPaper : folder name used to identify output data and figures

zr75 : identifier used to distinguish output sets in the same folder and allow parallel execution of different parameter estimations

parameters_zr75.csv : filename of output parameter estimation

PETimeCourses_zr75.csv : time course simulations for the first few 'best' parameter sets in the parameter estimation

DBTimeCourses_zr75.csv : time course simulations for the first few 'best' parameter sets in the parameter estimation done using the same model used for parameter estimation. used for debugging

zr75_data.xlsx : file name of experimental data source

elements.p : file containing data about the named elements composing the base model and derived models used for parameter estimation

antStringRedF.txt : filename for base model we are parameterising

antStringRedM.txt : filename for model with pre stabilisation period

antStringRedMP.txt : filename for model with multiple parallel submodes

DTCon.csv : file used to link names in experimental data to names in model

dataRelations.csv : file expressing relationships between variables in experimental data. Which measurement is a contributor to which

totalsRelations.csv : file expressing which names in experimental data relate to total protein levels.

parameter estimation settings:

200 : number of parameter sets in estimated output

200.0 : time allowed for model to stabilise

100.0 : upper bound for estimated parameters

0.01 : lower bound for estimated parameters

the method and parameters used in the first phase of parameter estimations
    method : particle_swarm
    swarm_size : 100
    iteration_limit : 5000
    weight_method : standard_deviation

the method and parameters used in the second phase of parameter estimations
    method : hooke_jeeves
    weight_method : standard_deviation

elments (parameters, variables etc) in the models:

model elements in the base model
    metabolites : Akt, AktpT308, FourEBP1, FourEBP1pT37_46, IRS1, IRS1pS636_639, PRAS40, PRAS40pT246, S6K, S6KpT389, TSC2, TSC2pT1462
    kineticParamsGlobal : AA, Insulin, kAktDephos, kAktPhos_kcat, kAktPhos_kcatB, kFourEBP1Dephos, kFourEBP1Phos, kFourEBP1PhosB, kIRS1Inact, kIRS1Phos, kPras40Dephos, kPras40Phos, kS6KDephos, kS6KPhos, kS6KPhosB, kTORPRAS, kTORTSC, kTSC2Dephos, kTSC2Phos
    compartments : Cell
    assignments : AktpT308_wb, IRS1_sum, IRS1pS636_639_wb, S6KpT389_wb, FourEBP1_sum, S6K_wb, S6K_sum, IRS1_wb, TSC2pT1462_wb, FourEBP1pT37_46_wb, PRAS40pT246_wb, PRAS40_sum, TSC2_wb, PRAS40_wb, Akt_wb, TSC2_sum, FourEBP1_wb, Akt_sum

model elements in the model with stabilisation period
    metabolites : Akt, AktpT308, FourEBP1, FourEBP1pT37_46, IRS1, IRS1pS636_639, PRAS40, PRAS40pT246, PT, S6K, S6KpT389, TSC2, TSC2pT1462
    kineticParamsGlobal : AA, Akt_T, Akt_wb_DVT0, Akt_wb_DVT1, Akt_wb_DVT2, Akt_wb_DVT3, Akt_wb_DVT4, Akt_wb_DVT5, AktpT308_wb_DVT0, AktpT308_wb_DVT1, AktpT308_wb_DVT2, AktpT308_wb_DVT3, AktpT308_wb_DVT4, AktpT308_wb_DVT5, FourEBP1_T, FourEBP1_wb_DVT0, FourEBP1_wb_DVT1, FourEBP1_wb_DVT2, FourEBP1_wb_DVT3, FourEBP1_wb_DVT4, FourEBP1_wb_DVT5, FourEBP1pT37_46_wb_DVT0, FourEBP1pT37_46_wb_DVT1, FourEBP1pT37_46_wb_DVT2, FourEBP1pT37_46_wb_DVT3, FourEBP1pT37_46_wb_DVT4, FourEBP1pT37_46_wb_DVT5, IRS1_T, IRS1_wb_DVT0, IRS1_wb_DVT1, IRS1_wb_DVT2, IRS1_wb_DVT3, IRS1_wb_DVT4, IRS1_wb_DVT5, IRS1pS636_639_wb_DVT0, IRS1pS636_639_wb_DVT1, IRS1pS636_639_wb_DVT2, IRS1pS636_639_wb_DVT3, IRS1pS636_639_wb_DVT4, IRS1pS636_639_wb_DVT5, Insulin, PRAS40_T, PRAS40_wb_DVT0, PRAS40_wb_DVT1, PRAS40_wb_DVT2, PRAS40_wb_DVT3, PRAS40_wb_DVT4, PRAS40_wb_DVT5, PRAS40pT246_wb_DVT0, PRAS40pT246_wb_DVT1, PRAS40pT246_wb_DVT2, PRAS40pT246_wb_DVT3, PRAS40pT246_wb_DVT4, PRAS40pT246_wb_DVT5, PTSpeed, S6K_T, S6K_wb_DVT0, S6K_wb_DVT1, S6K_wb_DVT2, S6K_wb_DVT3, S6K_wb_DVT4, S6K_wb_DVT5, S6KpT389_wb_DVT0, S6KpT389_wb_DVT1, S6KpT389_wb_DVT2, S6KpT389_wb_DVT3, S6KpT389_wb_DVT4, S6KpT389_wb_DVT5, TSC2_T, TSC2_wb_DVT0, TSC2_wb_DVT1, TSC2_wb_DVT2, TSC2_wb_DVT3, TSC2_wb_DVT4, TSC2_wb_DVT5, TSC2pT1462_wb_DVT0, TSC2pT1462_wb_DVT1, TSC2pT1462_wb_DVT2, TSC2pT1462_wb_DVT3, TSC2pT1462_wb_DVT4, TSC2pT1462_wb_DVT5, kAktDephos, kAktPhos_kcat, kAktPhos_kcatB, kFourEBP1Dephos, kFourEBP1Phos, kFourEBP1PhosB, kIRS1Inact, kIRS1Phos, kPras40Dephos, kPras40Phos, kS6KDephos, kS6KPhos, kS6KPhosB, kTORPRAS, kTORTSC, kTSC2Dephos, kTSC2Phos, mySubsA0, mySubsA1, mySubsA2, mySubsB0, mySubsB1, mySubsB2, mySubsB3, mySubsB4, mySubsC0, mySubsC1, mySubsC2, mySubsC3, mySubsC4, mySubsD0, mySubsD1, mySubsD2, mySubsD3, mySubsD4, mySubsE0, mySubsF0
    compartments : Cell, default_compartment
    assignments : AktpT308_wb, IRS1_sum, IRS1pS636_639_wb, S6KpT389_wb, FourEBP1_sum, S6K_wb, S6K_sum, IRS1_wb, TSC2pT1462_wb, PRAS40pT246_wb, FourEBP1pT37_46_wb, PRAS40_sum, TSC2_wb, PRAS40_wb, Akt_wb, TSC2_sum, FourEBP1_wb, Akt_sum

a more compleat human readable json of the data is saved in metaData.json