this file is auto generated from the comandRecord.p and elements.p files.

A parameter estimation using code to parameterise a single model to two cancer cell cell-lines (MCF7 and ZR75) in such a way that the same kinetic parameters fit both lines but with different initial conditions. The model is based on a medium sized network which includes interactions between mTOR RICTOR and RPTOR.

filenames and identifyers:

paper : folder name used to identify output data and figures

test : identifier used to distinguish output sets in the same folder and allow parallel execution of different parameter estimations

parameters.csv : filename of output parameter estimation

PETimeCourses.csv : time course simulations for the first few 'best' parameter sets in the parameter estimation

DBTimeCourses.csv : time course simulations for the first few 'best' parameter sets in the parameter estimation done using the same model used for parameter estimation. used for debugging

zr75_data.xlsx : file name of experimental data source

elements.p : file containing data about the named elements composing the base model and derived models used for parameter estimation

antString2.txt : filename for base model we are parameterising

antString2M.txt : filename for model with pre stabilisation period

antString2MP.txt : filename for model with multiple parallel submodes

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
    metabolites : Akt, AktpS473, AktpT308, AktpT308S473, FourEBP1, FourEBP1pT37_46, IRS1, IRS1a, IRS1pS636_639, PI3K, PRAS40, PRAS40pT246, RICTOR, RPTOR, S6K, S6KpT389, TSC2, TSC2pT1462, mTOR, mTORC1cyt, mTORC1lys, mTORC2, pPI3K, pmTORC1, pmTORC2
    kineticParamsGlobal : AA, AAB, Insulin, InsulinB, kAktByTor2Dephos, kAktByTor2Phos, kAktDephos, kAktPhos_kcat, kFourEBP1Dephos, kFourEBP1Phos, kIRS1Act, kIRS1Inact, kIRS1Phos, kPI3KDephos, kPI3KPhos, kPras40Dephos, kPras40Phos, kRAS, kS6KDephos, kS6KPhos, kTSC2Dephos, kTSC2Phos, kTSC2PhosBoost, kmTORC1Dephos, kmTORC1LysToCyt, kmTORC1Phos, kmTORC1_Dis, kmTORC1cytToLys, kmTORC2Dephos, kmTORC2DephosByS6K, kmTORC2_Dis, kmTOR_RIC_Comb, kmTOR_RPT_Comb
    compartments : Cell
    assignments : PRAS40_sum, S6K_sum, TSC2_wb, AktpT308_wb, TSC2pT1462_wb, IRS1_wb, IRS1pS636_639_wb, AktpS473_wb, FourEBP1pT37_46_wb, FourEBP1_wb, PRAS40_wb, Akt_wb, S6KpT389_wb, PI3K_sum, FourEBP1_sum, S6K_wb, IRS1_sum, PRAS40pT246_wb, Akt_sum, TSC2_sum

model elements in the model with stabilisation period
    metabolites : Akt, AktpS473, AktpT308, AktpT308S473, FourEBP1, FourEBP1pT37_46, IRS1, IRS1a, IRS1pS636_639, PI3K, PRAS40, PRAS40pT246, PT, RICTOR, RPTOR, S6K, S6KpT389, TSC2, TSC2pT1462, mTOR, mTORC1cyt, mTORC1lys, mTORC2, pPI3K, pmTORC1, pmTORC2
    kineticParamsGlobal : AA, AAB, Akt_T, Akt_wb_DVT0, Akt_wb_DVT1, Akt_wb_DVT2, Akt_wb_DVT3, Akt_wb_DVT4, Akt_wb_DVT5, AktpS473_wb_DVT0, AktpS473_wb_DVT1, AktpS473_wb_DVT2, AktpS473_wb_DVT3, AktpS473_wb_DVT4, AktpS473_wb_DVT5, AktpT308_wb_DVT0, AktpT308_wb_DVT1, AktpT308_wb_DVT2, AktpT308_wb_DVT3, AktpT308_wb_DVT4, AktpT308_wb_DVT5, FourEBP1_T, FourEBP1_wb_DVT0, FourEBP1_wb_DVT1, FourEBP1_wb_DVT2, FourEBP1_wb_DVT3, FourEBP1_wb_DVT4, FourEBP1_wb_DVT5, FourEBP1pT37_46_wb_DVT0, FourEBP1pT37_46_wb_DVT1, FourEBP1pT37_46_wb_DVT2, FourEBP1pT37_46_wb_DVT3, FourEBP1pT37_46_wb_DVT4, FourEBP1pT37_46_wb_DVT5, IRS1_T, IRS1_wb_DVT0, IRS1_wb_DVT1, IRS1_wb_DVT2, IRS1_wb_DVT3, IRS1_wb_DVT4, IRS1_wb_DVT5, IRS1pS636_639_wb_DVT0, IRS1pS636_639_wb_DVT1, IRS1pS636_639_wb_DVT2, IRS1pS636_639_wb_DVT3, IRS1pS636_639_wb_DVT4, IRS1pS636_639_wb_DVT5, Insulin, InsulinB, PI3K_T, PRAS40_T, PRAS40_wb_DVT0, PRAS40_wb_DVT1, PRAS40_wb_DVT2, PRAS40_wb_DVT3, PRAS40_wb_DVT4, PRAS40_wb_DVT5, PRAS40pT246_wb_DVT0, PRAS40pT246_wb_DVT1, PRAS40pT246_wb_DVT2, PRAS40pT246_wb_DVT3, PRAS40pT246_wb_DVT4, PRAS40pT246_wb_DVT5, PTSpeed, RICTOR_T, RPTOR_T, S6K_T, S6K_wb_DVT0, S6K_wb_DVT1, S6K_wb_DVT2, S6K_wb_DVT3, S6K_wb_DVT4, S6K_wb_DVT5, S6KpT389_wb_DVT0, S6KpT389_wb_DVT1, S6KpT389_wb_DVT2, S6KpT389_wb_DVT3, S6KpT389_wb_DVT4, S6KpT389_wb_DVT5, TSC2_T, TSC2_wb_DVT0, TSC2_wb_DVT1, TSC2_wb_DVT2, TSC2_wb_DVT3, TSC2_wb_DVT4, TSC2_wb_DVT5, TSC2pT1462_wb_DVT0, TSC2pT1462_wb_DVT1, TSC2pT1462_wb_DVT2, TSC2pT1462_wb_DVT3, TSC2pT1462_wb_DVT4, TSC2pT1462_wb_DVT5, kAktByTor2Dephos, kAktByTor2Phos, kAktDephos, kAktPhos_kcat, kFourEBP1Dephos, kFourEBP1Phos, kIRS1Act, kIRS1Inact, kIRS1Phos, kPI3KDephos, kPI3KPhos, kPras40Dephos, kPras40Phos, kRAS, kS6KDephos, kS6KPhos, kTSC2Dephos, kTSC2Phos, kTSC2PhosBoost, kmTORC1Dephos, kmTORC1LysToCyt, kmTORC1Phos, kmTORC1_Dis, kmTORC1cytToLys, kmTORC2Dephos, kmTORC2DephosByS6K, kmTORC2_Dis, kmTOR_RIC_Comb, kmTOR_RPT_Comb, mTOR_T, mySubsA0, mySubsA1, mySubsA2, mySubsB0, mySubsB1, mySubsB2, mySubsB3, mySubsB4, mySubsB5, mySubsB6, mySubsB7, mySubsC0, mySubsC1, mySubsC2, mySubsC3, mySubsC4, mySubsD0, mySubsD1, mySubsD10, mySubsD11, mySubsD12, mySubsD13, mySubsD14, mySubsD15, mySubsD16, mySubsD17, mySubsD18, mySubsD19, mySubsD2, mySubsD20, mySubsD21, mySubsD22, mySubsD23, mySubsD24, mySubsD25, mySubsD26, mySubsD27, mySubsD28, mySubsD3, mySubsD4, mySubsD5, mySubsD6, mySubsD7, mySubsD8, mySubsD9, mySubsE0, mySubsE1, mySubsE10, mySubsE11, mySubsE12, mySubsE13, mySubsE14, mySubsE15, mySubsE16, mySubsE17, mySubsE18, mySubsE19, mySubsE2, mySubsE20, mySubsE21, mySubsE22, mySubsE23, mySubsE24, mySubsE25, mySubsE26, mySubsE27, mySubsE28, mySubsE29, mySubsE3, mySubsE30, mySubsE31, mySubsE32, mySubsE33, mySubsE34, mySubsE35, mySubsE36, mySubsE37, mySubsE38, mySubsE39, mySubsE4, mySubsE40, mySubsE41, mySubsE42, mySubsE43, mySubsE44, mySubsE45, mySubsE46, mySubsE47, mySubsE48, mySubsE49, mySubsE5, mySubsE50, mySubsE51, mySubsE52, mySubsE53, mySubsE54, mySubsE55, mySubsE56, mySubsE57, mySubsE58, mySubsE59, mySubsE6, mySubsE60, mySubsE61, mySubsE62, mySubsE7, mySubsE8, mySubsE9, mySubsF0, mySubsF1, mySubsF10, mySubsF11, mySubsF12, mySubsF13, mySubsF14, mySubsF15, mySubsF16, mySubsF17, mySubsF18, mySubsF19, mySubsF2, mySubsF20, mySubsF21, mySubsF22, mySubsF23, mySubsF24, mySubsF25, mySubsF26, mySubsF27, mySubsF28, mySubsF29, mySubsF3, mySubsF30, mySubsF31, mySubsF32, mySubsF33, mySubsF34, mySubsF35, mySubsF36, mySubsF37, mySubsF38, mySubsF39, mySubsF4, mySubsF40, mySubsF41, mySubsF42, mySubsF43, mySubsF44, mySubsF45, mySubsF46, mySubsF47, mySubsF48, mySubsF49, mySubsF5, mySubsF50, mySubsF51, mySubsF52, mySubsF53, mySubsF54, mySubsF55, mySubsF56, mySubsF57, mySubsF58, mySubsF59, mySubsF6, mySubsF60, mySubsF61, mySubsF62, mySubsF63, mySubsF64, mySubsF65, mySubsF7, mySubsF8, mySubsF9, mySubsG0, mySubsG1, mySubsG2
    compartments : Cell, default_compartment
    assignments : PRAS40_sum, S6K_sum, TSC2_wb, TSC2pT1462_wb, AktpT308_wb, IRS1_wb, IRS1pS636_639_wb, AktpS473_wb, FourEBP1pT37_46_wb, FourEBP1_wb, PRAS40_wb, Akt_wb, S6KpT389_wb, PI3K_sum, FourEBP1_sum, S6K_wb, IRS1_sum, PRAS40pT246_wb, Akt_sum, TSC2_sum

a more compleat human readable json of the data is saved in metaData.json