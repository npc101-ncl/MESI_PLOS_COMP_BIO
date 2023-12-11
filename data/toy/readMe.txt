this file is auto generated from the comandRecord.p and elements.p files.



filenames and identifyers:

toy : folder name used to identify output data and figures

toy : identifier used to distinguish output sets in the same folder and allow parallel execution of different parameter estimations

parameters.csv : filename of output parameter estimation

PETimeCourses.csv : time course simulations for the first few 'best' parameter sets in the parameter estimation

DBTimeCourses.csv : time course simulations for the first few 'best' parameter sets in the parameter estimation done using the same model used for parameter estimation. used for debugging

toy_data.xlsx : file name of experimental data source

toyScriptSub.json : a JSON file used as an alternative and means to override command line arguments passed to the parameter estimation

elements.p : file containing data about the named elements composing the base model and derived models used for parameter estimation

antStringToy.txt : filename for base model we are parameterising

antStringToyM.txt : filename for model with pre stabilisation period

antStringToyMP.txt : filename for model with multiple parallel submodes

toyScriptSub.json : file used to link names in experimental data to names in model

toyScriptSub.json : file expressing relationships between variables in experimental data. Which measurement is a contributor to which

toyScriptSub.json : file expressing which names in experimental data relate to total protein levels.

parameter estimation settings:

100 : number of parameter sets in estimated output

200.0 : time allowed for model to stabilise

100.0 : upper bound for estimated parameters

0.01 : lower bound for estimated parameters

the method and parameters used in the first phase of parameter estimations
    method : particle_swarm
    swarm_size : 50
    iteration_limit : 2000

the method and parameters used in the second phase of parameter estimations
    method : hooke_jeeves

elments (parameters, variables etc) in the models:

model elements in the base model
    metabolites : A, Ap, B, Bp
    kineticParamsGlobal : Stim, k1, k2, k3, k4
    compartments : default_compartment
    assignments : Ap_wb, Bp_wb, B_wb, A_wb

model elements in the model with stabilisation period
    metabolites : A, Ap, B, Bp, PT
    kineticParamsGlobal : A_T, A_wb_DVT0, A_wb_DVT1, A_wb_DVT2, A_wb_DVT3, A_wb_DVT4, A_wb_DVT5, Ap_wb_DVT0, Ap_wb_DVT1, Ap_wb_DVT2, Ap_wb_DVT3, Ap_wb_DVT4, Ap_wb_DVT5, B_T, B_wb_DVT0, B_wb_DVT1, B_wb_DVT2, B_wb_DVT3, B_wb_DVT4, B_wb_DVT5, Bp_wb_DVT0, Bp_wb_DVT1, Bp_wb_DVT2, Bp_wb_DVT3, Bp_wb_DVT4, Bp_wb_DVT5, PTSpeed, Stim, k1, k2, k3, k4, mySubsA0, mySubsA1, mySubsA2
    compartments : default_compartment
    assignments : Ap_wb, Bp_wb, B_wb, A_wb

model elements in the model with multiple paralel sub models
    metabolites : PT_cellLineA, PT_cellLineB, cellLineA_A, cellLineA_Ap, cellLineA_B, cellLineA_Bp, cellLineB_A, cellLineB_Ap, cellLineB_B, cellLineB_Bp
    kineticParamsGlobal : A_T_cellLineA, A_T_cellLineB, A_wb_DVT0__cellLineA, A_wb_DVT0__cellLineB, A_wb_DVT0_cellLineA, A_wb_DVT0_cellLineB, A_wb_DVT1__cellLineA, A_wb_DVT1__cellLineB, A_wb_DVT1_cellLineA, A_wb_DVT1_cellLineB, A_wb_DVT2__cellLineA, A_wb_DVT2__cellLineB, A_wb_DVT2_cellLineA, A_wb_DVT2_cellLineB, A_wb_DVT3__cellLineA, A_wb_DVT3__cellLineB, A_wb_DVT3_cellLineA, A_wb_DVT3_cellLineB, A_wb_DVT4__cellLineA, A_wb_DVT4__cellLineB, A_wb_DVT4_cellLineA, A_wb_DVT4_cellLineB, A_wb_DVT5__cellLineA, A_wb_DVT5__cellLineB, A_wb_DVT5_cellLineA, A_wb_DVT5_cellLineB, A_wb_mean, Ap_wb_DVT0__cellLineA, Ap_wb_DVT0__cellLineB, Ap_wb_DVT0_cellLineA, Ap_wb_DVT0_cellLineB, Ap_wb_DVT1__cellLineA, Ap_wb_DVT1__cellLineB, Ap_wb_DVT1_cellLineA, Ap_wb_DVT1_cellLineB, Ap_wb_DVT2__cellLineA, Ap_wb_DVT2__cellLineB, Ap_wb_DVT2_cellLineA, Ap_wb_DVT2_cellLineB, Ap_wb_DVT3__cellLineA, Ap_wb_DVT3__cellLineB, Ap_wb_DVT3_cellLineA, Ap_wb_DVT3_cellLineB, Ap_wb_DVT4__cellLineA, Ap_wb_DVT4__cellLineB, Ap_wb_DVT4_cellLineA, Ap_wb_DVT4_cellLineB, Ap_wb_DVT5__cellLineA, Ap_wb_DVT5__cellLineB, Ap_wb_DVT5_cellLineA, Ap_wb_DVT5_cellLineB, Ap_wb_mean, B_T_cellLineA, B_T_cellLineB, B_wb_DVT0__cellLineA, B_wb_DVT0__cellLineB, B_wb_DVT0_cellLineA, B_wb_DVT0_cellLineB, B_wb_DVT1__cellLineA, B_wb_DVT1__cellLineB, B_wb_DVT1_cellLineA, B_wb_DVT1_cellLineB, B_wb_DVT2__cellLineA, B_wb_DVT2__cellLineB, B_wb_DVT2_cellLineA, B_wb_DVT2_cellLineB, B_wb_DVT3__cellLineA, B_wb_DVT3__cellLineB, B_wb_DVT3_cellLineA, B_wb_DVT3_cellLineB, B_wb_DVT4__cellLineA, B_wb_DVT4__cellLineB, B_wb_DVT4_cellLineA, B_wb_DVT4_cellLineB, B_wb_DVT5__cellLineA, B_wb_DVT5__cellLineB, B_wb_DVT5_cellLineA, B_wb_DVT5_cellLineB, B_wb_mean, Bp_wb_DVT0__cellLineA, Bp_wb_DVT0__cellLineB, Bp_wb_DVT0_cellLineA, Bp_wb_DVT0_cellLineB, Bp_wb_DVT1__cellLineA, Bp_wb_DVT1__cellLineB, Bp_wb_DVT1_cellLineA, Bp_wb_DVT1_cellLineB, Bp_wb_DVT2__cellLineA, Bp_wb_DVT2__cellLineB, Bp_wb_DVT2_cellLineA, Bp_wb_DVT2_cellLineB, Bp_wb_DVT3__cellLineA, Bp_wb_DVT3__cellLineB, Bp_wb_DVT3_cellLineA, Bp_wb_DVT3_cellLineB, Bp_wb_DVT4__cellLineA, Bp_wb_DVT4__cellLineB, Bp_wb_DVT4_cellLineA, Bp_wb_DVT4_cellLineB, Bp_wb_DVT5__cellLineA, Bp_wb_DVT5__cellLineB, Bp_wb_DVT5_cellLineA, Bp_wb_DVT5_cellLineB, Bp_wb_mean, cellLineA_PTSpeed, cellLineA_Stim, cellLineA_mySubsA0, cellLineA_mySubsA1, cellLineA_mySubsA2, cellLineB_PTSpeed, cellLineB_Stim, cellLineB_mySubsA0, cellLineB_mySubsA1, cellLineB_mySubsA2, k1, k2, k3, k4
    compartments : default_compartment
    assignments : cellLineA_A_wb, cellLineA_Bp_wb, cellLineA_Ap_wb, cellLineB_Bp_wb, cellLineA_B_wb, cellLineB_Ap_wb, cellLineB_A_wb, cellLineB_B_wb

a more compleat human readable json of the data is saved in metaData.json