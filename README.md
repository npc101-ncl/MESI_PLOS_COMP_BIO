# using this tool
This repository contains a set of command line tools used to run parameter estimations using copasi as a back end. The 2 important feature it provides are the ability to force copasi to allow a model to reach a stable state before running the simulations used for parameter estimation. The second feature is the ability to run multiple copies of a model concurrently within the parameterisation process and apply a normalisation step to the data produced by the simulations comprising the parameterisation process.

It stands alongside our research paper and can be used to recreate to process used to generate the results in that paper.

The simplest way to run the toy parameter estimation is to move into the directory with GPParamiteriser.py and toyScriptSub.json in and run the command line:

`python GPParamiteriser.py configFile:toyScriptSub.json copasiPath:<insert path to copasi>`

GPParamiteriser has the following package dependencies: pandas, json, sys, os, pickle, time, re, site, pycotools3, tellurium, logging, random, subprocess, math, gc

You will want to make sure these are installed first. It also only works under python 3.6 so you’ll need to make sure that is your current python environment.

Running this command will run a parameter estimation using the configuration options in the toyScriptSub.json file. This file can itself point to other files. toyScriptSub.json also points to the data file toy_data.xlsx, the model files antStringToy.txt, antStringToyM.txt and antStringToyMP.txt.

The toyScriptSub.json is actually a simplification of a prior system where you had to pass several configuration files and numerous other command line arguments. But this prior system is still supported and used for a lot of the work in this repository.

The base model file is antStringToy.txt that describes the model we actually want to parameterise. antStringToyM.txt is the model created by the stringWangler.py tool using the stabToy.json configuration file. antStringToyM.txt describes a model that has the built in, variable length, stabilisation period. antStringToyMP.txt is the model created from antStringToyM.txt by modelParaliser.py using the toyParaliser.json configuration file. The antStringToyMP.txt model adds multiple sub models for each cell line and adds the data normalisation step at the end of the model simulation.

You can find full and thorough documentation for the tools in the readmes: readme GPParamiteriser.rtf, readme GPVisualiser.rtf, readme modelParaliser.rtf, readme stringWangler.rtf.

This toy example is simple enough that it can be run on a single core locally. For practical use this system uses slurm to run the code on multiple HPC cores.
You can submit the shell file GPParamiteriserToy.sh to slurm using sbatch:

`sbatch GPParamiteriserToy.sh`

You will generally need to have copasi set up for use in your slurm environment. On a local machine the code will attempt to manually add copasi to the path before calling it but on a slurm system it will assume it is already on the path.

There is also an example of how to call GPParamiteriser.py to perform a profile likelihoods analysis. You can run this using in slurm using the shell file GPParamiteriserToyPL.sh which utilises the configuration file toyScriptSubPL.json.

`sbatch GPParamiteriserToyPL.sh`
## replicating figuers
The figures in the paper 5-8 were created with data using scripts 
*GPParamiteriserMCF7_P.sh
*GPParamiteriserPaper.sh
*GPParamiteriserPaperRed.sh
*GPParamiteriserZR75_P.sh
Respectively. Expect these scripts to use around 200 cores and take up to 2 days to run.
