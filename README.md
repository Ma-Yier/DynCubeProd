## About DynCubeProd
#### This is the DynCubeProd cource code. DynCubeProd is a tool to computing knockout strategies for genome-scale metabolic models by dynamic bilevel linear programming.

## Necessary Environment
#### An environment with MATLAB, COBRA Toolbox is required. And a solver assigned to [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/installation.html#solver-compatibility) is recommended (e.g. Gurobi).

## Example Code
#### To run the example, type the following in MATLAB command line:
#### >>test()
#### Then an example results of *'10-Formyltetrahydrofolate'* as target metabolites with input glucose=10, oxygen=10 and P=2 is in the *./results* file.

## How to use DynCubeProd
#### First init COBRA toolbox and load the data set via following command:
#### >>initCobraToolbox
#### >>load("iJO1366.mat") %example data set%
#### The example data set of *E.coli* iJO1366 is in the released DynCubeProd1.0. It can also be downloaded from [BiGG database](http://bigg.ucsd.edu/models/iJO1366).
#### After loading the data set, run the following code:
#### >>[minf,maxf,B,dynam,b,t,f,TMY,currentP,stoplabel,rxnList]=...
#### DynCubeProd(model,mets,glc_ex_reaction,o2_ex_reaction,biomass_reacion,"GUR",glc_uptake,"OUR",o2_uptake,"minGrowth",minGrowth,"P",P)
#### model is the COBRA structure model, mets is the target metabolites from model.mets. glc_ex_reaction and o2_ex_reaction are glucose exchange reaction and oxygen exchange reacion from model.rxns, respectively. biomass_reacion is the biomass growth reaction from model.rxns as well. The option variables are the following vars. "GUR", "OUR", "minGrowth" and "P" are strings to represent the var behind themselves. glc_uptake, o2_uptake and minGrowth are glucose uptake rate, oxygen uptake rate and minimum growth rate defined by user. P is the value of P which defines how many cubes should be constructed.
#### As an example, the following code for test:
#### >>[minf,maxf,B,dynam,b,t,f,TMY,currentP,stoplabel,rxnList]=...
#### DynCubeProd(model,model.mets(1,1),"EX_glc__D_e","EX_o2_e","BIOMASS_Ec_iJO1366_core_53p95M","GUR",10,"OUR",10,"minGrowth",0.01,"P",2)
#### This is to obtain the results of *'10-Formyltetrahydrofolate'* as target metabolites with input glucose=10, oxygen=10 and P=2.
