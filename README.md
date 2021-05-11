## About DynCubeProd
### This is the DynCubeProd cource code. DynCubeProd is a tool to computing knockout strategies for genome-scale metabolic models by dynamic bilevel linear programming.
## Necessary Environment
### An environment with MATLAB, COBRA Toolbox is required. And a solver assigned to COBRA Toolbox is recommended (e.g. Gurobi).
## Example Code
### The example data set of E.coli iJO1366 is in the released DynCubeProd1.0. You can also download iJO1366 data set from [BiGG database](http://bigg.ucsd.edu/models/iJO1366).
### Run $ \[ minf,maxf,B,dynam,b,t,f,TMY,currentP,stoplabel,rxnList \] =...DynCubeProd(model,model.mets(1,1),'EX_glc__D_e','EX_o2_e','BIOMASS_Ec_iJO1366_core_53p95M','GUR',10,'OUR',10,'minGrowth',0.01,'P',2)$ to obtain the results of $ \' 10-Formyltetrahydrofolate \' $ as target metabolites with input glucose 10, oxygen 10 and P=2.
