function [] = test()
%This is the test script of DynCubeProd by running an example of obtaining 
%the results of '10-Formyltetrahydrofolate' as target metabolites with 
%input glucose=10, oxygen=10 and P=2.
%
% May 12, 2021     Ma Yier
%

initCobraToolbox;
load("iJO1366.mat");
model=iJO1366;
[minf,maxf,B,dynam,b,t,f,TMY,currentP,stoplabel,rxnList]=...
DynCubeProd(model,model.mets(1,1),"EX_glc__D_e","EX_o2_e","BIOMASS_Ec_iJO1366_core_53p95M","GUR",10,"OUR",10,"minGrowth",0.01,"P",2)

filename=sprintf("results/rxnList_%s.mat",model.mets{1,1});
save(filename,"rxnList");
end