function[rxnNum]=CMet(model,glucoseRxnID,oxygenRxnID,targetMet,GUR,OUR)
%CMet is function to calculate the reaction ID of metabolite exchange
%reaction.
%
%INPUTS
% model    COBRA model structure
% glucoseRxnID    Reaction representing glucose uptake (e.g., EX_glc__D_e)
% oxygenRxnID    Reaction representing oxygen uptake (e.g., EX_o2_e)
% targetMet    target metabolites
% GUR    Glucose uptake ratio
% OUR    Oxygen uptake ratio
%
%OUTPUTS
% rxnNum reaction ID of target metabolite exchange reaction
%
% May 12, 2021    Ma Yier
%

modelc=model;
modelc.lb(glucoseRxnID)=-GUR;
modelc.lb(oxygenRxnID)=-OUR;
exName=join(['EX_',targetMet]);
dmName=join(['DM_',targetMet]);
rxnNum=0;

for k=1:size(modelc.rxns)
    if strcmp(modelc.rxns{k,1},exName)||strcmp(modelc.rxns{k,1},dmName)
        rxnNum=k;
        break;
    end
end

