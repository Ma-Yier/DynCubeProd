function [minf,maxf,opt5biomass,dynaml,rxnList] =...
    integrate(model2,model3,targetRID,biomassLB,biomassUB,targetLB,targetUB,n,biomassRxnID,lp,TMFlb,TMFub,minGrowth)
%'integrate' is a function of BTconstraintSearch.
%'integrate' employs linear programming twice in the designated sub-cube.
%In the first LP, the unused reaction candidates are obtained.
%In the second LP, the production rate is calculated when reactions obtained in the first LP are not used.
%
%INPUTS
% model2   the model that a transport reaction is added to the original
%          model.
% model3   the model that the objective function has been modified in model2.
% targetRID  ID of the transport reaction of the target metabolite in
%            model2 and model3.
% biomassLB  the lower bound for the value of the biomass objective
%            function in the first LP.
% biomassUB  the upper bound for the value of the biomass objective
%            function in the first LP.
% targetLB   the lower bound for the value of the target metabolite
%            production in the first LP.
% targetUB   the upper bound for the value of the target metabolite
%            production in the first LP.
% n          the number of reactions
% biomassRxnID  ID of the reaction of biomass objectice function in model2
%               and model3
% lp         an linear programming formalization
% TMFlb       lower bound of sum of absolute value of total fluxes for the
%             first LP.
% TMFub       upper bound of sum of absolute value of total fluxes for the
%             first LP.
% minGrowth    minimum rate of biomass reaction
%
%OUTPUTS
% minf    minimum flux of target reaction abtained by FVA
% maxf    maximum flux of target reaction abtained by FVA
% opt5biomass  the value of biomass objective function in the second LP
% dynam1    0 represents it is infeasible while 1 represents it is
%           feasible.
% rxnList    List of found unused reactions
%
% May 12, 2021   Tamura Takeyuki, Ma Yier
%

model4=model3;
lp.lb(biomassRxnID)=max(minGrowth,biomassLB);
lp.ub(biomassRxnID)=biomassUB;
lp.lb(targetRID)=targetLB;
lp.ub(targetRID)=targetUB;
lp.b1=[lp.b;-TMFlb;TMFub];

[opt4.x, opt4.f, opt4.stat] = gurobiLinProg(lp.f, lp.A, lp.b1, lp.Aeq, lp.beq, lp.lb, lp.ub);

if opt4.stat==1
    dynaml=1;
else
    dynaml=0;
end

if opt4.stat~=1
    opt5biomass=-99999;
    minf=0;
    maxf=0;
    return
end


usedRxns=find(abs(opt4.x)>=0.0000001);
blockedRxns=setdiff(([1:n])',usedRxns);
rxnList=model4.rxns(blockedRxns);
model5=changeRxnBounds(model2,rxnList,0,'b');


%FVA
opt=optimizeCbModel(model5);
opt5biomass=opt.f;
model6=model5;
model6.c(biomassRxnID)=0;
model6.c(targetRID)=1;
model6.lb(biomassRxnID)=opt5biomass;
model6.ub(biomassRxnID)=opt5biomass;
optmin=optimizeCbModel(model6,'min');
optmax=optimizeCbModel(model6,'max');
minf=optmin.f;
maxf=optmax.f;

return
end

