  function [minFlux,maxFlux,B,Brange,Trange,Frange,TMY,dynam1,currentP,stoplabel,rxnList]=...
    BTconstraintSearch(model,targetMet,GUR,OUR,minGrowth,glucoseRxnID,oxygenRxnID,biomassRxnID,Pvalue,rxnNum)
%BTconstraintSearch is a function of DynCubeProd.
%BTconstraintSearch obtains the constraints for the target production, 
%biomass objective function and sum of absolute values of total fluxes, 
%and returns the value of achieved target production.
%
%INPUTS
% model    COBRA model structure
% targetMet   target metabolite 
% GUR     glucose uptake ratio
% OUR     oxygen uptake ratio
% minGrowth    The minimum value of biomass objective fucnction that
%              DynCubeProd must satisfy.
% glucoseRxnID     ID of glucose reaction in the COBRA model
% oxygenRxnID      ID of oxygen reaction in the COBRA model
% biomassRxnID     ID of biomass reaction in the COBRA model
% Pvalue           Maximum cube size
% rxnNum           The raction ID of target metabolite exchange reaction
% 
%OUTPUTS
% minFlux    minimum flux of target reaction abtained by FVA
% maxFlux    maximum flux of target reaction abtained by FVA
% B    Biomass rate after applying candidate unused reactions
% Brange    The position of found knockout strategy on GR axis
% Trange    The position of found knockout strategy on PR axis
% Frange    The position of found knockout strategy on SF axis
% TMY    The theoretical maximum yield of target metabolite
% dynam1    A matrix with PXPXP size, in which the value of elements represent 
%         whether it is feasible to find a solution in this sub-cube. 
%         0 represents it is infeasible otherwise it is feasible when the 
%         element is 1.
% currentP    The value of P when candidate knockout strategy is found
% stoplabel   1 or 0. 1 means funciton ends before all cubes are computed,
%             0 means all cubes are computed.
% rxnList    List of found unused reactions
%
%  May 12, 2021   Tamura Takeyuki, Ma Yier
%

model.lb(glucoseRxnID)=-GUR;
model.lb(oxygenRxnID)=-OUR;
%model.ub(oxygenRxnID)=-OUR;
target=findMetIDs(model,targetMet);

if rxnNum==0
    [model2,rxnIDexists]=addReaction(model,'Transport',{targetMet},[-1]);
    m=size(model2.mets,1);
    n=size(model2.rxns,1);
    model2.S(target,n)=-1;
    model2.ub(n)=999999;
    model2.lb(n)=0;
    model2.rev(n)=0;
    targetRID=n;
else
    targetRID=rxnNum;
    model2=model;
    m=size(model2.mets,1);
    n=size(model2.rxns,1);
end

opt2=optimizeCbModel(model2);
MB=opt2.f;
model3=model2;
model3.c(biomassRxnID)=0;
model3.c(targetRID)=1;
model3.lb(biomassRxnID)=minGrowth;
opt2=optimizeCbModel(model3);
TMY=opt2.f;
TMF=2*sum(abs(opt2.v));

model3=model2;
model3.lb(biomassRxnID)=minGrowth;


lp.f=[zeros(n,1); ones(n,1)];
lp.A=[-eye(n) -eye(n);eye(n) -eye(n);zeros(1,n) -ones(1,n);zeros(1,n) ones(1,n)];
lp.b=zeros(2*n,1);
lp.Aeq=[model3.S zeros(m,n)];
lp.beq=zeros(m,1);
lp.lb=[model3.lb; zeros(n,1)];
lp.ub=[model3.ub; 999999*ones(n,1)];


expnum=log2(Pvalue);
dymli=ones(1,(expnum+1));
for i=1:expnum
    dymli(1,(i+1))=pow2(i);
end
dynam2=ones(1,1,1);
currentP=0;

for P=dymli
    if P ~= 1  
        if P==2
            s=[1,1,1];
        else
            s=size(dynam1);
        end
        dynam2 = zeros(P,P,P);
        for i=1:s(1)
            for j=1:s(2)
                for k=1:s(3)
                    dynam2((2*i-1),(2*j-1),(2*k-1))=dynam1(i,j,k);
                    dynam2((2*i-1),(2*j-1),2*k)=dynam1(i,j,k);
                    dynam2((2*i-1),2*j,(2*k-1))=dynam1(i,j,k);
                    dynam2((2*i-1),2*j,2*k)=dynam1(i,j,k);
                    dynam2(2*i,(2*j-1),(2*k-1))=dynam1(i,j,k);
                    dynam2(2*i,(2*j-1),2*k)=dynam1(i,j,k);
                    dynam2(2*i,2*j,(2*k-1))=dynam1(i,j,k);
                    dynam2(2*i,2*j,2*k)=dynam1(i,j,k);
                end
            end
        end
    end


    dynam1 = zeros(P,P,P);
    minf=zeros(P,P,P);
    maxf=zeros(P,P,P);
    opt5biomass = zeros(P,P,P);
    stoplabel=0;

    for i=1:P
        biomassLB=(MB/P)*(i-1);
        biomassUB=(MB/P)*i;
        for j=1:P
            targetLB=(TMY/P)*(j-1);
            targetUB=(TMY/P)*j;
            for k=1:P
                TMFlb=(TMF/P)*(k-1);
                TMFub=(TMF/P)*k;
                if dynam2(i,j,k)==1
                    [minf(i,j,k),maxf(i,j,k),opt5biomass(i,j,k),dynam1(i,j,k),rxnList]=...
                        integrate(model2,model3,targetRID,biomassLB,biomassUB,targetLB,targetUB,n,biomassRxnID,lp,TMFlb,TMFub,minGrowth);

                    if minf(i,j,k)>=TMY*0.1&&opt5biomass(i,j,k)>=minGrowth
                        stoplabel=1;
                        Brange=i;
                        Trange=j;
                        Frange=k;
                        currentP=P;
                        minFlux=minf(i,j,k);
                        maxFlux=maxf(i,j,k);
                        B=opt5biomass(i,j,k);

                    end
                    
                end
                if stoplabel==1
                    break;
                end
            end
            if stoplabel==1
                break;
            end
        end
        if stoplabel==1
            break;
        end
    end
    if stoplabel==1
        break;
    end
end


if stoplabel==0   
    table2=minf;
    table2(opt5biomass<minGrowth)=0;
    [M,I]=max(table2(:));
    [Brange,Trange,Frange] = ind2sub(size(table2),I);
    B=opt5biomass(Brange,Trange,Frange);
    minFlux=minf(Brange,Trange,Frange);
    maxFlux=maxf(Brange,Trange,Frange);
    currentP=Pvalue; 
    rxnList=0;
end

%filename=sprintf('results/BTconditionSearch_%d.mat',target);
%save(filename);

end

